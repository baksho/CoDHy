"""
llm_judge.py — LLM-as-a-Judge evaluation module for CoDHy

Purpose
-------
Independently scores the FINAL hypotheses produced by the "Full" pipeline
variant using a panel of external LLM judges, on a FIXED rubric/prompt
template. This is deliberately separate from hgenerator.py / hvalidator.py:
those files build a dynamically-changing prompt per (gene, disease, round);
this module builds ONE static rubric template per hypothesis and only fills
in the hypothesis's own fields (combination, mechanism, reasoning, etc.) —
the instructions/rubric text itself never changes across judges or items.

Judge panel / baseline panel (rationale — see accompanying chat writeup for
full history of how this evolved):
  - meta-llama/Llama-3.1-8B-Instruct is EXCLUDED as a judge/baseline on
    purpose: it is the model CoDHy itself (DrugScreeningAgent/ValidationAgent)
    uses to GENERATE hypotheses. Using a model to judge or compete against
    its own architecture's outputs is a documented self-preference bias risk.
  - BASELINE_MODELS (5): the zero-shot baseline generators. Each is prompted
    independently with no KG, embeddings, consistency filtering, or
    feedback/taboo context. All run locally via Ollama (see substitution
    notes below).
  - DEFAULT_JUDGE_PANEL (6): the same 5 BASELINE_MODELS, plus DeepSeek-R1 as
    a JUDGE-ONLY model that never generates baseline hypotheses of its own.
    Because 5 of the 6 judges are also baseline generators, a judge is
    excluded from scoring the specific hypothesis it generated itself (see
    _judge_is_self) — this is the self-preference-bias control. DeepSeek-R1
    never matches any hypothesis's source_model_id (it never generates), so
    it always scores every hypothesis, including every baseline's — and,
    like the other 5, it also always scores CoDHy's hypotheses, since CoDHy
    isn't any judge's own output either.
"""
import os
import re
import json
import time
import requests
from huggingface_hub import InferenceClient


# --------------------------------------------------------------------------- #
#  Model panels
# --------------------------------------------------------------------------- #
# All models run locally via Ollama after confirming none of the originally-
# specified checkpoints are reliably reachable via HF Inference Providers
# (exhausted PRO credits for some, "not a chat model"/"not supported by any
# provider" for others — verified from actual run errors, not assumed).
#
# Every entry below is a SUBSTITUTE for an originally-specified checkpoint
# except OpenBioLLM-8B (exact same checkpoint, GGUF-quantized only). State
# each substitution explicitly in the paper's methods section:
#
#   - Qwen3-8B ("qwen3:8b"): replaces Qwen3-235B-A22B-Instruct-2507 (and
#     replaces the intermediate qwen3:32b substitute used briefly before
#     this). Confirmed official Ollama library tag (5.2GB, dense). Qwen3-235B
#     is NOT used anywhere in this project — as either judge or baseline.
#   - JSL-MedLlama-3-8B-v9 ("hf.co/mradermacher/JSL-MedLlama-3-8B-v9-GGUF:Q4_K_M"):
#     v9, not the originally-specified v2.0 — same John Snow Labs lineage,
#     different release. Verified working via `ollama pull`.
#   - OpenBioLLM-8B ("hf.co/MoMonir/Llama3-OpenBioLLM-8B-GGUF:Q4_K_M"): SAME
#     checkpoint as originally specified (aaditya/Llama3-OpenBioLLM-8B),
#     just GGUF-quantized — no version substitution here.
#   - BioMistral-7B ("hf.co/BioMistral/BioMistral-7B-GGUF:Q4_K_M"): this is
#     the BASE BioMistral-7B model, not the -SLERP merge variant originally
#     specified — same lineage (published by the BioMistral org itself),
#     different specific checkpoint.
#   - MedGemma-1.5-4B ("medgemma1.5:4b"): newly added per the updated
#     evaluation plan, pulled locally — not a substitution for anything
#     previously specified, this is a new addition to the baseline set.
#   - DeepSeek-R1-8B-Distill ("deepseek-r1:8b"): a Qwen/Llama model
#     DISTILLED on DeepSeek-R1's reasoning traces, NOT the actual 671B
#     DeepSeek-R1 MoE model (the 671b tag is the real model but won't fit a
#     typical Colab session). Judge-only — see module docstring.
BASELINE_MODELS = [
    {"name": "Qwen3-8B",              "model_id": "qwen3:8b",                                               "backend": "ollama"},
    {"name": "MedGemma-1.5-4B",       "model_id": "medgemma1.5:4b",                                          "backend": "ollama"},
    {"name": "JSL-MedLlama-3-8B-v9",  "model_id": "hf.co/mradermacher/JSL-MedLlama-3-8B-v9-GGUF:Q4_K_M",     "backend": "ollama"},
    {"name": "OpenBioLLM-8B",         "model_id": "hf.co/MoMonir/Llama3-OpenBioLLM-8B-GGUF:Q4_K_M",          "backend": "ollama"},
    {"name": "BioMistral-7B",         "model_id": "hf.co/BioMistral/BioMistral-7B-GGUF:Q4_K_M",              "backend": "ollama"},
]

_DEEPSEEK_JUDGE_ONLY = {"name": "DeepSeek-R1-8B-Distill", "model_id": "deepseek-r1:8b", "backend": "ollama"}

# Judges = the 5 baseline generators (self-excluding on their own output) +
# DeepSeek-R1 as a 6th judge-only model. NOT reusing BASELINE_MODELS by
# reference — a plain list copy, so appending the judge-only model here
# never accidentally makes DeepSeek-R1 a baseline generator too.
DEFAULT_JUDGE_PANEL = BASELINE_MODELS + [_DEEPSEEK_JUDGE_ONLY]

# Scoring dimensions — 6 criteria, 1-5 scale. Clinical relevance was dropped
# per the updated evaluation plan ("model judgment cannot establish clinical
# utility"); "mechanistic_explanation" renamed to "mechanistic_clarity" to
# match the paper's column header. If you revert either of these, update
# build_judge_prompt/build_pairwise_prompt's text below to match — the field
# list and the prompt text are not auto-synced.
JUDGE_SCORE_FIELDS = [
    "novelty",
    "biological_plausibility",
    "evidence_grounding",
    "mechanistic_clarity",
    "specificity",
    "testability",
]
JUDGE_SCORE_MIN = 1
JUDGE_SCORE_MAX = 5

JUDGE_SYSTEM_PROMPT = (
    "You are an independent, skeptical biomedical reviewer evaluating an "
    "AI-generated drug-combination hypothesis for a research paper. You are "
    "NOT told which system produced this hypothesis, and you did not "
    "generate it yourself — evaluate it strictly and honestly, exactly as a "
    "critical, blinded peer reviewer would."
)


def build_judge_prompt(hypothesis: dict) -> str:
    """
    Builds the FIXED rubric prompt for one hypothesis. This template string
    is identical for every judge, every hypothesis, and every generating
    system — only the hypothesis's own content fields are substituted. This
    is what keeps judging BLIND: nothing here ever reveals which system
    (CoDHy or a baseline) produced the hypothesis, so callers must never
    pass system_name/source_model_id into this prompt.
    """
    return f"""
[HYPOTHESIS UNDER REVIEW]
Combination: {hypothesis.get('combination', 'N/A')}
Target Disease: {hypothesis.get('target_disease', 'N/A')}
Focus Gene: {hypothesis.get('focus_gene', 'N/A')}
Proposed Mechanism: {hypothesis.get('mechanism', 'N/A')}
Stated Reasoning: {hypothesis.get('reasoning', 'N/A')}
Rationale Type Claimed: {hypothesis.get('rationale_type', 'N/A')}
Supporting Evidence Cited: {hypothesis.get('supporting_evidence', 'N/A')}

[TASK]
Score this hypothesis on each dimension below using an INTEGER from 1 to 5
(1 = very poor, 5 = excellent):

1. novelty — does this go beyond the single most obvious/well-known combination for this gene/disease?
2. biological_plausibility — is the underlying biology of this combination sound?
3. evidence_grounding — is the cited evidence real, relevant, and specific (not vague or fabricated)?
4. mechanistic_clarity — does the proposed mechanism logically hold together end-to-end?
5. specificity — is this a concrete, well-defined proposal (vs. a vague/generic statement)?
6. testability — could this realistically be turned into a specific, falsifiable experiment or trial?

Then write a 2-3 sentence justification explaining your scores, and explicitly
flag anything that looks fabricated, vague, or scientifically implausible.

HARD RULES:
- Do not reward the hypothesis for confident-sounding language alone.
- If supporting evidence looks vague or generic ("general biological inference"
  with no real citation), evidence_grounding should be low.
- Base your judgment only on what is written above plus your own biomedical
  knowledge — do not assume missing details in the hypothesis's favor.

OUTPUT FORMAT (strict JSON only, no other text, no markdown fences):
{{
  "novelty": 0,
  "biological_plausibility": 0,
  "evidence_grounding": 0,
  "mechanistic_clarity": 0,
  "specificity": 0,
  "testability": 0,
  "justification": "..."
}}
"""


def _extract_json_dict(text: str) -> dict:
    """Same robust-extraction pattern used elsewhere in this project
    (hvalidator.py / feedback_agent.py) for consistency."""
    if not text:
        return {}
    if '</think>' in text.lower():
        idx = text.lower().find('</think>') + len('</think>')
        text = text[idx:]
    text = text.strip()
    text = re.sub(r'```json', '', text)
    text = re.sub(r'```', '', text).strip()
    try:
        data = json.loads(text)
        if isinstance(data, dict):
            return data
    except Exception:
        pass
    try:
        match = re.search(r'\{.*\}', text, re.DOTALL)
        if match:
            return json.loads(match.group())
    except Exception:
        pass
    return {}


def _empty_score_dict(reason: str) -> dict:
    d = {field: None for field in JUDGE_SCORE_FIELDS}
    d["justification"] = f"(judge call failed: {reason})"
    d["_parse_failed"] = True
    return d


# JSON schemas for Ollama's `format` parameter (structured-output, constrained
# decoding — documented at https://ollama.com/blog/structured-outputs). This
# is the SAME fix that was applied to hgenerator.py's 404/markdown-output
# problem earlier, but it was never applied to judge scoring calls until now
# — that gap is why weaker/"thinking" local models (MedGemma, DeepSeek-R1
# distill) were producing preamble text, commentary, or truncated JSON
# instead of valid score objects.
JUDGE_SCORE_JSON_SCHEMA = {
    "type": "object",
    "properties": {
        **{f: {"type": "integer", "minimum": JUDGE_SCORE_MIN, "maximum": JUDGE_SCORE_MAX} for f in JUDGE_SCORE_FIELDS},
        "justification": {"type": "string"},
    },
    "required": JUDGE_SCORE_FIELDS + ["justification"],
}

PAIRWISE_JSON_SCHEMA = {
    "type": "object",
    "properties": {
        "winner": {"type": "string", "enum": ["A", "B", "Tie"]},
        "justification": {"type": "string"},
    },
    "required": ["winner", "justification"],
}


def _call_hf_judge(model_id: str, system_prompt: str, user_prompt: str,
                    max_tokens: int = 1000, temperature: float = 0.0,
                    hf_token: str | None = None) -> str:
    """Calls a judge model via HF Inference Providers (router)."""
    client = InferenceClient(model_id, token=hf_token or os.environ.get("HF_TOKEN"))
    response = client.chat_completion(
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        max_tokens=max_tokens,
        temperature=temperature,
    )
    return response.choices[0].message.content


def _call_ollama_judge(model_id: str, system_prompt: str, user_prompt: str,
                        ollama_base_url: str = "http://localhost:11434/api/generate",
                        max_tokens: int = 2000, temperature: float = 0.0,
                        response_schema: dict = None) -> str:
    """Calls a judge model via a local Ollama server, matching the pattern
    used in hgenerator.py/hvalidator.py/feedback_agent.py. response_schema,
    when given, is passed as Ollama's `format` parameter for constrained
    decoding — pass JUDGE_SCORE_JSON_SCHEMA or PAIRWISE_JSON_SCHEMA rather
    than relying on prompt instructions alone. max_tokens default raised
    from 1000 to 2000: reasoning-style local models (DeepSeek-R1 distill,
    MedGemma) spend a meaningful chunk of their budget on preamble/thinking
    tokens BEFORE the constrained JSON portion even starts, so schema
    constraint alone doesn't guarantee the response fits in a small budget."""
    payload = {
        "model": model_id,
        "prompt": system_prompt + " Output ONLY the JSON object, nothing else — no preamble, no commentary, no markdown fences.\n\n" + user_prompt,
        "stream": False,
        "options": {"temperature": temperature, "num_predict": max_tokens},
    }
    if response_schema:
        payload["format"] = response_schema
    resp = requests.post(ollama_base_url, json=payload, timeout=300)
    resp.raise_for_status()
    return resp.json().get("response", "")


def judge_single_hypothesis(hypothesis: dict, judge_cfg: dict,
                             ollama_base_url: str = "http://localhost:11434/api/generate",
                             hf_token: str | None = None,
                             max_retries: int = 2) -> dict:
    """
    Runs ONE judge on ONE hypothesis. Returns a dict with the judge's scores
    plus bookkeeping fields (judge name/model, hypothesis id).
    """
    user_prompt = build_judge_prompt(hypothesis)
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            if judge_cfg["backend"] == "hf":
                raw = _call_hf_judge(
                    judge_cfg["model_id"], JUDGE_SYSTEM_PROMPT, user_prompt, hf_token=hf_token
                )
            elif judge_cfg["backend"] == "ollama":
                raw = _call_ollama_judge(
                    judge_cfg["model_id"], JUDGE_SYSTEM_PROMPT, user_prompt,
                    ollama_base_url=ollama_base_url, response_schema=JUDGE_SCORE_JSON_SCHEMA,
                )
            else:
                raise ValueError(f"Unknown backend '{judge_cfg['backend']}' for judge {judge_cfg['name']}")

            parsed = _extract_json_dict(raw)
            if parsed and all(f in parsed for f in JUDGE_SCORE_FIELDS):
                result = {f: parsed.get(f) for f in JUDGE_SCORE_FIELDS}
                result["justification"] = parsed.get("justification", "")
                result["_parse_failed"] = False
                break
            else:
                last_error = f"unparseable/incomplete JSON (attempt {attempt+1}): {raw[:200]}"
                result = _empty_score_dict(last_error)
        except Exception as e:
            last_error = str(e)
            result = _empty_score_dict(last_error)
            time.sleep(1.5 * (attempt + 1))

    if result["_parse_failed"]:
        # Previously silent — judge failures (bad HF router responses, model
        # not chat-completion-served, credits exhausted, etc.) never printed
        # anything, unlike baseline_generators.py which does. That asymmetry
        # made it easy to assume judging succeeded when it silently didn't.
        print(f"  [FAILED] {judge_cfg['name']} on hypothesis {hypothesis.get('id')}: {last_error}")

    result.update({
        "judge_name": judge_cfg["name"],
        "judge_model_id": judge_cfg["model_id"],
        "hypothesis_id": hypothesis.get("id"),
        "combination": hypothesis.get("combination", ""),
        "focus_gene": hypothesis.get("focus_gene", ""),
        "target_disease": hypothesis.get("target_disease", ""),
        # Bookkeeping only — NEVER passed into build_judge_prompt(), so
        # judging stays blind. Used downstream for self-exclusion and for
        # grouping scores by system in judge_metrics/objective_metrics.
        "system_name": hypothesis.get("system_name", "CoDHy-Full"),
        "source_model_id": hypothesis.get("source_model_id"),
    })
    return result


def _judge_is_self(judge_cfg: dict, hypothesis: dict) -> bool:
    """
    True if this judge would be scoring a hypothesis generated by (the same
    checkpoint as) itself — e.g. a DeepSeek-R1 zero-shot baseline hypothesis
    being judged by the DeepSeek-R1 judge. Compares on source_model_id,
    which baseline_generators.py stamps onto every hypothesis it produces
    using the SAME model_id strings as DEFAULT_JUDGE_PANEL, so this is an
    exact string match, not a heuristic.
    """
    return bool(hypothesis.get("source_model_id")) and hypothesis.get("source_model_id") == judge_cfg["model_id"]


def run_judge_panel(hypotheses: list[dict], judge_panel: list[dict] = None,
                     only_method: str | None = "Full",
                     n_repeats: int = 1,
                     ollama_base_url: str = "http://localhost:11434/api/generate",
                     hf_token: str | None = None,
                     verbose: bool = True) -> "pandas.DataFrame":
    """
    Runs the full judge panel over a list of hypothesis dicts and returns a
    long-format DataFrame: one row per (hypothesis, judge, repeat), with one
    column per score field.

    only_method: if set (default "Full"), hypotheses whose 'method' field
    doesn't match are skipped before judging. Set to None when judging a
    mix of CoDHy + baseline hypotheses that don't share a 'method' value.

    n_repeats: run each (hypothesis, judge) pair this many times (supervisor
    point 5 — repeat the evaluation for robustness). Rows are tagged with
    'repeat_idx' so judge_metrics.repeat_consistency() can report test-retest
    agreement separately from inter-judge agreement.

    Self-exclusion (bias control): if a hypothesis's source_model_id matches
    a judge's model_id (i.e. that judge would be scoring its own baseline
    output), that judge is SKIPPED for that hypothesis and the skip is
    logged, rather than silently producing a biased score.
    """
    import pandas as pd

    judge_panel = judge_panel or DEFAULT_JUDGE_PANEL
    active_panel = [j for j in judge_panel if j["model_id"] != "TODO_FIND_GGUF_OR_DROP"]
    skipped = [j["name"] for j in judge_panel if j not in active_panel]
    if skipped and verbose:
        print(f"Skipping unconfigured judges (still TODO placeholders): {skipped}")
    judge_panel = active_panel

    if only_method:
        filtered = [h for h in hypotheses if h.get("method") == only_method]
        if verbose:
            print(f"Filtered {len(hypotheses)} -> {len(filtered)} hypotheses with method == '{only_method}'.")
    else:
        filtered = hypotheses

    rows = []
    total = len(filtered) * len(judge_panel) * n_repeats
    done = 0
    for hypo in filtered:
        for judge_cfg in judge_panel:
            if _judge_is_self(judge_cfg, hypo):
                if verbose:
                    print(f"  [skip] {judge_cfg['name']} excluded from judging hypothesis "
                          f"{hypo.get('id')} — same model generated it (source_model_id="
                          f"{hypo.get('source_model_id')}).")
                continue
            for rep in range(n_repeats):
                if verbose:
                    done += 1
                    print(f"  [{done}/{total}] {judge_cfg['name']} judging hypothesis {hypo.get('id')} "
                          f"({hypo.get('combination', '')[:60]}) [repeat {rep+1}/{n_repeats}]...")
                row = judge_single_hypothesis(
                    hypo, judge_cfg, ollama_base_url=ollama_base_url, hf_token=hf_token
                )
                row["repeat_idx"] = rep
                rows.append(row)

    df = pd.DataFrame(rows)
    return df


# --------------------------------------------------------------------------- #
#  Pairwise A/B comparison (supervisor point 4)
# --------------------------------------------------------------------------- #
PAIRWISE_SYSTEM_PROMPT = (
    "You are an independent, skeptical biomedical reviewer comparing two "
    "candidate drug-combination hypotheses for the SAME gene and disease. "
    "You are not told which system produced either one. Pick the stronger "
    "hypothesis, or say it's a genuine tie — do not default to 'tie' to "
    "avoid a decision; only use it when the two are truly comparable."
)


def build_pairwise_prompt(hyp_a: dict, hyp_b: dict) -> str:
    """
    Fixed comparison template. hyp_a/hyp_b are shown ONLY as "Hypothesis A"
    and "Hypothesis B" — no system_name or source_model_id ever appears
    here, keeping the comparison blind. Caller is responsible for
    randomizing which underlying system is "A" vs "B" per call (see
    run_pairwise_panel) to avoid position bias.
    """
    def block(h):
        return (
            f"Combination: {h.get('combination', 'N/A')}\n"
            f"Target Disease: {h.get('target_disease', 'N/A')}\n"
            f"Focus Gene: {h.get('focus_gene', 'N/A')}\n"
            f"Proposed Mechanism: {h.get('mechanism', 'N/A')}\n"
            f"Stated Reasoning: {h.get('reasoning', 'N/A')}\n"
            f"Supporting Evidence Cited: {h.get('supporting_evidence', 'N/A')}"
        )

    return f"""
[HYPOTHESIS A]
{block(hyp_a)}

[HYPOTHESIS B]
{block(hyp_b)}

[TASK]
Both hypotheses address the same gene and disease. Considering biological
plausibility, evidence grounding, mechanistic clarity, specificity, and
testability together, which is the stronger hypothesis?

OUTPUT FORMAT (strict JSON only, no other text, no markdown fences):
{{
  "winner": "A" | "B" | "Tie",
  "justification": "..."
}}
"""


def judge_pairwise(hyp_a: dict, hyp_b: dict, judge_cfg: dict,
                    ollama_base_url: str = "http://localhost:11434/api/generate",
                    hf_token: str | None = None, max_retries: int = 2) -> dict:
    """Runs ONE judge on ONE (hyp_a, hyp_b) pair. Returns winner + bookkeeping."""
    user_prompt = build_pairwise_prompt(hyp_a, hyp_b)
    result = {"winner": None, "justification": "", "_parse_failed": True}
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            if judge_cfg["backend"] == "hf":
                raw = _call_hf_judge(judge_cfg["model_id"], PAIRWISE_SYSTEM_PROMPT, user_prompt,
                                      max_tokens=400, hf_token=hf_token)
            elif judge_cfg["backend"] == "ollama":
                raw = _call_ollama_judge(judge_cfg["model_id"], PAIRWISE_SYSTEM_PROMPT, user_prompt,
                                          ollama_base_url=ollama_base_url, max_tokens=800,
                                          response_schema=PAIRWISE_JSON_SCHEMA)
            else:
                raise ValueError(f"Unknown backend '{judge_cfg['backend']}'")

            parsed = _extract_json_dict(raw)
            if parsed.get("winner") in ("A", "B", "Tie"):
                result = {"winner": parsed["winner"], "justification": parsed.get("justification", ""),
                          "_parse_failed": False}
                break
            last_error = f"unparseable/invalid winner (attempt {attempt+1}): {raw[:200]}"
        except Exception as e:
            last_error = str(e)
            time.sleep(1.5 * (attempt + 1))

    if result["_parse_failed"]:
        result["justification"] = f"(pairwise judge call failed: {last_error})"
    return result


def run_pairwise_panel(hyps_system_a: list[dict], hyps_system_b: list[dict],
                        system_a_name: str, system_b_name: str,
                        judge_panel: list[dict] = None,
                        ollama_base_url: str = "http://localhost:11434/api/generate",
                        hf_token: str | None = None,
                        random_seed: int = 0,
                        verbose: bool = True) -> "pandas.DataFrame":
    """
    Pairs hyps_system_a[i] with hyps_system_b[i] (same gene/disease context
    and same generation index — per supervisor point 2, both systems should
    already have been given the same biomarker/disease/hypo_count, so index
    alignment is meaningful). For each pair and each judge (excluding any
    judge that generated either hypothesis — see _judge_is_self), the A/B
    position is randomized per call to control for position bias, then the
    win is attributed back to the correct underlying system regardless of
    which slot it was shown in.
    """
    import pandas as pd
    import random

    rng = random.Random(random_seed)
    judge_panel = judge_panel or DEFAULT_JUDGE_PANEL
    judge_panel = [j for j in judge_panel if j["model_id"] != "TODO_FIND_GGUF_OR_DROP"]
    n = min(len(hyps_system_a), len(hyps_system_b))
    if len(hyps_system_a) != len(hyps_system_b) and verbose:
        print(f"Warning: system lengths differ ({len(hyps_system_a)} vs {len(hyps_system_b)}); "
              f"comparing only the first {n} matched pairs.")

    rows = []
    for i in range(n):
        h_a, h_b = hyps_system_a[i], hyps_system_b[i]
        for judge_cfg in judge_panel:
            if _judge_is_self(judge_cfg, h_a) or _judge_is_self(judge_cfg, h_b):
                if verbose:
                    print(f"  [skip] {judge_cfg['name']} excluded from pair {i} — generated one side of it.")
                continue

            # Randomize which underlying system is shown as "A".
            if rng.random() < 0.5:
                shown_a, shown_b = h_a, h_b
                a_is_system_a = True
            else:
                shown_a, shown_b = h_b, h_a
                a_is_system_a = False

            if verbose:
                print(f"  {judge_cfg['name']} comparing pair {i}: {system_a_name} vs {system_b_name}...")
            res = judge_pairwise(shown_a, shown_b, judge_cfg, ollama_base_url=ollama_base_url, hf_token=hf_token)

            if res["winner"] == "Tie" or res["winner"] is None:
                winning_system = "Tie" if res["winner"] == "Tie" else None
            elif (res["winner"] == "A") == a_is_system_a:
                winning_system = system_a_name
            else:
                winning_system = system_b_name

            rows.append({
                "pair_index": i,
                "judge_name": judge_cfg["name"],
                "system_a": system_a_name,
                "system_b": system_b_name,
                "focus_gene": h_a.get("focus_gene", ""),
                "target_disease": h_a.get("target_disease", ""),
                "raw_winner_slot": res["winner"],
                "winning_system": winning_system,
                "justification": res["justification"],
                "_parse_failed": res["_parse_failed"],
            })

    return pd.DataFrame(rows)

"""
baseline_generators.py — Zero-shot baseline hypothesis generation for CoDHy

Purpose (supervisor's evaluation plan, points 1 & 2)
-----------------------------------------------------
Produces drug-combination hypotheses for the SAME (gene, disease, hypo_count)
CoDHy was run on, using a plain zero-shot prompt — no knowledge graph, no
embeddings, no VanillaConsistencyChecker, no taboo/feedback context. This is
what "a general-purpose LLM / biomedical LLM baseline" concretely means: the
same underlying models you're already using as judges, but here asked to
generate directly from their own pretrained knowledge instead of scoring.

Every hypothesis produced here is stamped with:
  - system_name        e.g. "DeepSeek-R1-ZeroShot"  (human-readable, for tables)
  - source_model_id     e.g. "deepseek-ai/DeepSeek-R1"  (exact judge-panel
                         model_id string, used by llm_judge._judge_is_self()
                         to stop that same model from judging its own
                         baseline output — see llm_judge.py docstring)
  - method              "Baseline-ZeroShot"  (distinct from CoDHy's own
                         "Full"/"No-Embed"/"LLM-only" method tags, so you can
                         filter cleanly)

Output schema matches CoDHy's own hypothesis dicts (id, method,
target_disease, focus_gene, reasoning, combination, mechanism,
rationale_type, supporting_evidence) so it flows straight into
llm_judge.run_judge_panel() / run_pairwise_panel() and
objective_metrics.py unchanged.

Extending to non-panel baselines (e.g. GPT-4o, Gemini):
  Add an entry to BASELINE_MODELS with a new backend value (e.g. "openai",
  "gemini") and a matching branch in _call_baseline_model(). Not implemented
  here since it needs API keys/SDKs I haven't confirmed you have — this
  module is structured so adding one is a small, isolated change.
"""
import os
import re
import json
import time
from huggingface_hub import InferenceClient

# BASELINE_MODELS (5 models) is the source of truth for baseline generation.
# It does NOT include DeepSeek-R1 — DeepSeek-R1 is judge-only (see
# llm_judge.py's DEFAULT_JUDGE_PANEL / _DEEPSEEK_JUDGE_ONLY). Importing the
# narrower BASELINE_MODELS list here (rather than DEFAULT_JUDGE_PANEL, which
# now also includes DeepSeek-R1) is what keeps that separation enforced —
# this module can never accidentally generate a DeepSeek-R1 baseline.
from llm_judge import BASELINE_MODELS, _call_hf_judge, _call_ollama_judge

BASELINE_SYSTEM_PROMPT = (
    "You are an expert AI biologist proposing candidate drug-combination "
    "hypotheses. Rely ONLY on your own pretrained biomedical knowledge — "
    "you have NOT been given any external database, literature search, or "
    "knowledge graph for this task."
)


def _build_baseline_prompt(cancer_type: str, focus_gene: str, hypo_count: int) -> str:
    """
    Fixed zero-shot prompt template, deliberately parallel in structure to
    hgenerator.py's OUTPUT FORMAT block so downstream code (parsing, judge
    prompts, objective metrics) can treat CoDHy and baseline hypotheses
    identically. No KG/embedding context is given, per supervisor point 1.
    """
    return f"""
Cancer: {cancer_type}
Biomarker: {focus_gene}

TASK: Propose exactly {hypo_count} distinct drug-combination hypotheses that
could suppress {cancer_type} progression targeting {focus_gene}. Use your own
biomedical knowledge only.

INSTRUCTIONS:
1. Generate EXACTLY {hypo_count} hypotheses.
2. Each hypothesis MUST use the focus gene + a different secondary drug/target.
3. Briefly explain the biological rationale for each (mechanism).
4. Supporting evidence: describe the general basis for the claim in your own
   words. DO NOT invent PubMed IDs, DOIs, or URLs of any kind — if you don't
   have a specific citation, write "General biological inference" instead.

OUTPUT FORMAT (strict JSON list of {hypo_count} items, no markdown fences,
no text before or after the list):
[
  {{
    "id": 1,
    "method": "Baseline-ZeroShot",
    "target_disease": "{cancer_type}",
    "focus_gene": "{focus_gene}",
    "reasoning": "Step-by-step biological reasoning for this combination.",
    "combination": "Drug A + Drug B",
    "mechanism": "...",
    "rationale_type": "LLM Internal Knowledge",
    "supporting_evidence": "Explanation... (no invented citations)"
  }}
]
"""


def _baseline_json_schema_for(hypo_count: int) -> dict:
    """
    Builds the array schema WITH minItems (and maxItems) set to the exact
    requested count. The previous module-level BASELINE_JSON_SCHEMA had no
    minItems, so an empty array [] technically satisfied it — that's exactly
    the bug behind BioMistral-7B repeatedly returning "[ ]" as a "valid"
    baseline response. min/max both set to hypo_count enforces the count the
    prompt already asks for, not just "at least one item".
    """
    return {
        "type": "array",
        "minItems": hypo_count,
        "maxItems": hypo_count,
        "items": {
            "type": "object",
            "properties": {
                "id": {"type": "integer"},
                "method": {"type": "string"},
                "target_disease": {"type": "string"},
                "focus_gene": {"type": "string"},
                "reasoning": {"type": "string"},
                "combination": {"type": "string"},
                "mechanism": {"type": "string"},
                "rationale_type": {"type": "string"},
                "supporting_evidence": {"type": "string"},
            },
            "required": [
                "id", "method", "target_disease", "focus_gene", "reasoning",
                "combination", "mechanism", "rationale_type", "supporting_evidence",
            ],
        },
    }


def _extract_json_list(text: str) -> list:
    """Same robust-extraction pattern used across this project."""
    if not text:
        return []
    text = re.sub(r'```json', '', text)
    text = re.sub(r'```', '', text).strip()
    try:
        data = json.loads(text)
        if isinstance(data, list):
            return [x for x in data if isinstance(x, dict)]
    except Exception:
        pass
    try:
        match = re.search(r'\[.*\]', text, re.DOTALL)
        if match:
            data = json.loads(match.group())
            if isinstance(data, list):
                return [x for x in data if isinstance(x, dict)]
    except Exception:
        pass
    return []


def generate_baseline_hypotheses(model_cfg: dict, cancer_type: str, focus_gene: str,
                                  hypo_count: int = 4,
                                  ollama_base_url: str = "http://localhost:11434/api/generate",
                                  hf_token: str | None = None,
                                  max_retries: int = 2) -> list[dict]:
    """
    Generates `hypo_count` zero-shot baseline hypotheses from ONE model.
    Returns [] (with a printed warning) if the model fails to produce
    parseable JSON after retries — callers should check for this before
    handing an empty list to the judge/pairwise panels.
    """
    user_prompt = _build_baseline_prompt(cancer_type, focus_gene, hypo_count)
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            if model_cfg["backend"] == "hf":
                raw = _call_hf_judge(model_cfg["model_id"], BASELINE_SYSTEM_PROMPT, user_prompt,
                                      max_tokens=4000, temperature=0.2, hf_token=hf_token)
            elif model_cfg["backend"] == "ollama":
                # NOTE: /api/generate's `format` (JSON-schema constrained decoding)
                # only helps for Ollama-backed models — the "hf" branch above
                # relies on prompt instructions + retries instead, since HF's
                # router chat_completion() doesn't expose the same parameter
                # for every provider. Verify this against your Ollama version
                # if you rely on it for a baseline model.
                raw = _call_ollama_judge_with_format(
                    model_cfg["model_id"], BASELINE_SYSTEM_PROMPT, user_prompt,
                    ollama_base_url=ollama_base_url, max_tokens=4000, temperature=0.2,
                    response_schema=_baseline_json_schema_for(hypo_count),
                )
            else:
                raise ValueError(f"Unknown backend '{model_cfg['backend']}'")

            hyps = _extract_json_list(raw)
            if hyps:
                for h in hyps:
                    h["system_name"] = f"{model_cfg['name']}-ZeroShot"
                    h["source_model_id"] = model_cfg["model_id"]
                    h.setdefault("method", "Baseline-ZeroShot")
                return hyps
            last_error = f"unparseable/empty JSON (attempt {attempt+1}): {raw[:200]}"
        except Exception as e:
            last_error = str(e)
            time.sleep(1.5 * (attempt + 1))

    print(f"  Baseline generation FAILED for {model_cfg['name']} on {focus_gene}/{cancer_type}: {last_error}")
    return []


def _call_ollama_judge_with_format(model_id, system_prompt, user_prompt, ollama_base_url,
                                    max_tokens=4000, temperature=0.2, response_schema=None):
    """Ollama call with structured-output `format` schema for baseline generation
    (see the fix applied to hgenerator.py for the same reasoning/citation)."""
    import requests
    payload = {
        "model": model_id,
        "prompt": system_prompt + " Output ONLY valid JSON, nothing else.\n\n" + user_prompt,
        "stream": False,
        "format": response_schema,
        "options": {"temperature": temperature, "num_predict": max_tokens},
    }
    resp = requests.post(ollama_base_url, json=payload, timeout=300)
    resp.raise_for_status()
    return resp.json().get("response", "")


def generate_all_baselines(cancer_type: str, focus_gene: str, hypo_count: int = 4,
                            models: list[dict] = None,
                            ollama_base_url: str = "http://localhost:11434/api/generate",
                            hf_token: str | None = None,
                            verbose: bool = True) -> dict:
    """
    Runs generate_baseline_hypotheses() for every model in `models` (default:
    all 5 judge-panel models). Returns {system_name: [hypothesis, ...]}, one
    entry per model, skipping (with a printed note) any model that produced
    nothing after retries.
    """
    models = models or BASELINE_MODELS
    results = {}
    for model_cfg in models:
        if verbose:
            print(f"Generating baseline hypotheses: {model_cfg['name']} on {focus_gene}/{cancer_type}...")
        hyps = generate_baseline_hypotheses(
            model_cfg, cancer_type, focus_gene, hypo_count=hypo_count,
            ollama_base_url=ollama_base_url, hf_token=hf_token,
        )
        if hyps:
            results[f"{model_cfg['name']}-ZeroShot"] = hyps
        elif verbose:
            print(f"  -> skipped {model_cfg['name']} (no valid hypotheses returned).")
    return results

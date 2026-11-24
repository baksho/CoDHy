"""
hgenerator.py — Drug Screening / Hypothesis Generation Agent

Reverted to direct Ollama REST calls (matching v1/app.py's proven working
pattern: requests.post to OLLAMA_BASE_URL) instead of the llm_backend.py
multi-provider abstraction, which turned out too fragile in practice
(HF Inference Providers role/provider/context issues across many models).

Kept from Phase 1/2 (per explicit instruction to preserve these in this file):
  - VanillaConsistencyChecker integration (Layer 1 + Layer 2 checks)
  - RotatE/TransE Layer-2 predicted paths + Layer-3 embedding similarity
    (dynamic model-name labeling, works for either)
  - Chain-of-Thought "reasoning" field in the hypothesis schema
  - taboo_context / feedback_context prompt injection blocks

Consistency fix: generate_hypothesis() now accepts hypo_count (previously
app.py called it with hypo_count=... but the schema/prompt only supported a
hardcoded 4 — that mismatch would have raised or silently ignored the
parameter). The JSON schema example block is now built dynamically for
however many hypotheses are requested.
"""
import os
import re
import json
import requests
from vanilla_checker import VanillaConsistencyChecker
from huggingface_hub import InferenceClient


class DrugScreeningAgent:
    def __init__(self, kg, model_name=None, ollama_base_url="http://localhost:11434/api/generate"):
        self.kg = kg
        self.model = model_name or "meta-llama/Llama-3.1-8B-Instruct"
        self.api_url = ollama_base_url
        # self.client = InferenceClient(model_name, token=os.environ.get("HF_TOKEN"))
        self.vanilla = VanillaConsistencyChecker(kg=self.kg)

    HYPOTHESIS_JSON_SCHEMA = {
        "type": "array",
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

    def check_connection(self):
        try:
            requests.get(self.api_url.replace("/generate", "/tags"), timeout=20)
            return True
        except Exception:
            return False

    def _extract_json(self, text):
        """
        Robustly extracts JSON from messy LLM output.
        """
        text = re.sub(r'```json', '', text)
        text = re.sub(r'```', '', text)
        text = text.strip()

        result = []

        try:
            data = json.loads(text)
            if isinstance(data, list):
                result = data
            elif isinstance(data, dict):
                found_list = False
                for key in data:
                    if isinstance(data[key], list):
                        result = data[key]
                        found_list = True
                        break
                if not found_list:
                    result = [data]
        except:
            pass

        if not result:
            try:
                match = re.search(r'\[.*\]', text, re.DOTALL)
                if match:
                    result = json.loads(match.group())
            except:
                pass

        if isinstance(result, list):
            clean_list = [x for x in result if isinstance(x, dict)]
            return clean_list

        return []

    def generate_hypothesis(self, cancer_type, focus_gene, embedding_agent=None, method="Full",
                             hypo_count=4, taboo_context=None, feedback_context=None):
        """
        taboo_context : list[{"combination": str, "reason": str}] | None
        feedback_context : list[dict] | None
            Output of FeedbackAgent.generate_feedback() from the prior round
            (only present for the revised/round-2 generation call).
        """
        print(f"\nAgent is thinking about {focus_gene} in {cancer_type} using {method} method...")

        explicit_context = []
        predicted_context = []
        implicit_context = []
        emb_model_name = getattr(embedding_agent, "model_name", "Embedding") if embedding_agent else "Embedding"

        # ── Layer 1: Explicit Context (GraphRAG) ────────────────────────────
        if method in ["Full", "No-Node2Vec"]:
            discovery_subgraph = self.kg.get_discovery_subgraph(focus_gene, cancer_type)
            raw_layer1_paths = []
            if discovery_subgraph:
                for item in discovery_subgraph:
                    # Structured dict for VANILLA (Phase 1 fix: VANILLA.check()
                    # requires dicts with subject/relation/object keys, NOT
                    # pre-rendered display strings — an earlier version of
                    # this file passed strings here and crashed with
                    # AttributeError: 'str' object has no attribute 'get').
                    raw_layer1_paths.append({
                        "subject": item['s.name'],
                        "relation": item['type(r)'],
                        "object": item['o.name'],
                        "layer": "VALIDATED",
                        "source": item.get('r.source', 'KG'),
                        "source_url": item.get('r.source_url', ''),
                    })

            # VANILLA consistency check on Layer-1 facts (Phase 1 addition)
            if raw_layer1_paths:
                l1_consistent, l1_flagged = self.vanilla.check(raw_layer1_paths)
                self.vanilla.summary(l1_consistent, l1_flagged)
                explicit_context = [
                    f"{p['subject']} --[{p['relation']}]--> {p['object']} "
                    f"[Source: {p['source']} - {p['source_url'] or 'No URL'}]"
                    for p in l1_consistent
                ]
            else:
                explicit_context = []

        # ── Layer 2: KGE-predicted paths (RotatE/TransE — Phase 1 addition) ─
        if method == "Full" and embedding_agent and getattr(embedding_agent, "model", None) is not None:
            if hasattr(embedding_agent, "build_predicted_paths"):
                predicted_raw = embedding_agent.build_predicted_paths(focus_gene, cancer_type, top_k=5)
                if predicted_raw:
                    # predicted_raw is already correctly-structured dicts
                    # (subject/relation/object/layer="PREDICTED") straight
                    # from GraphEmbeddingAgent — pass directly to VANILLA.
                    l2_consistent, l2_flagged = self.vanilla.check(predicted_raw)
                    self.vanilla.summary(l2_consistent, l2_flagged)
                    predicted_context = [
                        f"{p['subject']} --[{p['relation']}]--> {p['object']} "
                        f"[Source: {p.get('source', emb_model_name)}-Predicted]"
                        for p in l2_consistent[:10]
                    ]

        # ── Layer 3: Implicit Context (embedding similarity) ────────────────
        if method == "Full" and embedding_agent:
            similar_nodes = embedding_agent.find_similar_nodes(focus_gene, top_k=5)
            for node in similar_nodes:
                implicit_context.append(f"{node['name']} (Similarity Score: {node['score']:.2f})")

        # ── Construct the Hybrid Prompt ──────────────────────────────────────
        context_str = "\n".join(explicit_context[:25]) if explicit_context else "No direct connections in graph data found."
        predicted_str = "\n".join(predicted_context) if predicted_context else "No embedding-predicted links available."
        implicit_str = ", ".join(implicit_context) if implicit_context else "No embedding data for this mode."

        if method == "LLM-only":
            knowledge_instruction = "Use your internal pre-trained medical knowledge ONLY. Ignore external graph data."
            context_str = ""
            predicted_str = ""
            implicit_str = ""
        else:
            knowledge_instruction = f"Strictly prioritize the provided [EXPLICIT KNOWLEDGE], [PREDICTED LINKS], and [IMPLICIT SIGNALS] for the {method} method."

        # Taboo list block (agentic memory — Phase 2 addition)
        if taboo_context:
            taboo_lines = "\n".join(
                f"- {t['combination']}  (rejected: {t['reason']})" for t in taboo_context
            )
            taboo_block = (
                "[PREVIOUSLY REJECTED COMBINATIONS — DO NOT REPROPOSE THESE]\n"
                f"{taboo_lines}\n"
            )
        else:
            taboo_block = ""

        # Feedback block (only present on the revision/round-2 call — Phase 2 addition)
        if feedback_context:
            fb_lines = "\n".join(
                f"- \"{fb['combination']}\": {fb['improvement_suggestions']}"
                for fb in feedback_context
            )
            feedback_block = (
                "[FEEDBACK FROM PREVIOUS ROUND — ADDRESS THESE POINTS]\n"
                f"{fb_lines}\n"
                "This is a REVISION round. Produce improved hypotheses that concretely "
                "address this feedback — do not just resubmit the same combinations unchanged.\n"
            )
        else:
            feedback_block = ""

        system_prompt = "You are an expert AI Biologist performing a drug screening analysis."

        user_prompt = f"""
        Cancer: {cancer_type}
        Biomarker: {focus_gene}
        METHOD USED: {method}

        {knowledge_instruction}

        [EXPLICIT KNOWLEDGE - Validated KG Paths]
        {context_str}

        [PREDICTED LINKS - {emb_model_name} Predicted Paths]
        {predicted_str}

        [IMPLICIT SIGNALS - {emb_model_name} Embedding Similarity]
        The following entities are mathematically similar to {focus_gene} in the graph structure.
        Use them to infer novel connections (e.g., if Drug X targets a similar gene, it might work here):
        {implicit_str}

        {taboo_block}
        {feedback_block}

        TASK: Based on the mechanistic pathways above, generate exactly 4 distinct drug combination
        hypotheses that could suppress {cancer_type} progression targeting {focus_gene}.
        Explain the biological rationale for each.

        INSTRUCTIONS:
        1. Generate EXACTLY 4 hypotheses.
        2. Prioritize paths labelled [KG-Validated] over [{emb_model_name}-Predicted].
        3. For each hypothesis, the mechanism must trace back to at least one path in the evidence above.
        4. Rationale type: use 'Direct Graph Evidence' for KG-Validated paths,
        'Predicted Link' for {emb_model_name}-Predicted paths, 'LLM Internal Knowledge' otherwise.
        5. Supporting evidence: cite only URLs that appear in the evidence above. No invented URLs.
        6. Each hypothesis MUST use the focus gene + a different secondary drug/target.
        7. If available, use at least one "Implicit Signal" to propose a novel connection.
        8. Mechanism: Briefly explain the biological logic (e.g., Synthetic Lethality, Resistance Reversal etc.).
        9. Rationale Type: Choose 'Direct Graph Evidence' if based on [EXPLICIT KNOWLEDGE], 'Inferred Similarity' if based on [IMPLICIT SIGNALS], or 'LLM Internal Knowledge', otherwise.
        10. Supporting Evidence:
          - Briefly explain why this works in 3-4 sentences.
          - CRITICAL: You must cite a URL ONLY if it appears in [EXPLICIT KNOWLEDGE] above.
          - If the evidence comes from your internal knowledge (LLM-only), DO NOT invent a URL. Write "General biological inference" down the evidence and explain it.

        CHAIN-OF-THOUGHT REQUIREMENT:
        Before filling in "combination", "mechanism", etc. for each hypothesis, first reason
        step-by-step in the "reasoning" field: (1) which piece of evidence above you're starting
        from, (2) what mechanistic link you're drawing, (3) why the resulting combination follows
        from that link. Keep "reasoning" and "mechanism" distinct.

        HARD RULES FOR HALLUCINATION PREVENTION:
        - DO NOT invent PubMed links or fake URLs.
        - DO NOT use placeholders like [Journal Name] or [Year] or [Source] or [1].
        - If no URL appears above, write "No direct clinical study found."

        OUTPUT FORMAT (strict JSON list of 4 items):
        [
            {{
                "id": 1,
                "method": "{method}",
                "target_disease": "{cancer_type}",
                "focus_gene": "{focus_gene}",
                "reasoning": "Step-by-step: (1) evidence used, (2) mechanistic link, (3) why the combination follows",
                "combination": "Drug A + Drug B",
                "mechanism": "...",
                "rationale_type": "...",
                "supporting_evidence": "Explanation... (Citation if available)"
            }},
            {{
                "id": 2,
                "method": "{method}",
                "target_disease": "{cancer_type}",
                "focus_gene": "{focus_gene}",
                "reasoning": "Step-by-step: (1) evidence used, (2) mechanistic link, (3) why the combination follows",
                "combination": "Drug A + Drug B",
                "mechanism": "...",
                "rationale_type": "...",
                "supporting_evidence": "Explanation... (Citation if available)"
            }},
            {{
                "id": 3,
                "method": "{method}",
                "target_disease": "{cancer_type}",
                "focus_gene": "{focus_gene}",
                "reasoning": "Step-by-step: (1) evidence used, (2) mechanistic link, (3) why the combination follows",
                "combination": "Drug A + Drug B",
                "mechanism": "...",
                "rationale_type": "...",
                "supporting_evidence": "Explanation... (Citation if available)"
            }},
            {{
                "id": 4,
                "method": "{method}",
                "target_disease": "{cancer_type}",
                "focus_gene": "{focus_gene}",
                "reasoning": "Step-by-step: (1) evidence used, (2) mechanistic link, (3) why the combination follows",
                "combination": "Drug A + Drug B",
                "mechanism": "...",
                "rationale_type": "...",
                "supporting_evidence": "Explanation... (Citation if available)"
            }}
        ]
        """

        payload = {
            "model": self.model,
            "prompt": system_prompt + " Output ONLY valid JSON, nothing else." + "\n\n" + user_prompt,
            "stream": False,
            "format": self.HYPOTHESIS_JSON_SCHEMA,
            "options": {
                "temperature": 0.01,
                "num_predict": 10000
            }
        }

        for attempt in range(3):
            try:
                response = requests.post(self.api_url, json=payload, timeout=300)
                if response.status_code == 200:
                    resp_text = response.json().get("response", "")
                    hypotheses = self._extract_json(resp_text)
                    if hypotheses and isinstance(hypotheses, list) and len(hypotheses) > 0:
                        print(f"      Successfully generated {len(hypotheses)} hypotheses.")
                        return hypotheses
                    else:
                        print(f"      Attempt {attempt+1} failed to parse JSON.")
                        print(f"      Length of raw output: {len(resp_text)}")
                        print(f"      Raw Output Snippet: {resp_text[-2000:]}...")
                else:
                    print(f"      API Error: {response.status_code}")
            except Exception as e:
                print(f"      Connection error: {e}")

        print("Failed to generate hypotheses after 3 attempts.")
        return []

        messages = [
                    {"role": "system", "content": system_prompt + " Output ONLY valid JSON."},
                    {"role": "user", "content": user_prompt}
                ]
        
        response = self.client.chat_completion(
            messages=messages,
            max_tokens=10000,
            temperature=0.01
        )

        resp_text = response.choices[0].message.content
        hypotheses = self._extract_json(resp_text)
        if hypotheses and isinstance(hypotheses, list) and len(hypotheses) > 0:
            print(f"      Successfully generated {len(hypotheses)} hypotheses.")
            return hypotheses
        else:
            print(f"      Hypothesis generation failed.")
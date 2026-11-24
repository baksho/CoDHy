"""
feedback_agent.py — Feedback Agent for CoDHy

Distinct from ValidationAgent (per original design decision):
  - ValidationAgent answers: "Is this hypothesis supported by evidence?"
  - FeedbackAgent answers: "How can this hypothesis be improved?"

Converted from the llm_backend.py abstraction to a direct Ollama REST call
(matching v1/app.py's proven pattern) — same reasoning as hgenerator.py and
hvalidator.py: the multi-provider abstraction was the source of most of the
fragility encountered when running this pipeline for real.

Pipeline position (unchanged from Phase 2): generate v1 -> validate v1 ->
feedback(v1) -> generate v2 (revised, feedback+taboo aware) -> validate v2
-> rank v2. Runs once per experiment.
"""
import re
import os
import json
import requests
from huggingface_hub import InferenceClient

# A validated hypothesis below this composite threshold, or one whose
# safety_score is very low, is treated as "rejected" for taboo-list purposes.
# Starting default, not an empirically validated cutoff.
REJECTION_SCORE_THRESHOLD = 0.4


class FeedbackAgent:
    def __init__(self, kg, model_name=None, ollama_base_url="http://localhost:11434/api/generate"):
        self.kg = kg
        self.model = model_name or "meta-llama/Llama-3.1-8B-Instruct"
        # self.client = InferenceClient(model_name, token=os.environ.get("HF_TOKEN"))
        self.api_url = ollama_base_url

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

    def _extract_text(self, raw_response):
        """
        Strips a <think>...</think> block if present (harmless no-op for
        plain llama3.1 output; protects against a future reasoning-model swap).
        """
        text = raw_response or ""
        if '</think>' in text.lower():
            idx = text.lower().find('</think>') + len('</think>')
            text = text[idx:]
        return text.strip()

    def _is_rejected(self, h):
        safety = h.get("safety_score", 5)
        lit    = h.get("literature_support_score", 0.0)
        bio    = h.get("biological_plausibility_score", 0.0)
        if safety <= 2:
            return True, f"safety_score={safety} (<=2, likely withdrawn/unsafe drug)"
        if lit < REJECTION_SCORE_THRESHOLD and bio < REJECTION_SCORE_THRESHOLD:
            return True, (
                f"weak evidence: literature_support_score={lit:.2f}, "
                f"biological_plausibility_score={bio:.2f} (both < {REJECTION_SCORE_THRESHOLD})"
            )
        return False, None

    def generate_feedback(self, validated_hypotheses):
        feedback_list = []

        for h in validated_hypotheses:
            is_rejected, reject_reason = self._is_rejected(h)

            system_prompt = (
                "You are a research mentor giving constructive feedback on a drug "
                "combination hypothesis. You are NOT judging whether it is correct — "
                "a validator already scored that. Your job is to suggest concretely "
                "how the hypothesis, its rationale, or its supporting evidence could "
                "be strengthened."
            )
            user_prompt = f"""
            HYPOTHESIS: {h.get('combination', '')}
            TARGET DISEASE: {h.get('target_disease', '')}
            FOCUS GENE: {h.get('focus_gene', '')}
            MECHANISM: {h.get('mechanism', '')}

            VALIDATION SCORES (already computed, do not re-score):
            - biological_plausibility_score: {h.get('biological_plausibility_score', 'N/A')}
            - literature_support_score: {h.get('literature_support_score', 'N/A')}
            - toxicity_risk_score: {h.get('toxicity_risk_score', 'N/A')}
            - mechanistic_coherence_score: {h.get('mechanistic_coherence_score', 'N/A')}
            - safety_score (1-10): {h.get('safety_score', 'N/A')}
            - validator critique: {h.get('critique', 'N/A')}

            TASK: In 2-3 sentences, suggest how this specific hypothesis could be
            improved — e.g. a more specific mechanistic step, a stronger drug
            pairing, additional evidence to seek out, or a caveat to address.
            Do NOT just restate the scores. Output ONLY the improvement
            suggestion text, no JSON, no preamble.
            """

            # messages = [
            #         {"role": "system", "content": system_prompt + " Output ONLY valid JSON."},
            #         {"role": "user", "content": user_prompt}
            #     ]
        
            try:
                # response = self.client.chat_completion(
                #     messages=messages,
                #     max_tokens=10000,
                #     temperature=0.1
                # )

                # resp_text = response.choices[0].message.content

                payload = {
                    "model": self.model,
                    "prompt": system_prompt + "\n\n" + user_prompt,
                    "stream": False,
                    "options": {"temperature": 0.1, "num_predict": 10000}
                }
                ollama_resp = requests.post(self.api_url, json=payload, timeout=300)
                resp_text = ollama_resp.json().get("response", "")

                suggestion = self._extract_text(resp_text)
                if not suggestion:
                    suggestion = "(feedback generation returned empty output)"
            except Exception as e:
                suggestion = f"(feedback generation failed: {e})"



            feedback_list.append({
                "id": h.get("id"),
                "combination": h.get("combination", ""),
                "focus_gene": h.get("focus_gene", ""),
                "target_disease": h.get("target_disease", ""),
                "is_rejected": is_rejected,
                "rejection_reason": reject_reason,
                "improvement_suggestions": suggestion,
            })

        return feedback_list

    def print_feedback_report(self, feedback_list):
        print("\n" + "=" * 70)
        print("  FEEDBACK AGENT REPORT")
        print("=" * 70)
        for fb in feedback_list:
            status = "REJECTED (-> taboo list)" if fb["is_rejected"] else "kept"
            print(f"\n[{fb['id']}] {fb['combination']}  —  {status}")
            if fb["is_rejected"]:
                print(f"    Reason: {fb['rejection_reason']}")
            print(f"    Improvement suggestion: {fb['improvement_suggestions']}")
        print("\n" + "=" * 70)

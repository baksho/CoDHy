"""
hvalidator.py — Validation Agent

Kept v1's direct Ollama REST call pattern (it was already using
requests.post to OLLAMA_BASE_URL — no llm_backend.py abstraction needed
here at all).

Changed vs. v1, per explicit instruction to preserve Phase 1/2 updates in
this file:
  - OUTPUT FORMAT restored to the richer 5-field scoring schema
    (biological_plausibility_score, literature_support_score,
    toxicity_risk_score, mechanistic_coherence_score, safety_score) instead
    of v1's 3-field schema (safety_score, plausibility, combination_toxicity_risk).
    This is required for FeedbackAgent's rejection logic, which reads these
    specific field names — v1's schema alone would leave FeedbackAgent's
    checks silently defaulting to 0.0/5 for every hypothesis.
  - JSON parsing made robust: v1's version was
        val_data = json.loads(r.json().get("response", "{}"))
    with a bare except that fell back to a fixed error dict on ANY parse
    failure — no bracket-search fallback, no markdown-fence stripping, no
    <think>-block handling. Replaced with the same robust _extract_json_dict
    logic developed in Phase 2 (harmless no-op for plain llama3.1 output,
    but protects against a future model swap to something that wraps its
    answer in commentary or a reasoning block).
"""
import os
import re
import json
import requests
from Bio import Entrez
from huggingface_hub import InferenceClient


class ValidationAgent:
    def __init__(self, kg, model_name=None, ollama_base_url="http://localhost:11434/api/generate",
                 email="you@example.com"):
        self.kg = kg
        self.model = model_name or "meta-llama/Llama-3.1-8B-Instruct"
        # self.client = InferenceClient(model_name, token=os.environ.get("HF_TOKEN"))
        self.api_url = ollama_base_url
        self.email = email
        Entrez.email = self.email

    VALIDATION_JSON_SCHEMA = {
        "type": "object",
        "properties": {
            "biological_plausibility_score": {"type": "number"},
            "literature_support_score": {"type": "number"},
            "toxicity_risk_score": {"type": "number"},
            "mechanistic_coherence_score": {"type": "number"},
            "safety_score": {"type": "number"},
            "critique": {"type": "string"},
            "supporting_evidence": {"type": "string"},
        },
        "required": [
            "biological_plausibility_score", "literature_support_score",
            "toxicity_risk_score", "mechanistic_coherence_score",
            "safety_score", "critique", "supporting_evidence",
        ],
    }

    @staticmethod
    def _extract_json_dict(text):
        """Robustly extracts a JSON dict from messy LLM output (see module docstring)."""
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

    def check_combination_evidence(self, drugs, cancer_type):
        """Performs a 'Just-in-Time' PubMed search for the specific combination."""
        if len(drugs) < 2:
            # Pre-existing bug fixed: this branch previously omitted
            # ne_hits/nc_hits/is_novel_exact/is_novel_context, which
            # validate_batch() unconditionally reads later — a malformed
            # combination string (parsing to <2 drugs) would KeyError here.
            return {
                "status": "Single Agent", "urls": [], "note": "Single drug",
                "ne_hits": 0, "nc_hits": 0,
                "is_novel_exact": True, "is_novel_context": True,
            }

        clean_drugs = [f'"{d}"' for d in drugs]
        query_ne = f"({' AND '.join(clean_drugs)})"
        query_nc = f"({' AND '.join(clean_drugs)}) AND {cancer_type}"

        ne_hits, ne_ids = self._get_pubmed_count(query_ne)
        nc_hits, nc_ids = self._get_pubmed_count(query_nc)

        found_urls = []
        if ne_hits > 0 or nc_hits > 0:
            all_ids = list(set(ne_ids + nc_ids))
            for pmid in all_ids:
                found_urls.append(f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")
            return {
                "status": "Validated (Literature Exists)", "urls": found_urls,
                "ne_hits": ne_hits, "nc_hits": nc_hits,
                "is_novel_exact": ne_hits == 0, "is_novel_context": nc_hits == 0,
            }

        return {
            "status": "Inferred (Novel Combination)", "urls": [],
            "ne_hits": ne_hits, "nc_hits": nc_hits,
            "is_novel_exact": ne_hits == 0, "is_novel_context": nc_hits == 0,
        }

    def _get_pubmed_count(self, query):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
            record = Entrez.read(handle)
            return int(record["Count"]), record["IdList"]
        except Exception:
            return 0, []

    def _sanitize_string(self, text):
        """Removes or escapes control characters that break JSON payloads."""
        if not isinstance(text, str):
            return str(text)
        # Replace actual newlines/tabs with escaped versions and remove non-printable chars
        text = text.replace("\n", "\\n").replace("\t", "\\t").replace("\r", "\\r")
        return "".join(char for char in text if char.isprintable() or char == "\\")

    def validate_batch(self, hypotheses_list):
        print("\nValidation Agent is reviewing candidates...")
        validated_results = []

        if not hypotheses_list:
            print("   No hypotheses to validate.")
            return []

        for hypo in hypotheses_list:
            if not isinstance(hypo, dict):
                print(f"   Skipping invalid data format: {hypo}")
                continue

            combo_str = hypo.get('combination', '')
            drugs = [d.strip().split()[0] for d in combo_str.replace("+", ",").split(",") if d.strip()]
            focus_gene = hypo.get('focus_gene')
            cancer_type = hypo.get('target_disease')

            print(f"   -> Reviewing Hypothesis {hypo.get('id', '?')}: {combo_str}")

            combo_check = self.check_combination_evidence(drugs, cancer_type)

            hypo['civic_urls'] = []
            hypo['chembl_urls'] = []
            hypo['trial_urls'] = []
            hypo['reactome_urls'] = []
            safety_context = []

            for drug in drugs:
                validated_facts = self.kg.get_validated_subgraph(
                    drug_name=drug, focus_gene=focus_gene, cancer_type=cancer_type
                ) if hasattr(self.kg, "get_validated_subgraph") else []
                fact_entries = []
                for record in validated_facts:
                    source = record.get('r.source', 'Unknown')
                    rel = record.get('type(r)', 'interacts with')
                    obj = record.get('o.name', 'Unknown')
                    fact_entries.append(f"({source} Fact): {drug} {rel} {obj}")

                query = """
                MATCH (d:Drug)
                WHERE toLower(d.name) CONTAINS toLower($name)
                  OR  toLower($name) CONTAINS toLower(d.name)

                OPTIONAL MATCH (g:Gene)
                WHERE toLower(g.name) CONTAINS toLower($focus_gene)
                  OR toLower($focus_gene) CONTAINS toLower(g.name)

                WITH d, g
                OPTIONAL MATCH (d)-[r1]-(:Phase)
                OPTIONAL MATCH (d)-[r2]-(g)

                WITH d, g, (collect(DISTINCT r1.source_url) + collect(DISTINCT r2.source_url)) as phase_urls
                OPTIONAL MATCH (d)-[r_civic]-(x)
                WHERE r_civic.source = "CIViC"

                WITH d, g, phase_urls, collect(DISTINCT r_civic.source_url) as civic_urls
                OPTIONAL MATCH (t:ClinicalTrial)-[r_trial1]-(d)
                OPTIONAL MATCH (t:ClinicalTrial)-[r_trial2]-(g)

                WITH d, g, phase_urls, civic_urls, (collect(DISTINCT r_trial1.source_url) + collect(DISTINCT r_trial2.source_url)) as trial_urls, collect(DISTINCT t.nct_id) as nct_ids
                OPTIONAL MATCH (g)-[r_pathway]-(p:Pathway)

                WITH d, g, phase_urls, civic_urls, trial_urls, nct_ids, collect(DISTINCT r_pathway.source_url) as reactome_urls
                OPTIONAL MATCH (d)-[:CAUSES_SIDE_EFFECT]->(se:SideEffect)

                RETURN
                    d.name as drug_name,
                    d.max_phase as phase,
                    d.is_withdrawn as withdrawn,
                    g.name as gene_name,
                    phase_urls,
                    civic_urls,
                    trial_urls,
                    reactome_urls,
                    collect(DISTINCT se.name)[..15] as side_effects
                """
                with self.kg.driver.session() as session:
                    results = session.run(query, name=drug, focus_gene=focus_gene).data()
                    if results:
                        res = results[0]
                        d_name = res['drug_name']
                        phase = res['phase']
                        withdrawn = res['withdrawn']
                        se_list = res['side_effects']

                        def clean_urls(url_list):
                            return [u for u in url_list if isinstance(u, str) and u.startswith("http")]

                        hypo['chembl_urls'].extend(clean_urls(res['phase_urls']))
                        hypo['civic_urls'].extend(clean_urls(res['civic_urls']))
                        hypo['trial_urls'].extend(clean_urls(res['trial_urls']))
                        hypo['reactome_urls'].extend(clean_urls(res['reactome_urls']))

                        status_str = f"Phase {phase}" if phase else "Unknown Status"
                        if withdrawn: status_str += " (WITHDRAWN)"
                        se_str = ", ".join(se_list) if se_list else "No data"
                        facts_str = f" | Biological Facts: {'; '.join(fact_entries)}" if fact_entries else ""
                        safety_context.append(f"DRUG: {d_name} | STATUS: {status_str} | SIDE EFFECTS: {se_str} | Biological Facts: {facts_str}")
                    else:
                        safety_context.append(f"DRUG: {drug} | Not found in Knowledge Graph.")

            hypo['civic_urls'] = list(set(hypo['civic_urls']))
            hypo['chembl_urls'] = list(set(hypo['chembl_urls']))
            hypo['trial_urls'] = list(set(hypo['trial_urls']))
            hypo['reactome_urls'] = list(set(hypo['reactome_urls']))

            safety_str = "\n".join(safety_context)

            system_prompt = "You are a clinical auditor. You need to validate the drug combination hypothesis for the given disease."
            user_prompt = f"""
            HYPOTHESIS: {combo_str}
            DISEASE: {cancer_type}

            [EVIDENCE CHECK]
            Status: {combo_check['status']}

            [INDIVIDUAL DRUG PROFILES]
            {safety_str}

            TASK:
            1. Score biological_plausibility_score (0.0-1.0): how well-supported is the underlying mechanism?
            2. Score literature_support_score (0.0-1.0): based on [EVIDENCE CHECK] above.
            3. Score toxicity_risk_score (0.0-1.0): combined toxicity risk, higher = riskier.
               If a drug is "WITHDRAWN", this MUST be >= 0.8.
            4. Score mechanistic_coherence_score (0.0-1.0): does the proposed mechanism make internal sense?
            5. Assign a Safety Score (1-10). If a drug is "WITHDRAWN", score must be < 3.
            6. Write a short critique.
            7. Provide supporting evidence:
              - For every claim, explain the finding.
              - For citation, you MUST output the "Real PubMed Hits" listed above in [EVIDENCE CHECK].
              - If the list is empty, WRITE "No direct clinical study found." and cite no URL.
              - DO NOT invent new URLs.

            HARD RULES FOR HALLUCINATION PREVENTION:
            - DO NOT use numeric placeholders like [1], [2].
            - DO NOT use placeholders like [Journal Name], [Year], or [Source] etc.
            - DO NOT generate fake PubMed links (e.g., pubmed.ncbi.nlm.nih.gov/12345678).
            - DO NOT make up Source Names or IDs.
            - If you do not see a "http..." link in the context provided above, DO NOT WRITE A URL.

            OUTPUT FORMAT (JSON only, no other text):
            {{
                "biological_plausibility_score": 0.0,
                "literature_support_score": 0.0,
                "toxicity_risk_score": 0.0,
                "mechanistic_coherence_score": 0.0,
                "safety_score": 8,
                "critique": "...",
                "supporting_evidence": "..."
            }}
            """
            hypo.update({
                'source_urls': list(dict.fromkeys(combo_check['urls'])),
                'evidence_status': combo_check['status'],
                'ne_hits': combo_check['ne_hits'],
                'nc_hits': combo_check['nc_hits'],
                'is_novel_exact': combo_check['is_novel_exact'],
                'is_novel_context': combo_check['is_novel_context']
            })

            user_prompt = self._sanitize_string(user_prompt)
            
            try:
                # messages = [
                #     {"role": "system", "content": system_prompt},
                #     {"role": "user", "content": user_prompt}
                # ]

                # response = self.client.chat_completion(
                #     messages=messages,
                #     max_tokens=4000,
                #     # response_format={"type": "json_object"}
                #     # temperature=0.01
                # )
                
                # raw_content = response.choices[0].message.content

                payload = {
                    "model": self.model,
                    "prompt": system_prompt + "\n\n" + user_prompt,
                    "stream": False,
                    "format": self.VALIDATION_JSON_SCHEMA,
                    "options": {"temperature": 0.01, "num_predict": 4000}
                }
                ollama_resp = requests.post(self.api_url, json=payload, timeout=300)
                ollama_resp.raise_for_status()
                raw_content = ollama_resp.json().get("response", "")

                val_data = self._extract_json_dict(raw_content)

 
                if not val_data:
                    print(f"⚠️  Validation agent failed to output valid JSON. Raw content: {raw_content[:200]}...")
                    val_data = {
                        "biological_plausibility_score": 0.0, "literature_support_score": 0.0,
                        "toxicity_risk_score": 0.5, "mechanistic_coherence_score": 0.0,
                        "safety_score": 0, "critique": f"Validation agent failed to output valid JSON. non-JSON text: {raw_content[:50]}...",
                        "supporting_evidence": "",
                    }

                hypo.update(val_data)

                validated_results.append(hypo)

            except Exception as e:
                print(f"Validation API Error: {e}")
                hypo.update({
                    "biological_plausibility_score": 0.0, "literature_support_score": 0.0,
                    "toxicity_risk_score": 0.5, "mechanistic_coherence_score": 0.0,
                    "safety_score": 0, "critique": f"Validation agent failed to run. No response received.",
                    "supporting_evidence": "",
                })
                validated_results.append(hypo)

        return validated_results

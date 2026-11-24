class RankingAgent:
    def __init__(self, kg):
        self.kg = kg

    def calculate_graph_support(self, hypothesis):
        """
        Queries the graph to see if the genes/drugs involved have VALIDATED or EXPERIMENTAL edges.
        """
        skip_words = {'a', 'an', 'and', 'the', 'by', 'of', 'from', 'on', 'in', 'to', 'with', 'activator', 'inhibitor', 'inducer', '(e.g.,'}
        keywords = hypothesis.get('combination', '').replace("+", " ").split()
        keywords = [k.strip() for k in keywords if len(k) > 2 and k.lower() not in skip_words]

        score = 0.0
        focus_gene = hypothesis.get('focus_gene', '')
        cancer_type = hypothesis.get('target_disease', '')

        query = """
        MATCH path = (d:Drug)-[r*1..2]-(target)
        WHERE toLower(d.name) IN $kws
          AND (toLower(target.name) CONTAINS toLower($gene)
              OR toLower(target.name) CONTAINS toLower($disease))
        UNWIND relationships(path) AS rel
        RETURN type(rel) AS rel_type, rel.layer AS layer //, length(path) AS hop_count, rel.source as source
        LIMIT 50
        """
        with self.kg.driver.session() as session:
            kws_lower = [k.lower() for k in keywords]
            results = session.run(query, kws=kws_lower, gene=focus_gene, disease=cancer_type).data()

            for res in results:
                if res['layer'] == 'VALIDATED':
                    score += 2.0
                elif res['layer'] == 'EXPERIMENTAL':
                    score += 1
                else:
                    score += 0.6

        return min(score, 10.0)

    def rank(self, hypotheses):
        print("\nRanking Agent is sorting hypotheses...")

        for h in hypotheses:
            # Pull structured scores from ValidationAgent (0-1 scale each)
            bio_plausibility  = float(h.get("biological_plausibility_score", 0.0))
            lit_support       = float(h.get("literature_support_score", 0.0))
            toxicity_risk     = float(h.get("toxicity_risk_score", 0.5))  # higher = worse
            mech_coherence    = float(h.get("mechanistic_coherence_score", 0.0))
            graph_score       = self.calculate_graph_support(h)
            graph_score_norm  = graph_score / 10.0   # normalise to 0-1

            # Composite formula — fully deterministic, no LLM involvement
            # toxicity_risk is inverted: low toxicity = higher score
            composite = (
                0.35 * graph_score_norm +
                0.25 * lit_support +
                0.20 * (1.0 - toxicity_risk) +    # invert: low toxicity is good
                0.10 * ((bio_plausibility + mech_coherence) / 2.0) +
                0.10 * (h.get("safety_score", 5) / 10.0)  # NOT inverted — validator already scores high=safe
            )
            h["graph_evidence_score"] = round(graph_score, 2)
            h["ranking_score"]        = round(composite * 10, 2)  # back to 0-10 scale

            method = h["method"]

            if method == "LLM-only":
                h['rationale_type'] = "LLM Internal Knowledge"
            elif method == "KG-Only" or graph_score >= 5:
                h['rationale_type'] = "Direct Graph Evidence"
            elif "Inferred" in h.get('rationale_type', '') or graph_score < 5:
                 h['rationale_type'] = "Inferred Similarity (Embeddings)"


            ne_hits  = h.get('ne_hits', 0)
            nc_hits  = h.get('nc_hits', 0)
            score    = h["ranking_score"]
            withdrawn = h.get("is_withdrawn", False) or h.get("safety_score", 10) <= 2

            if withdrawn:
                calculated_verdict = "Reject"
            elif score >= 8.0:
                calculated_verdict = "Proceed"
            elif ne_hits == 0 and nc_hits == 0:
                if score >= 5.0:
                    calculated_verdict = "Proceed (Novel - Caution)"
                else:
                    calculated_verdict = "Caution (Novel - Low Evidence)"
            elif score >= 5.0:
                calculated_verdict = "Caution"
            else:
                calculated_verdict = "Reject"

            h['verdict'] = calculated_verdict

        hypotheses_sorted = sorted(hypotheses, key=lambda x: x["ranking_score"], reverse=True)

        return hypotheses_sorted
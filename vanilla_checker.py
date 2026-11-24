"""
VANILLA-style logical consistency checker.

Based on: Purohit et al. (2025) "VANILLA: Validated Knowledge Graph Completion —
A Normalization-based framework for Integrity, Link Prediction, and Logical Accuracy."
Knowledge Based Systems, 325. https://doi.org/10.1016/j.knosys.2025.113939

Purpose: Filter predicted paths (from Layer 2) that contain internal logical
contradictions before they are passed to the LLM.

VANILLA defines "plausibility" as logical consistency between relations in a
multi-hop chain. We implement the core idea: if two relations in the same path
are semantically contradictory (e.g. Drug INHIBITS Gene AND Drug ACTIVATES Gene),
the path is inconsistent and should be removed or flagged.
"""


# ---------------------------------------------------------------------------
# Semantic consistency rules
# Each entry: (relation_A, relation_B) → these two relations on the same pair
# of nodes are contradictory.
# ---------------------------------------------------------------------------
CONTRADICTION_PAIRS = {
    frozenset({"INHIBITS",  "ACTIVATES"}),
    frozenset({"TREATS",    "CAUSES"}),
    frozenset({"PROMOTES",  "SUPPRESSES"}),
    frozenset({"UPREGULATES", "DOWNREGULATES"}),
    frozenset({"INDUCES",   "BLOCKS"}),
}

# ---------------------------------------------------------------------------
# Transitivity rules  (head_rel, tail_rel) → implied composite relation
# Used to check if a predicted chain is directionally coherent.
# E.g. Drug INHIBITS Gene + Gene ACTIVATES Pathway → Drug net-SUPPRESSES Pathway
# ---------------------------------------------------------------------------
TRANSITIVITY_MAP = {
    ("INHIBITS",  "ACTIVATES"):    "SUPPRESSES",    # inhibit→activate = suppress downstream
    ("INHIBITS",  "PROMOTES"):     "SUPPRESSES",
    ("ACTIVATES", "PROMOTES"):     "PROMOTES",
    ("ACTIVATES", "INHIBITS"):     "SUPPRESSES",
    ("INHIBITS",  "INHIBITS"):     "ACTIVATES",     # double inhibition = activation
    ("TREATS",    "PROMOTES"):     "CONTRADICTS",   # treat then promote = incoherent
}


class VanillaConsistencyChecker:
    """
    Applies VANILLA-style plausibility rules to a list of paths
    (each path is a dict from _build_raw_paths or build_predicted_paths).

    Returns:
        consistent_paths   — paths that pass all checks
        flagged_paths      — paths with inconsistency details added
    """

    def __init__(self, kg=None):
        # Optionally pass the KG to check against validated facts already in the graph
        self.kg = kg

    def check(self, paths):
        consistent = []
        flagged    = []
        for path in paths:
            issues = self._check_single(path)
            if issues:
                path = dict(path)   # don't mutate original
                path["vanilla_issues"] = issues
                path["vanilla_passed"] = False
                flagged.append(path)
            else:
                path = dict(path)
                path["vanilla_passed"] = True
                consistent.append(path)
        return consistent, flagged

    def _check_single(self, path):
        issues = []
        rel = path.get("relation", "").upper()

        # Rule 1: Direct contradiction between the primary edge and any
        # supplementary relation stored in the path dict.
        alt_rel = path.get("alt_relation", "").upper()
        if alt_rel:
            pair = frozenset({rel, alt_rel})
            if pair in CONTRADICTION_PAIRS:
                issues.append(
                    f"Contradiction: {rel} and {alt_rel} on "
                    f"{path['subject']} → {path['object']}"
                )

        # Rule 2: Multi-hop transitivity check for predicted paths.
        # Predicted paths carry 'intermediate' and 'pathway' fields.
        if path.get("layer") == "PREDICTED":
            hop2_rel = path.get("hop2_relation", "").upper()
            composite = TRANSITIVITY_MAP.get((rel, hop2_rel))
            if composite == "CONTRADICTS":
                issues.append(
                    f"Transitivity incoherence: "
                    f"{path['subject']} {rel} {path['intermediate']} "
                    f"then {path['intermediate']} {hop2_rel} {path.get('pathway', '?')} "
                    f"yields an incoherent composite effect."
                )

        # Rule 3: Cross-check against validated KG facts (if KG available).
        # If the KG already has an edge that CONTRADICTS the predicted edge, flag it.
        if self.kg and path.get("layer") == "PREDICTED":
            contradiction = self._kg_contradiction_check(
                path["subject"], path["relation"], path["object"]
            )
            if contradiction:
                issues.append(
                    f"KG contradiction: graph contains "
                    f"'{contradiction}' between {path['subject']} "
                    f"and {path['object']}, which contradicts predicted {rel}."
                )

        return issues

    def _kg_contradiction_check(self, subject, relation, obj):
        """
        Checks if the KG already stores a relation between subject and object
        that is semantically contradictory to 'relation'.
        Returns the contradicting relation type if found, else None.
        """
        query = """
        MATCH (s)-[r]-(o)
        WHERE toLower(s.name) CONTAINS toLower($subject)
          AND toLower(o.name) CONTAINS toLower($object)
        RETURN type(r) AS rel_type
        LIMIT 20
        """
        try:
            with self.kg.driver.session() as session:
                results = session.run(query, subject=subject, object=obj).data()
            for row in results:
                existing = row["rel_type"].upper()
                pair = frozenset({relation.upper(), existing})
                if pair in CONTRADICTION_PAIRS:
                    return existing
        except Exception:
            pass
        return None

    def summary(self, consistent, flagged):
        total = len(consistent) + len(flagged)
        print(f"   VANILLA check: {len(consistent)}/{total} paths passed, "
              f"{len(flagged)} flagged for inconsistency.")
        for p in flagged[:3]:   # print first 3 issues
            print(f"   ⚠ {p.get('vanilla_issues', [])}")
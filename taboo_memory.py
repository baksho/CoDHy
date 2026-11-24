"""
taboo_memory.py — Agentic memory (Taboo list) for CoDHy

Purpose (per project agenda item 5): remember drug-combination hypotheses
that were previously rejected — either because ValidationAgent scored them
below a threshold, or because FeedbackAgent flagged them as weak — so the
generator doesn't keep re-proposing the same rejected ideas across runs.

Design notes / assumptions (flagged explicitly, please confirm):
- Persists to a JSON file on disk, keyed by (focus_gene, target_disease).
  This is a deliberate choice over Neo4j: it's simpler, needs no schema
  change to build_kg.py, and "memory" for an agentic-memory feature usually
  means surviving across separate runs/sessions — a plain file does that.
  If you'd rather this live in Neo4j as its own node/relationship type
  instead, that's a straightforward swap — let me know.
- Matching is by normalized combination string (case-insensitive, drug
  order-independent — "DrugA + DrugB" and "DrugB + DrugA" are treated as
  the same taboo entry) rather than fuzzy matching, to avoid accidentally
  blocking legitimately different combinations that merely share a word.
"""
import os
import json
from datetime import datetime, timezone


class TabooMemory:
    def __init__(self, storage_path="taboo_memory.json"):
        self.storage_path = storage_path
        self._data = self._load()

    # ------------------------------------------------------------------ #
    #  Persistence                                                        #
    # ------------------------------------------------------------------ #
    def _load(self):
        if os.path.exists(self.storage_path):
            try:
                with open(self.storage_path, "r") as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError) as e:
                print(f"⚠️  Taboo memory file corrupt/unreadable ({e}); starting fresh.")
        return {}

    def _save(self):
        with open(self.storage_path, "w") as f:
            json.dump(self._data, f, indent=2)

    # ------------------------------------------------------------------ #
    #  Key normalization                                                 #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _context_key(focus_gene, target_disease):
        return f"{focus_gene.strip().lower()}::{target_disease.strip().lower()}"

    @staticmethod
    def _normalize_combination(combination):
        """
        "DrugA + DrugB" and "DrugB + DrugA" -> the same canonical key,
        so re-ordered re-proposals are still caught.
        """
        parts = [p.strip().lower() for p in combination.split("+")]
        return " + ".join(sorted(parts))

    # ------------------------------------------------------------------ #
    #  Write                                                             #
    # ------------------------------------------------------------------ #
    def add_rejected(self, focus_gene, target_disease, combination, reason, score=None):
        ctx_key = self._context_key(focus_gene, target_disease)
        combo_key = self._normalize_combination(combination)

        self._data.setdefault(ctx_key, {})
        self._data[ctx_key][combo_key] = {
            "original_combination": combination,
            "reason": reason,
            "score": score,
            "rejected_at": datetime.now(timezone.utc).isoformat(),
        }
        self._save()

    # ------------------------------------------------------------------ #
    #  Read                                                              #
    # ------------------------------------------------------------------ #
    def is_taboo(self, focus_gene, target_disease, combination):
        ctx_key = self._context_key(focus_gene, target_disease)
        combo_key = self._normalize_combination(combination)
        return combo_key in self._data.get(ctx_key, {})

    def get_taboo_list(self, focus_gene, target_disease):
        """
        Returns a list of {"combination": ..., "reason": ...} dicts for the
        given gene/disease context — used to build the "do not repropose"
        block in DrugScreeningAgent's prompt.
        """
        ctx_key = self._context_key(focus_gene, target_disease)
        entries = self._data.get(ctx_key, {})
        return [
            {"combination": v["original_combination"], "reason": v["reason"]}
            for v in entries.values()
        ]

    def count(self, focus_gene=None, target_disease=None):
        if focus_gene is not None and target_disease is not None:
            ctx_key = self._context_key(focus_gene, target_disease)
            return len(self._data.get(ctx_key, {}))
        return sum(len(v) for v in self._data.values())

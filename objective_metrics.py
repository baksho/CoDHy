"""
objective_metrics.py — Table 3 objective metrics for CoDHy's ICDM evaluation

Supervisor's point 6, Table 3: "Evidence Coverage, Novel Combination Rate,
Diversity, Citation Precision, and Runtime" across systems.

Where each metric comes from:
  - Evidence Coverage & Novel Combination Rate: computed directly from
    fields hvalidator.py ALREADY writes onto every hypothesis dict
    (ne_hits, nc_hits, is_novel_exact, is_novel_context). No new
    instrumentation needed for CoDHy's own hypotheses. Baseline hypotheses
    from baseline_generators.py do NOT go through hvalidator.py (they're
    zero-shot, no KG lookup), so these two metrics are only meaningful for
    hypotheses that were actually validated — see the `requires_validation`
    note in build_objective_metrics_table().
  - Diversity: computed from the `combination` field text alone — works for
    both CoDHy and baseline hypotheses.
  - Citation Precision: cross-checks any URL-like substring found in
    `supporting_evidence` against the URL lists hvalidator.py already
    attaches (civic_urls/chembl_urls/trial_urls/reactome_urls/source_urls).
    Only meaningful for validated (CoDHy) hypotheses for the same reason as
    above — a baseline hypothesis has no such "known-good" URL list to check
    against, so by construction any URL it cites can't be verified this way.
    If you want Citation Precision for baselines too, you'd need a live
    lookup (e.g. query PubMed for each cited PMID) — not implemented here,
    flagging rather than guessing at what "verified" should mean for that
    case.
  - Runtime: NOT derivable from the hypothesis dicts — needs timing
    instrumentation added to run_full_experiment() in the notebook (small
    diff provided separately). Pass it in as `runtime_seconds` per system.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd

_URL_RE = re.compile(r'https?://\S+')


def evidence_coverage(hypotheses: list[dict]) -> float:
    """Fraction of hypotheses with at least one real literature hit
    (ne_hits or nc_hits > 0, as recorded by hvalidator.check_combination_evidence).
    Requires hypotheses to have gone through ValidationAgent.validate_batch()."""
    if not hypotheses:
        return np.nan
    covered = sum(1 for h in hypotheses if (h.get("ne_hits", 0) or 0) > 0 or (h.get("nc_hits", 0) or 0) > 0)
    return covered / len(hypotheses)


def novel_combination_rate(hypotheses: list[dict], novelty_field: str = "is_novel_exact") -> float:
    """
    Fraction of hypotheses flagged as novel by ValidationAgent. Default
    novelty_field="is_novel_exact" (no exact-combination literature hit —
    the stricter definition). Use "is_novel_context" for the looser
    definition (no hit even when co-occurring with the disease term).
    """
    if not hypotheses:
        return np.nan
    novel = sum(1 for h in hypotheses if h.get(novelty_field) is True)
    return novel / len(hypotheses)


def diversity(hypotheses: list[dict]) -> float:
    """
    Unique combination strings / total hypotheses. 1.0 = every hypothesis
    proposes a distinct combination; lower values indicate repeated/
    near-duplicate proposals. Case-insensitive, whitespace-normalized.
    """
    if not hypotheses:
        return np.nan
    combos = [re.sub(r'\s+', ' ', h.get("combination", "")).strip().lower() for h in hypotheses]
    combos = [c for c in combos if c]
    if not combos:
        return np.nan
    return len(set(combos)) / len(combos)


def citation_precision(hypotheses: list[dict]) -> float:
    """
    Of all URL-like substrings appearing in supporting_evidence across the
    given hypotheses, what fraction match a URL already present in that same
    hypothesis's own civic_urls/chembl_urls/trial_urls/reactome_urls/
    source_urls (i.e. a URL hvalidator.py actually retrieved from the KG /
    PubMed — not one the LLM invented)? Hypotheses with zero cited URLs are
    excluded from both numerator and denominator (vacuous — nothing to
    verify), matching the pipeline's own "no direct clinical study found"
    convention for uncited claims.
    """
    total_cited = 0
    total_verified = 0
    for h in hypotheses:
        known_urls = set()
        for field in ("civic_urls", "chembl_urls", "trial_urls", "reactome_urls", "source_urls"):
            known_urls.update(h.get(field, []) or [])
        cited = _URL_RE.findall(h.get("supporting_evidence", "") or "")
        for url in cited:
            total_cited += 1
            if url.rstrip('.,);') in known_urls or any(url.startswith(k) for k in known_urls):
                total_verified += 1
    if total_cited == 0:
        return np.nan
    return total_verified / total_cited


def build_objective_metrics_table(system_to_hyps: dict[str, list[dict]],
                                   runtime_seconds: dict[str, float] = None) -> pd.DataFrame:
    """
    Table 3. system_to_hyps: {system_name: [hypothesis dicts]}, e.g.
        {"CoDHy-Full": ranked_hypotheses, "DeepSeek-R1-ZeroShot": [...], ...}
    runtime_seconds: optional {system_name: seconds}, merged in as a column.

    Evidence Coverage / Novel Combination Rate / Citation Precision are only
    meaningful for hypotheses that went through ValidationAgent — they'll
    show NaN for baseline systems that skip validation by design, which is
    expected, not a bug. Report that limitation explicitly in the paper
    rather than silently filling zeros.
    """
    runtime_seconds = runtime_seconds or {}
    rows = []
    for system_name, hyps in system_to_hyps.items():
        validated = any(("ne_hits" in h or "civic_urls" in h) for h in hyps) if hyps else False
        rows.append({
            "system": system_name,
            "n_hypotheses": len(hyps),
            "evidence_coverage": evidence_coverage(hyps) if validated else np.nan,
            "novel_combination_rate": novel_combination_rate(hyps) if validated else np.nan,
            "diversity": diversity(hyps),
            "citation_precision": citation_precision(hyps) if validated else np.nan,
            "runtime_seconds": runtime_seconds.get(system_name, np.nan),
            "validated_pipeline": validated,
        })
    return pd.DataFrame(rows)

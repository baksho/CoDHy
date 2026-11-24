"""
judge_metrics.py — Inter-judge agreement & reporting for CoDHy's LLM-as-Judge
evaluation.

Given the long-format DataFrame produced by llm_judge.run_judge_panel()
(one row per hypothesis x judge, with columns from JUDGE_SCORE_FIELDS), this
module computes the statistics an ICDM reviewer will expect from an
LLM-as-Judge evaluation:

  1. Krippendorff's alpha (interval level) per scoring dimension — the
     standard chance-corrected agreement statistic for >2 raters on
     ordinal/interval data. Threshold conventions from the reliability
     literature (Krippendorff himself): alpha >= 0.800 = reliable,
     0.667 <= alpha < 0.800 = tentative conclusions only, below that =
     unreliable. These are conventional cutoffs, not universal law — state
     them as such if you cite them.
  2. ICC — intraclass correlation (pingouin), reported as both ICC(C,k)
     ("consistency", average of k raters — appropriate here since your
     panel is a fixed, specific set of 5 named judges rather than a random
     sample of possible judges) and ICC(A,k) ("absolute agreement") for
     completeness. This is a genuine judgment call in the stats literature
     (fixed vs random effects); note the choice explicitly in your paper.
  3. Per-judge bias check: mean/std of each judge's scores per dimension,
     to flag if one judge is a systematic outlier (e.g. always scores high).
  4. Pairwise Spearman correlation between judges, per dimension.
  5. Self-vs-external comparison: correlates the ValidationAgent's own
     (in-pipeline) scores against the external judge panel's MEAN score, to
     report how much the self-validation step agrees with independent judges.

pip install krippendorff pingouin pandas numpy scipy   (all confirmed
importable at the time this module was written).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    import krippendorff
except ImportError as e:
    raise ImportError(
        "krippendorff package not found. Install with: pip install krippendorff"
    ) from e

try:
    import pingouin as pg
except ImportError as e:
    raise ImportError(
        "pingouin package not found. Install with: pip install pingouin"
    ) from e


# Conventional Krippendorff's-alpha interpretation bands (Krippendorff, as
# commonly cited in reliability-methods papers). State as convention, not law.
ALPHA_RELIABLE = 0.800
ALPHA_TENTATIVE = 0.667


def diagnose_judge_failures(judge_df: pd.DataFrame) -> pd.DataFrame:
    """
    Per-judge call failure rate (from the `_parse_failed` flag llm_judge.py
    sets on every row) plus, for reference, the per-judge mean of
    overall_quality among only the successful calls. Run this FIRST whenever
    a metric looks degenerate (all-NaN, single unique value, a crash) — it
    tells you whether a judge model is systematically failing to respond
    (e.g. not actually served by HF Inference Providers) rather than the
    hypotheses genuinely being scored identically.
    """
    if "_parse_failed" not in judge_df.columns:
        raise ValueError("judge_df has no '_parse_failed' column — was it produced by llm_judge.run_judge_panel()?")

    summary = judge_df.groupby("judge_name").agg(
        n_calls=("_parse_failed", "size"),
        n_failed=("_parse_failed", "sum"),
    )
    summary["failure_rate"] = (summary["n_failed"] / summary["n_calls"]).round(2)
    ok = judge_df[judge_df["_parse_failed"] == False]
    if "overall_quality" in judge_df.columns and not ok.empty:
        summary["mean_overall_quality_on_success"] = ok.groupby("judge_name")["overall_quality"].mean().round(2)
    return summary.reset_index()


def krippendorff_alpha_by_dimension(judge_df: pd.DataFrame, dimensions: list[str]) -> pd.DataFrame:
    """
    Computes Krippendorff's alpha (interval level) for each scoring
    dimension. judge_df must be long-format with columns:
    ['hypothesis_id', 'judge_name', <dimension columns...>].
    """
    results = []
    hyp_ids = sorted(judge_df["hypothesis_id"].dropna().unique())
    judges = sorted(judge_df["judge_name"].dropna().unique())

    for dim in dimensions:
        # Build (M judges x N hypotheses) reliability matrix, np.nan for missing/failed.
        matrix = np.full((len(judges), len(hyp_ids)), np.nan)
        for j_idx, judge in enumerate(judges):
            for h_idx, hid in enumerate(hyp_ids):
                cell = judge_df[(judge_df["judge_name"] == judge) & (judge_df["hypothesis_id"] == hid)]
                if not cell.empty:
                    val = cell.iloc[0][dim]
                    if pd.notna(val):
                        matrix[j_idx, h_idx] = float(val)

        n_valid = np.isnan(matrix).sum()
        distinct_values = np.unique(matrix[~np.isnan(matrix)])

        if len(distinct_values) < 2:
            # krippendorff.alpha() raises ValueError here rather than
            # returning e.g. NaN — guard explicitly instead of crashing the
            # whole report. This happens when every non-missing score for
            # this dimension is identical (no disagreement to measure) OR
            # when so many judge calls failed that too little real data
            # remains. Which one it is matters a lot for interpretation —
            # see diagnose_judge_failures() below to tell them apart.
            results.append({
                "dimension": dim, "krippendorff_alpha": np.nan,
                "interpretation": "undefined (fewer than 2 distinct non-missing scores — "
                                   "either genuine unanimous agreement or too many failed judge calls; "
                                   "run diagnose_judge_failures(judge_df) to check)",
            })
            continue

        alpha = krippendorff.alpha(reliability_data=matrix, level_of_measurement="interval")
        band = (
            "reliable" if alpha >= ALPHA_RELIABLE
            else "tentative" if alpha >= ALPHA_TENTATIVE
            else "unreliable"
        )
        results.append({"dimension": dim, "krippendorff_alpha": alpha, "interpretation": band})

    return pd.DataFrame(results)


def icc_by_dimension(judge_df: pd.DataFrame, dimensions: list[str]) -> pd.DataFrame:
    """
    Computes ICC(C,k) and ICC(A,k) per dimension via pingouin. Rows with a
    missing score for any judge/hypothesis pair are omitted listwise
    (nan_policy='omit') rather than raising.
    """
    rows = []
    for dim in dimensions:
        sub = judge_df[["hypothesis_id", "judge_name", dim]].dropna()
        sub = sub.rename(columns={dim: "score"})
        if sub["hypothesis_id"].nunique() < 2 or sub["judge_name"].nunique() < 2:
            rows.append({"dimension": dim, "ICC(C,k)": np.nan, "ICC(A,k)": np.nan,
                         "note": "insufficient data (need >=2 hypotheses and >=2 judges with non-null scores)"})
            continue
        icc_table = pg.intraclass_corr(
            data=sub, targets="hypothesis_id", raters="judge_name", ratings="score",
            nan_policy="omit",
        )
        ck = icc_table.loc[icc_table["Type"] == "ICC(C,k)", "ICC"]
        ak = icc_table.loc[icc_table["Type"] == "ICC(A,k)", "ICC"]
        rows.append({
            "dimension": dim,
            "ICC(C,k)": float(ck.iloc[0]) if not ck.empty else np.nan,
            "ICC(A,k)": float(ak.iloc[0]) if not ak.empty else np.nan,
            "note": "",
        })
    return pd.DataFrame(rows)


def judge_bias_report(judge_df: pd.DataFrame, dimensions: list[str]) -> pd.DataFrame:
    """Mean/std per judge per dimension — flags systematically lenient/harsh judges.
    This is Table 2 (tab:judge-diagnostics) in the paper: "each judge's own
    scoring tendency (mean +/- std across every hypothesis it judged)"."""
    return (
        judge_df.groupby("judge_name")[dimensions]
        .agg(["mean", "std"])
        .round(2)
    )


def system_median_scores(judge_df: pd.DataFrame, dimensions: list[str]) -> pd.DataFrame:
    """
    Median score per SYSTEM per dimension, across whichever judges actually
    scored each hypothesis (self-excluded judges never produced a row for
    their own baseline's hypotheses in the first place, so a plain groupby
    here already respects that — no extra filtering needed). This is Table 1
    (tab:judge-results) in the paper: CoDHy is scored by all 6 judges;
    each baseline's own hypotheses are scored by the remaining 5 (the 6th,
    that baseline's own generating model, abstains).
    """
    return judge_df.groupby("system_name")[dimensions].median().round(2)


def build_table1(judge_df: pd.DataFrame, dimensions: list[str]) -> pd.DataFrame:
    """
    Table 1 exactly as captioned: system_median_scores() with a final
    "Krippendorff's alpha" row appended (one alpha value per dimension,
    computed once across the full judge x hypothesis matrix — alpha is a
    property of inter-judge agreement overall, not per-system, so this row
    is the same regardless of which system's rows sit above it).
    """
    medians = system_median_scores(judge_df, dimensions)
    alpha_df = krippendorff_alpha_by_dimension(judge_df, dimensions)
    alpha_row = pd.DataFrame(
        [alpha_df.set_index("dimension")["krippendorff_alpha"].round(3)],
        index=["Krippendorff's alpha"],
    )
    return pd.concat([medians, alpha_row[medians.columns]])


def pairwise_spearman(judge_df: pd.DataFrame, dimension: str) -> pd.DataFrame:
    """
    Pairwise Spearman rank correlation between every pair of judges, for one
    scoring dimension. Useful as a supplementary/robustness check alongside
    Krippendorff's alpha and ICC (rank-based rather than magnitude-based).
    """
    wide = judge_df.pivot_table(index="hypothesis_id", columns="judge_name", values=dimension)
    judges = wide.columns.tolist()
    rows = []
    for i in range(len(judges)):
        for j in range(i + 1, len(judges)):
            pair = wide[[judges[i], judges[j]]].dropna()
            if len(pair) < 3:
                rho, p = np.nan, np.nan
            else:
                rho, p = spearmanr(pair[judges[i]], pair[judges[j]])
            rows.append({"judge_a": judges[i], "judge_b": judges[j], "dimension": dimension,
                         "spearman_rho": rho, "p_value": p, "n_pairs": len(pair)})
    return pd.DataFrame(rows)


def compare_to_self_validation(judge_df: pd.DataFrame, hypotheses: list[dict],
                                field_map: dict = None) -> pd.DataFrame:
    """
    Correlates the ValidationAgent's own in-pipeline scores (already present
    on each hypothesis dict from hvalidator.py) against the judge panel's
    MEAN external score, per mapped dimension. Reports how well self-scoring
    agrees with independent judges — a useful sanity/validity check to
    report in the paper.

    field_map: {self_validation_field: judge_dimension}. Defaults to the
    fields that exist in both hvalidator.py's OUTPUT FORMAT and
    JUDGE_SCORE_FIELDS in llm_judge.py. NOTE: ValidationAgent's
    biological_plausibility_score/literature_support_score are on a 0.0-1.0
    scale while judge scores are 0-10 — this function rescales the
    self-validation score to 0-10 before correlating, since Spearman
    correlation is scale-invariant anyway but keeping ranges consistent
    makes the plotted comparison meaningful too.
    """
    # NOTE: the judge rubric now uses the supervisor-specified 7 dimensions
    # on a 1-5 scale (see llm_judge.JUDGE_SCORE_FIELDS); "safety_reasonableness"
    # and "overall_quality" from the earlier 0-10 rubric no longer exist on
    # the judge side, so only biological_plausibility has a direct
    # self-validation counterpart. Add more mappings here if you reintroduce
    # matching dimensions to the judge rubric.
    field_map = field_map or {
        "biological_plausibility_score": "biological_plausibility",
    }

    hyp_by_id = {h.get("id"): h for h in hypotheses}
    judge_means = judge_df.groupby("hypothesis_id")[list(set(field_map.values()))].mean()

    rows = []
    for self_field, judge_field in field_map.items():
        pairs = []
        for hid, jrow in judge_means.iterrows():
            h = hyp_by_id.get(hid)
            if h is None or self_field not in h:
                continue
            self_val = h[self_field]
            if self_val is None:
                continue
            # Rescale 0.0-1.0 self-validation fields to the judge rubric's
            # 1-5 scale for comparability (Spearman rho is scale-invariant,
            # but this keeps any plotted comparison on matching axes).
            if self_field.endswith("_score") and self_field != "safety_score" and 0 <= self_val <= 1:
                self_val = 1 + self_val * 4
            pairs.append((self_val, jrow[judge_field]))

        if len(pairs) >= 3:
            self_vals, judge_vals = zip(*pairs)
            rho, p = spearmanr(self_vals, judge_vals)
        else:
            rho, p = np.nan, np.nan

        rows.append({
            "self_validation_field": self_field,
            "judge_dimension": judge_field,
            "n_hypotheses": len(pairs),
            "spearman_rho": rho,
            "p_value": p,
        })
    return pd.DataFrame(rows)


def repeat_consistency(judge_df: pd.DataFrame, dimensions: list[str]) -> pd.DataFrame:
    """
    Test-retest reliability (supervisor point 5): for judge_df produced with
    n_repeats > 1, this is agreement ACROSS REPEATS of the SAME judge on the
    SAME hypothesis — distinct from krippendorff_alpha_by_dimension, which
    measures agreement ACROSS JUDGES. Low repeat-consistency for a judge
    means its scores are noisy/unstable even before you compare it to other
    judges. Reports std-across-repeats averaged over (hypothesis, judge)
    pairs, per dimension — lower is more consistent (0 = perfectly stable).

    No-op (returns a note) if judge_df has no repeat_idx column, i.e. it was
    produced with the default n_repeats=1.
    """
    if "repeat_idx" not in judge_df.columns or judge_df["repeat_idx"].nunique() < 2:
        return pd.DataFrame([{"note": "no repeat_idx with >=2 repeats found — "
                                       "re-run run_judge_panel(..., n_repeats>=2) to compute this."}])

    rows = []
    for dim in dimensions:
        std_per_cell = (
            judge_df.groupby(["hypothesis_id", "judge_name"])[dim]
            .std()
            .dropna()
        )
        rows.append({
            "dimension": dim,
            "mean_std_across_repeats": round(float(std_per_cell.mean()), 3) if not std_per_cell.empty else np.nan,
            "n_hypothesis_judge_pairs": len(std_per_cell),
        })
    return pd.DataFrame(rows)


def pairwise_report(pairwise_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates llm_judge.run_pairwise_panel()'s long-format output into
    win/tie/loss percentages per (system_a, system_b) pair (supervisor
    point 4, Table 2). "Win" is always reported from system_a's perspective
    per row of the input df's own (system_a, system_b) labeling — the
    caller controls that labeling by how they invoke run_pairwise_panel.
    Rows where the judge call failed (_parse_failed=True) are excluded from
    the percentages but counted separately so you can see the failure rate.
    """
    rows = []
    for (sys_a, sys_b), grp in pairwise_df.groupby(["system_a", "system_b"]):
        valid = grp[~grp["_parse_failed"]]
        n_valid = len(valid)
        n_failed = len(grp) - n_valid
        if n_valid == 0:
            rows.append({"system_a": sys_a, "system_b": sys_b, "n_comparisons": 0,
                         "n_failed": n_failed, "win_pct": np.nan, "tie_pct": np.nan, "loss_pct": np.nan})
            continue
        wins = (valid["winning_system"] == sys_a).sum()
        losses = (valid["winning_system"] == sys_b).sum()
        ties = (valid["winning_system"] == "Tie").sum()
        rows.append({
            "system_a": sys_a, "system_b": sys_b, "n_comparisons": n_valid, "n_failed": n_failed,
            "win_pct": round(100 * wins / n_valid, 1),
            "tie_pct": round(100 * ties / n_valid, 1),
            "loss_pct": round(100 * losses / n_valid, 1),
        })
    return pd.DataFrame(rows)


def build_full_report(judge_df: pd.DataFrame, hypotheses: list[dict],
                       dimensions: list[str] = None) -> dict:
    """
    Convenience wrapper: runs all of the above and returns a dict of
    DataFrames, plus markdown-ready tables for pasting into a paper draft.
    Requires the 'tabulate' package for .to_markdown() (pip install tabulate)
    — falls back to plain string repr if unavailable.
    """
    from llm_judge import JUDGE_SCORE_FIELDS
    dimensions = dimensions or [d for d in JUDGE_SCORE_FIELDS if d != "justification"]

    report = {
        "judge_failure_diagnostics": diagnose_judge_failures(judge_df),
        "table1_system_median_scores": build_table1(judge_df, dimensions),
        "krippendorff_alpha": krippendorff_alpha_by_dimension(judge_df, dimensions),
        "icc": icc_by_dimension(judge_df, dimensions),
        "judge_bias": judge_bias_report(judge_df, dimensions),
        "repeat_consistency": repeat_consistency(judge_df, dimensions),
        "spearman_pairwise_judges": pd.concat(
            [pairwise_spearman(judge_df, d) for d in dimensions if d in judge_df.columns],
            ignore_index=True,
        ),
        "self_vs_external": compare_to_self_validation(judge_df, hypotheses),
    }

    print("\n" + "=" * 70)
    print("  LLM-AS-JUDGE EVALUATION REPORT")
    print("=" * 70)
    for name, df in report.items():
        print(f"\n--- {name} ---")
        try:
            print(df.to_markdown(index=(name == "judge_bias")))
        except ImportError:
            print(df)
    print("\n" + "=" * 70)

    return report

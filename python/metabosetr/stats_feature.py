"""Approach A: per-metabolite group testing (port of ``stats-feature.R``).

DESIGN.md §6.1. Each metabolite is an individual feature. Default engine is
the empirical-Bayes moderated t via InMoose (a drop-in Python port of R
limma, near-identical results); Wilcoxon rank-sum via SciPy is the support
option (MIGRATION.md §2 decision A). Concentrations are log2-transformed
with the adaptive per-feature pseudocount before testing; outlier samples
are excluded.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._padjust import p_adjust
from .transform import log2_with_pseudocount


def _resolve_group(group, prep):
    """Resolve ``group`` to a length-``n_samples`` array (before outlier drop).

    Accepts a single ``sample_meta`` column name or a length-n vector.
    """
    meta = prep.sample.sample_meta
    n = prep.sample.assay.shape[0]
    if isinstance(group, str):
        if group not in meta.columns:
            raise ValueError(
                f"`group` column '{group}' not in sample_meta. "
                f"Available: {', '.join(map(str, meta.columns))}"
            )
        return meta[group].to_numpy()
    group = np.asarray(group, dtype=object)
    if group.shape[0] != n:
        raise ValueError(
            "`group` must be a single sample_meta column name or a vector of "
            "length nrow(sample.assay)."
        )
    return group


def _two_levels(group):
    """Sorted, dropped-empty factor levels; error unless exactly two."""
    vals = [g for g in group if g is not None and not (isinstance(g, float) and np.isnan(g))]
    levels = sorted(set(vals))
    if len(levels) != 2:
        raise ValueError(
            f"test requires exactly two non-empty group levels; got {len(levels)}."
        )
    return levels


def _add_padjust(out: pd.DataFrame, padjust) -> pd.DataFrame:
    for m in padjust:
        out[f"padj_{m}"] = p_adjust(out["P.Value"].to_numpy(), m)
    return out


def _test_wilcoxon(assay_log, group, padjust, metabolite_names) -> pd.DataFrame:
    from scipy.stats import mannwhitneyu

    levels = _two_levels(group)
    group = np.asarray(group, dtype=object)
    g1 = group == levels[0]
    g2 = group == levels[1]
    n = assay_log.shape[0]
    pvals = np.full(n, np.nan)
    fcs = np.full(n, np.nan)
    for i in range(n):
        x1 = assay_log[i, g1]
        x2 = assay_log[i, g2]
        x1f = x1[np.isfinite(x1)]
        x2f = x2[np.isfinite(x2)]
        if x1f.size < 1 or x2f.size < 1:
            continue
        try:
            res = mannwhitneyu(
                x2f, x1f, alternative="two-sided", use_continuity=True, method="asymptotic"
            )
            pvals[i] = res.pvalue
        except ValueError:
            pvals[i] = np.nan
        fcs[i] = x2f.mean() - x1f.mean()
    out = pd.DataFrame(
        {"metabolite": metabolite_names, "log2FC": fcs, "P.Value": pvals}
    )
    return _add_padjust(out, padjust)


def _test_limma(assay_log, group, padjust, metabolite_names) -> pd.DataFrame:
    from inmoose.limma import eBayes, lmFit, topTable

    levels = _two_levels(group)
    group = np.asarray(group, dtype=object)
    # design: intercept + indicator for the second (reference-sorted) level,
    # matching R model.matrix(~ group) with alphabetically-sorted levels.
    g2 = (group == levels[1]).astype(float)
    design = np.column_stack([np.ones(len(group)), g2])

    expr = pd.DataFrame(assay_log, index=metabolite_names)
    fit = eBayes(lmFit(expr, design))
    # InMoose labels the design columns column0/column1; the group effect is
    # column1. sort_by="none" is unsupported, so sort then restore feature
    # order via reindex to match R topTable(sort.by="none").
    tt = topTable(fit, coef="column1", number=assay_log.shape[0], sort_by="p")
    tt = tt.reindex(metabolite_names)
    out = pd.DataFrame(
        {
            "metabolite": metabolite_names,
            "log2FC": tt["log2FoldChange"].to_numpy(),
            "AveExpr": tt["AveExpr"].to_numpy(),
            "t": tt["stat"].to_numpy(),
            "P.Value": tt["pvalue"].to_numpy(),
        }
    )
    return _add_padjust(out, padjust)


def test_metabolites(prep, group, method="limma", padjust=("BH", "BY")) -> pd.DataFrame:
    """Feature-level group testing (Approach A).

    ``method`` is ``"limma"`` (default, InMoose moderated t) or
    ``"wilcoxon"``. Returns one row per metabolite.
    """
    if method not in ("limma", "wilcoxon"):
        raise ValueError("method must be 'limma' or 'wilcoxon'.")
    group = _resolve_group(group, prep)
    keep = ~prep.sample_outliers
    group = group[keep]
    assay = prep.sample.assay[keep, :]
    assay_log = log2_with_pseudocount(assay.T)  # rows = metabolites
    names = list(prep.sample.metabolite_names)
    if method == "limma":
        return _test_limma(assay_log, group, padjust, names)
    return _test_wilcoxon(assay_log, group, padjust, names)

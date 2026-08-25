"""Approach B: GSEA over curated pathway sets (port of ``stats-gsea.R``).

DESIGN.md §6.2. Rank statistic = per-metabolite log2FC (group2 vs group1)
over log2-transformed concentrations; outlier samples excluded. The R
version uses ``fgsea``; this port uses ``gseapy`` prerank (MIGRATION.md §2
decision B). Results are not numerically identical to fgsea, but are
deterministic given ``seed``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .pathways import _check_domain, get_pathway_sets
from .stats_feature import _resolve_group, _two_levels
from .transform import log2_with_pseudocount


def _metabolite_log2fc(prep, group) -> pd.Series:
    """Per-metabolite log2FC = mean(g2) - mean(g1); finite entries only."""
    keep = ~prep.sample_outliers
    if keep.sum() < 2:
        raise ValueError(
            "Fewer than two non-outlier samples available; cannot compute "
            "metabolite log2FC."
        )
    group = np.asarray(group, dtype=object)[keep]
    assay = prep.sample.assay[keep, :]
    metabs = list(prep.sample.metabolite_names)
    metab_log = log2_with_pseudocount(assay.T)  # rows = metabolites

    levels = _two_levels(group)
    g1 = group == levels[0]
    g2 = group == levels[1]
    with np.errstate(invalid="ignore"):
        m1 = np.nanmean(metab_log[:, g1], axis=1)
        m2 = np.nanmean(metab_log[:, g2], axis=1)
    fc = m2 - m1
    s = pd.Series(fc, index=metabs)
    return s[np.isfinite(s)]


def gsea_pathway_sets(prep, group, domain, seed: int = 0, **kwargs):
    """GSEA over the built-in curated pathway sets for ``domain``.

    ``seed`` pins the permutation RNG for reproducibility (decision B).
    Extra keyword args pass through to ``gseapy.prerank`` (e.g. ``min_size``,
    ``max_size``, ``permutation_num``). Returns the gseapy ``res2d``
    DataFrame.
    """
    import gseapy as gp

    _check_domain(domain)
    group = _resolve_group(group, prep)
    pathways = get_pathway_sets(domain)
    fc = _metabolite_log2fc(prep, group)

    present = set(fc.index)
    pathways = {sn: [m for m in members if m in present] for sn, members in pathways.items()}
    pathways = {sn: members for sn, members in pathways.items() if members}
    if not pathways:
        raise ValueError(
            "No pathway has any member overlapping the rank vector. Check that "
            f"`prep` covers the metabolites of domain '{domain}'."
        )

    rnk = fc.sort_values(ascending=False)
    params = dict(
        rnk=rnk,
        gene_sets=pathways,
        seed=seed,
        min_size=1,
        max_size=len(rnk),
        permutation_num=1000,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    params.update(kwargs)
    pre = gp.prerank(**params)
    return pre.res2d

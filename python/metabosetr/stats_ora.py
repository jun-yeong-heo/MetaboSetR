"""Approach C: over-representation analysis (port of ``stats-ora.R``).

DESIGN.md §6.6. Threshold-based complement to GSEA. One-sided
hypergeometric test per curated pathway set; background is the measured
panel (``universe``). ``scipy.stats.hypergeom.sf`` reproduces R
``phyper(..., lower.tail = FALSE)`` exactly.
"""

from __future__ import annotations

import warnings

import pandas as pd
from scipy.stats import hypergeom

from ._padjust import p_adjust
from .pathways import _check_domain, get_pathway_sets


def _unique_preserve(seq):
    seen = set()
    out = []
    for v in seq:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            continue
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _intersect_preserve(a, keep):
    keep = set(keep)
    return [v for v in a if v in keep]


def ora_pathway_sets(sig, universe, domain, padjust=("BH", "BY")) -> pd.DataFrame:
    """Over-representation analysis of ``sig`` against a pathway-set domain.

    ``sig`` and ``universe`` are metabolite ``short_name``s. Returns one row
    per pathway with at least one member in ``universe``, sorted by ``pval``.
    """
    _check_domain(domain)
    sig = list(sig)
    universe = list(universe)

    universe = _unique_preserve(universe)
    if len(universe) == 0:
        raise ValueError("`universe` is empty; supply the measured feature names.")
    sig = _unique_preserve(sig)
    off = [s for s in sig if s not in set(universe)]
    if off:
        warnings.warn(
            f"{len(off)} of `sig` not in `universe`; dropping them. "
            "`sig` should be a subset of the measured panel.",
            stacklevel=2,
        )
    sig = _intersect_preserve(sig, universe)

    sets = get_pathway_sets(domain)
    N = len(universe)
    n = len(sig)
    uni_set = set(universe)
    sig_set = set(sig)

    rows = []
    for sn, members in sets.items():
        members_in = _intersect_preserve(members, uni_set)
        M = len(members_in)
        if M == 0:
            continue
        hits = _intersect_preserve(sig, set(members_in) & sig_set)
        k = len(hits)
        pval = float(hypergeom.sf(k - 1, N, M, n))  # P(X >= k)
        rows.append(
            {
                "set_name": sn,
                "set_size": M,
                "overlap": k,
                "expected": n * M / N,
                "pval": pval,
                "overlap_members": ";".join(hits),
            }
        )
    if not rows:
        raise ValueError(
            f"No pathway in domain '{domain}' has any member in `universe`."
        )

    out = pd.DataFrame(rows)
    for m in padjust:
        out[f"padj_{m}"] = p_adjust(out["pval"].to_numpy(), m)

    cols = (
        ["set_name", "set_size", "overlap", "expected", "pval"]
        + [f"padj_{m}" for m in padjust]
        + ["overlap_members"]
    )
    out = out[cols].sort_values("pval", kind="stable").reset_index(drop=True)
    return out

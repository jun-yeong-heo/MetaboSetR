"""Multiple-testing correction (port of R ``stats::p.adjust``).

Implements the BH and BY methods used across the package (decisions.md #29)
without pulling in ``statsmodels`` -- matching the R semantics exactly,
including NA handling (NAs are excluded from ranking and reinserted, and the
default ``n`` is the full length including NAs).
"""

from __future__ import annotations

import numpy as np


def p_adjust(p, method: str) -> np.ndarray:
    """Adjust p-values by ``method`` in {"BH"/"fdr", "BY"}.

    Mirrors ``stats::p.adjust``: work on the non-NA subset, use ``n`` equal
    to the full input length, and put NAs back in place.
    """
    p = np.asarray(p, dtype=float)
    nna = ~np.isnan(p)
    pv = p[nna]
    lp = pv.size
    # R's default `n = length(p)` is a lazy promise evaluated only after
    # `p <- p[!is.na(p)]`, so n is the non-NA count, not the full length.
    n = lp
    out = np.full(p.shape, np.nan)
    if lp == 0:
        return out

    # order by decreasing p (stable), and the inverse permutation.
    o = np.argsort(-pv, kind="stable")
    ro = np.argsort(o, kind="stable")
    i = np.arange(lp, 0, -1)  # lp, lp-1, ..., 1

    if method in ("BH", "fdr"):
        factor = n / i
    elif method == "BY":
        q = np.sum(1.0 / np.arange(1, n + 1))
        factor = q * n / i
    else:
        raise ValueError(f"Unsupported p.adjust method: {method}")

    adj = np.minimum(1.0, np.minimum.accumulate(factor * pv[o]))[ro]
    out[nna] = adj
    return out

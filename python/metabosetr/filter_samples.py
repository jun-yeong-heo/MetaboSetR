"""Sample-level outlier detection (port of ``filter-samples.R``).

DESIGN.md §4.2. Per-sample ``<LOD`` / ``>ULOQ`` rates drive outlier
flagging via IQR / percentile / both. ``numpy`` linear quantile matches R's
default ``quantile`` type 7.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _apply_outlier_method(rates: np.ndarray, method: str, percentile_cutoff: float) -> np.ndarray:
    rates = np.asarray(rates, dtype=float)
    finite = ~np.isnan(rates)

    def iqr_outlier(x):
        q1, q3 = np.nanquantile(x, [0.25, 0.75])  # type 7 (linear)
        threshold = q3 + 1.5 * (q3 - q1)
        return finite & (x > threshold)

    def pct_outlier(x):
        threshold = np.nanquantile(x, percentile_cutoff)
        return finite & (x > threshold)

    if method == "IQR":
        return iqr_outlier(rates)
    if method == "percentile":
        return pct_outlier(rates)
    if method == "both":
        return iqr_outlier(rates) | pct_outlier(rates)
    raise ValueError(f"Unknown outlier method: {method}")


def _detect_sample_outliers(
    sample_status: np.ndarray,
    metric: str = "below_lod",
    method: str = "IQR",
    percentile_cutoff: float = 0.95,
) -> dict:
    """Per-sample outlier flags plus a per-sample log frame."""
    if metric not in ("below_lod", "above_uloq", "either"):
        raise ValueError(f"Unknown metric: {metric}")
    status = np.asarray(sample_status, dtype=object)
    n_metab = status.shape[1]
    if n_metab == 0:
        raise ValueError("sample_status has zero metabolite columns.")

    lod_rate = (status == "< LOD").sum(axis=1) / n_metab
    uloq_rate = (status == "> ULOQ").sum(axis=1) / n_metab

    outlier_lod = _apply_outlier_method(lod_rate, method, percentile_cutoff)
    outlier_uloq = _apply_outlier_method(uloq_rate, method, percentile_cutoff)

    if metric == "below_lod":
        is_outlier = outlier_lod
    elif metric == "above_uloq":
        is_outlier = outlier_uloq
    else:
        is_outlier = outlier_lod | outlier_uloq

    log = pd.DataFrame(
        {
            "lod_rate": lod_rate,
            "uloq_rate": uloq_rate,
            "outlier_below_lod": outlier_lod,
            "outlier_above_uloq": outlier_uloq,
            "is_outlier": is_outlier,
        }
    )
    return {"outliers": np.asarray(is_outlier, dtype=bool), "log": log}

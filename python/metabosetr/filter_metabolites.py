"""Two-step metabolite filter (port of ``filter-metabolites.R``).

DESIGN.md §4.3. CV first (per-kit AND, or pooled), then LOD-rate. CV uses
the *sample* standard deviation (``ddof=1``) to match R ``sd()`` -- this is
the key numerical-fidelity point (numpy defaults to ``ddof=0``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _cv_pct(x: np.ndarray) -> float:
    """CV(%) = 100 * sd / mean; NA when mean is 0 or input is degenerate."""
    x = np.asarray(x, dtype=float)
    finite = x[~np.isnan(x)]
    if finite.size < 2:  # R sd() is NA for length < 2
        return np.nan
    m = finite.mean()
    if m == 0:
        return np.nan
    s = finite.std(ddof=1)  # sample sd, matches R sd()
    return 100 * s / m


def _cv_by_columns(assay: np.ndarray) -> np.ndarray:
    return np.array([_cv_pct(assay[:, j]) for j in range(assay.shape[1])])


def _cv_filter(
    qc_assay: np.ndarray,
    qc_sample_meta,
    cv_threshold: float = 20,
    cv_mode: str = "per_kit",
) -> dict:
    """CV filter on the QC matrix. Returns pass vector + long log frame."""
    if cv_mode not in ("per_kit", "pooled"):
        raise ValueError(f"Unknown cv_mode: {cv_mode}")
    n_metab = qc_assay.shape[1]

    if cv_mode == "pooled":
        cv_values = _cv_by_columns(qc_assay)
        cv_pass = (~np.isnan(cv_values)) & (cv_values <= cv_threshold)
        log = pd.DataFrame({"cv_pooled": cv_values, "cv_pass": cv_pass})
        return {"pass": np.asarray(cv_pass, dtype=bool), "log": log}

    if "kit_id" not in qc_sample_meta.columns:
        raise ValueError("CV per_kit mode requires `kit_id` in QC sample metadata.")
    kit_col = qc_sample_meta["kit_id"]
    kit_ids = sorted(k for k in kit_col.dropna().unique())
    if len(kit_ids) < 1:
        raise ValueError("No kits found in QC sample metadata.")

    cv_by_kit = np.column_stack(
        [
            _cv_by_columns(qc_assay[(kit_col == k).to_numpy(), :])
            if (kit_col == k).any()
            else np.full(n_metab, np.nan)
            for k in kit_ids
        ]
    )
    cv_pass = np.array(
        [
            bool(np.all((~np.isnan(row)) & (row <= cv_threshold)))
            for row in cv_by_kit
        ]
    )
    log_data = {f"cv_{k}": cv_by_kit[:, i] for i, k in enumerate(kit_ids)}
    log_data["cv_pass"] = cv_pass
    log = pd.DataFrame(log_data)
    return {"pass": np.asarray(cv_pass, dtype=bool), "log": log}


def _lod_rate_filter(sample_status: np.ndarray, lod_rate_threshold: float = 0.5) -> dict:
    """Per-metabolite ``<LOD`` rate on the (already-filtered) status matrix."""
    status = np.asarray(sample_status, dtype=object)
    n = status.shape[0]
    if n == 0:
        raise ValueError("Cannot compute LOD rate from a zero-row status matrix.")
    lod_rate = (status == "< LOD").sum(axis=0) / n
    lod_pass = (~np.isnan(lod_rate)) & (lod_rate <= lod_rate_threshold)
    return {"pass": np.asarray(lod_pass, dtype=bool), "lod_rate": lod_rate}

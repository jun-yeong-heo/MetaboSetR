"""Preprocessing pipeline (port of ``preprocess.R``).

DESIGN.md §4.1. Sample-level outlier detection, then the two-step
metabolite filter (CV first, then LOD rate). Raw assays are not mutated;
filter results are recorded as flags. Optional ULOQ winsorization produces
a copied sample assay.
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np

from .classes import ConcentrationData, PreprocessedData, QCData, ThresholdData
from .filter_metabolites import _cv_filter, _lod_rate_filter
from .filter_samples import _detect_sample_outliers


def _winsorize_uloq(conc: ConcentrationData, thresholds: ThresholdData) -> ConcentrationData:
    """Replace ``> ULOQ`` cells with the kit ULOQ where available."""
    assay = conc.assay.copy()
    kit_ids = list(thresholds.kit_ids)
    kit_col = conc.sample_meta["kit_id"].tolist()
    for i in range(assay.shape[0]):
        kit = kit_col[i]
        if kit is None or kit not in kit_ids:
            continue
        kit_idx = kit_ids.index(kit)
        uloq_row = thresholds.uloq[kit_idx, :]
        cells = (conc.status[i, :] == "> ULOQ") & ~np.isnan(uloq_row)
        assay[i, cells] = uloq_row[cells]
    return ConcentrationData(
        assay=assay,
        status=conc.status,
        sample_meta=conc.sample_meta,
        metabolite_names=conc.metabolite_names,
    )


def preprocess(
    conc: ConcentrationData,
    qc: QCData,
    thresholds: Optional[ThresholdData] = None,
    sample_outlier_method: str = "IQR",
    outlier_metric: str = "below_lod",
    percentile_cutoff: float = 0.95,
    cv_threshold: float = 20,
    cv_mode: str = "per_kit",
    lod_rate_threshold: float = 0.5,
    apply_lod_rate_filter: bool = True,
    winsorize_uloq: bool = False,
) -> PreprocessedData:
    """Run the preprocessing pipeline, returning a :class:`PreprocessedData`."""
    if sample_outlier_method not in ("IQR", "percentile", "both"):
        raise ValueError(f"Unknown sample_outlier_method: {sample_outlier_method}")
    if outlier_metric not in ("below_lod", "above_uloq", "either"):
        raise ValueError(f"Unknown outlier_metric: {outlier_metric}")
    if cv_mode not in ("per_kit", "pooled"):
        raise ValueError(f"Unknown cv_mode: {cv_mode}")

    if conc.metabolite_names != qc.metabolite_names:
        raise ValueError("`conc` and `qc` must share the same metabolite_names.")
    if winsorize_uloq and thresholds is None:
        warnings.warn(
            "winsorize_uloq=True but thresholds=None; skipping winsorization.",
            stacklevel=2,
        )
        winsorize_uloq = False

    # 1. Sample-level filter.
    sample_res = _detect_sample_outliers(
        conc.status,
        metric=outlier_metric,
        method=sample_outlier_method,
        percentile_cutoff=percentile_cutoff,
    )

    # 2a. CV filter (always, on raw QC values).
    cv_res = _cv_filter(
        qc.assay, qc.sample_meta, cv_threshold=cv_threshold, cv_mode=cv_mode
    )

    # 2b. LOD-rate filter (conditional, on samples that passed step 1).
    metab_log = cv_res["log"].copy()
    if apply_lod_rate_filter:
        keep = ~sample_res["outliers"]
        if keep.sum() == 0:
            raise ValueError(
                "All samples flagged as outliers; cannot compute LOD-rate "
                "filter. Loosen sample_outlier_method or set "
                "apply_lod_rate_filter=False."
            )
        lod_res = _lod_rate_filter(
            conc.status[keep, :], lod_rate_threshold=lod_rate_threshold
        )
        metabolite_pass = cv_res["pass"] & lod_res["pass"]
        metab_log["lod_rate"] = lod_res["lod_rate"]
        metab_log["lod_pass"] = lod_res["pass"]
    else:
        metabolite_pass = cv_res["pass"]
        metab_log["lod_rate"] = np.nan
        metab_log["lod_pass"] = None
    metab_log["pass"] = metabolite_pass

    # 3. Optional ULOQ winsorization.
    sample = _winsorize_uloq(conc, thresholds) if winsorize_uloq else conc

    return PreprocessedData(
        sample=sample,
        qc=qc,
        thresholds=thresholds,
        sample_outliers=np.asarray(sample_res["outliers"], dtype=bool),
        metabolite_pass=np.asarray(metabolite_pass, dtype=bool),
        sample_filter_log=sample_res["log"],
        metabolite_filter_log=metab_log,
        params={
            "sample_outlier_method": sample_outlier_method,
            "outlier_metric": outlier_metric,
            "percentile_cutoff": percentile_cutoff,
            "cv_threshold": cv_threshold,
            "cv_mode": cv_mode,
            "lod_rate_threshold": lod_rate_threshold,
            "apply_lod_rate_filter": apply_lod_rate_filter,
            "winsorize_uloq": winsorize_uloq,
        },
    )

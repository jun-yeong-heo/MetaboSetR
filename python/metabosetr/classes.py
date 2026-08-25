"""Data containers (ports of the R S4 classes).

Layout convention is samples x metabolites, matching the WebIDQ export and
the R classes (``ConcentrationData.R`` etc.). The statistics layer
transposes to feature x sample when tests run.

Validity checks in ``__post_init__`` mirror the R ``validity=`` functions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class ConcentrationData:
    """Sample concentrations plus the parallel WebIDQ status matrix."""

    assay: np.ndarray  # float, sample x metabolite
    status: np.ndarray  # str/object, same shape as assay
    sample_meta: pd.DataFrame  # one row per sample
    metabolite_names: list[str]

    def __post_init__(self) -> None:
        _check_conc_like(self)

    def __repr__(self) -> str:
        n_s, n_m = self.assay.shape
        cols = ", ".join(map(str, self.sample_meta.columns))
        return f"<ConcentrationData> {n_s} samples x {n_m} metabolites\n  sample_meta columns: {cols}"


@dataclass
class QCData:
    """Pooled-QC concentrations plus status. Same layout as ConcentrationData."""

    assay: np.ndarray  # float, qc_sample x metabolite
    status: np.ndarray
    sample_meta: pd.DataFrame
    metabolite_names: list[str]

    def __post_init__(self) -> None:
        _check_conc_like(self)

    def __repr__(self) -> str:
        n_s, n_m = self.assay.shape
        cols = ", ".join(map(str, self.sample_meta.columns))
        return f"<QCData> {n_s} QC samples x {n_m} metabolites\n  sample_meta columns: {cols}"


@dataclass
class ThresholdData:
    """Optional per-kit LOD / LLOQ / ULOQ holder (kit x metabolite)."""

    lod: np.ndarray
    lloq: np.ndarray
    uloq: np.ndarray
    kit_ids: list[str]
    metabolite_names: list[str]

    def __post_init__(self) -> None:
        for nm in ("lod", "lloq", "uloq"):
            mat = getattr(self, nm)
            if not np.issubdtype(np.asarray(mat).dtype, np.number):
                raise ValueError(f"@{nm} must be numeric")
        if not (self.lod.shape == self.lloq.shape == self.uloq.shape):
            raise ValueError("@lod, @lloq, @uloq must have identical dimensions")
        if len(self.kit_ids) != self.lod.shape[0]:
            raise ValueError("len(kit_ids) must equal nrow(lod)")
        if len(self.metabolite_names) != self.lod.shape[1]:
            raise ValueError("len(metabolite_names) must equal ncol(lod)")

    def __repr__(self) -> str:
        n_k, n_m = self.lod.shape
        return (
            f"<ThresholdData> {n_k} kits x {n_m} metabolites (lod/lloq/uloq)\n"
            f"  kits: {', '.join(self.kit_ids)}"
        )


@dataclass
class PreprocessedData:
    """Output of :func:`preprocess`. Raw assays are retained unmutated."""

    sample: ConcentrationData
    qc: QCData
    thresholds: Optional[ThresholdData]
    sample_outliers: np.ndarray  # bool, length nrow(sample.assay)
    metabolite_pass: np.ndarray  # bool, length ncol(sample.assay)
    sample_filter_log: pd.DataFrame
    metabolite_filter_log: pd.DataFrame
    params: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not (self.thresholds is None or isinstance(self.thresholds, ThresholdData)):
            raise ValueError("thresholds must be None or a ThresholdData")
        n_samp, n_metab = self.sample.assay.shape
        if len(self.sample_outliers) != n_samp:
            raise ValueError("len(sample_outliers) must equal nrow(sample.assay)")
        if len(self.metabolite_pass) != n_metab:
            raise ValueError("len(metabolite_pass) must equal ncol(sample.assay)")
        if len(self.sample_filter_log) != n_samp:
            raise ValueError("nrow(sample_filter_log) must equal nrow(sample.assay)")
        if len(self.metabolite_filter_log) != n_metab:
            raise ValueError("nrow(metabolite_filter_log) must equal ncol(sample.assay)")
        if self.sample.assay.shape[1] != self.qc.assay.shape[1]:
            raise ValueError("sample and qc must share the same metabolite count")

    def __repr__(self) -> str:
        s = self.sample.assay
        q = self.qc.assay
        return (
            "<PreprocessedData>\n"
            f"  sample: {s.shape[0]} x {s.shape[1]}, qc: {q.shape[0]} x {q.shape[1]}\n"
            f"  sample outliers: {int(np.nansum(self.sample_outliers))} / {len(self.sample_outliers)}\n"
            f"  metabolite pass: {int(np.nansum(self.metabolite_pass))} / {len(self.metabolite_pass)}\n"
            f"  thresholds: {'none' if self.thresholds is None else 'attached'}"
        )


def _check_conc_like(obj) -> None:
    assay = np.asarray(obj.assay)
    if not np.issubdtype(assay.dtype, np.number):
        raise ValueError("@assay must be numeric")
    status = np.asarray(obj.status)
    if status.dtype.kind not in ("U", "S", "O"):
        raise ValueError("@status must be character")
    if assay.shape != status.shape:
        raise ValueError("@assay and @status must have identical dimensions")
    if len(obj.sample_meta) != assay.shape[0]:
        raise ValueError("nrow(sample_meta) must equal nrow(assay)")
    if len(obj.metabolite_names) != assay.shape[1]:
        raise ValueError("len(metabolite_names) must equal ncol(assay)")

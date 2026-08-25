"""QC diagnostics (port of ``qc-diagnostics.R``).

DESIGN.md §4.4. Diagnostic only -- not called by :func:`preprocess`. These
are best-effort ports (CLAUDE.md §3): PCA via NumPy SVD (R ``prcomp``),
PERMANOVA via scikit-bio (R ``vegan::adonis2``), histograms via Matplotlib
(R ``ggplot2``). Plot functions return a Matplotlib ``Figure``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .transform import log2_with_pseudocount


def _qc_log_matrix(assay: np.ndarray) -> np.ndarray:
    """log2 (per-metabolite pseudocount), drop non-finite and zero-var cols."""
    mat_log = log2_with_pseudocount(assay.T).T  # back to sample x metabolite
    finite_keep = np.isfinite(mat_log).all(axis=0)
    mat_log = mat_log[:, finite_keep]
    sd_keep = mat_log.std(axis=0, ddof=1) > 0
    return mat_log[:, sd_keep]


def qc_pca_plot(prep, mode: str = "qc_only"):
    """PCA score plot of pooled QC (or QC + samples). Returns a Figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if mode not in ("qc_only", "with_samples"):
        raise ValueError("mode must be 'qc_only' or 'with_samples'.")
    qc_meta = prep.qc.sample_meta.copy()
    qc_meta["type"] = "QC"
    if mode == "with_samples":
        s_meta = prep.sample.sample_meta.copy()
        s_meta["type"] = "Sample"
        common = [c for c in qc_meta.columns if c in s_meta.columns]
        mat = np.vstack([prep.qc.assay, prep.sample.assay])
        meta = pd.concat([qc_meta[common], s_meta[common]], ignore_index=True)
    else:
        mat = prep.qc.assay
        meta = qc_meta

    mat_log = _qc_log_matrix(mat)
    if mat_log.shape[1] < 2:
        raise ValueError(
            "PCA needs at least two metabolite columns with all-finite "
            f"values; only {mat_log.shape[1]} survive filtering."
        )
    # center + scale, then SVD (matches prcomp(center=TRUE, scale.=TRUE)).
    x = (mat_log - mat_log.mean(0)) / mat_log.std(0, ddof=1)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    var = s**2 / np.sum(s**2) * 100

    fig, ax = plt.subplots()
    kit = meta["kit_id"].astype(str).to_numpy()
    for k in pd.unique(kit):
        sel = kit == k
        ax.scatter(scores[sel, 0], scores[sel, 1], label=str(k), s=40)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.legend(title="kit_id")
    return fig


def qc_permanova(prep):
    """Kit-effect PERMANOVA on Euclidean distance over log values.

    Returns the scikit-bio ``permanova`` result (a pandas Series).
    """
    from scipy.spatial.distance import pdist, squareform
    from skbio.stats.distance import DistanceMatrix, permanova

    qc_meta = prep.qc.sample_meta
    kit = qc_meta["kit_id"]
    if kit.dropna().nunique() < 2:
        raise ValueError(
            "PERMANOVA needs at least two distinct kit_id values in QC metadata."
        )
    mat_log = _qc_log_matrix(prep.qc.assay)
    dm = DistanceMatrix(squareform(pdist(mat_log)))
    return permanova(dm, grouping=list(kit.astype(str)))


def qc_cv_plot(prep):
    """Histogram of per-metabolite CV(%) with the threshold line. Figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    log = prep.metabolite_filter_log
    threshold = prep.params["cv_threshold"]
    cv_cols = [c for c in log.columns if c.startswith("cv_") and c != "cv_pass"]
    if not cv_cols:
        raise ValueError("No `cv_*` columns found in metabolite_filter_log.")
    fig, ax = plt.subplots()
    for col in cv_cols:
        ax.hist(log[col].dropna(), bins=30, alpha=0.5, label=col.replace("cv_", ""))
    ax.axvline(threshold, linestyle="--", color="red")
    ax.set_xlabel("CV (%)")
    ax.set_ylabel("Metabolite count")
    ax.set_title(f"Metabolite CV (threshold = {threshold:g}%)")
    ax.legend()
    return fig


def qc_lod_rate_plot(prep, by: str = "sample"):
    """Histogram of per-sample or per-metabolite ``<LOD`` rates. Figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if by not in ("sample", "metabolite"):
        raise ValueError("by must be 'sample' or 'metabolite'.")
    rates = (
        prep.sample_filter_log["lod_rate"]
        if by == "sample"
        else prep.metabolite_filter_log["lod_rate"]
    )
    rates = np.asarray(rates, dtype=float)
    if rates.size == 0 or np.all(np.isnan(rates)):
        extra = " (apply_lod_rate_filter may have been False)" if by == "metabolite" else ""
        raise ValueError(f"No `<LOD` rate data available{extra}.")
    fig, ax = plt.subplots()
    ax.hist(rates[~np.isnan(rates)], bins=20)
    ax.set_xlabel(f"<LOD rate (per {by})")
    ax.set_ylabel("Count")
    if by == "metabolite":
        ax.axvline(prep.params["lod_rate_threshold"], linestyle="--", color="red")
    return fig


def qc_filter_summary(prep) -> pd.DataFrame:
    """Per-axis total / pass / fail counts."""
    n_samp = prep.sample.assay.shape[0]
    n_metab = prep.sample.assay.shape[1]
    return pd.DataFrame(
        {
            "item": ["samples", "metabolites"],
            "total": [n_samp, n_metab],
            "pass": [
                int(np.sum(~prep.sample_outliers)),
                int(np.sum(prep.metabolite_pass)),
            ],
            "fail": [
                int(np.sum(prep.sample_outliers)),
                int(np.sum(~prep.metabolite_pass)),
            ],
        }
    )

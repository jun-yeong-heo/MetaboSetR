"""MetaboSetR preprocessing layer (Python port).

Preprocessing-first port of the R package (see ../docs/MIGRATION.md).
Covers WebIDQ reading, the two-step metabolite filter, sample-outlier
detection, and the PreprocessedData container. Statistics / GSEA / ORA are
not yet ported.
"""

from __future__ import annotations

from .annotate import annotate_metabolites
from .classes import ConcentrationData, PreprocessedData, QCData, ThresholdData
from .pathways import (
    get_pathway_sets,
    get_pathway_sets_meta,
    list_pathway_sets,
    pathway_set_gmt_path,
)
from .preprocess import preprocess
from .qc_diagnostics import (
    qc_cv_plot,
    qc_filter_summary,
    qc_lod_rate_plot,
    qc_pca_plot,
    qc_permanova,
)
from .read_thresholds import read_thresholds
from .read_webidq import read_webidq, read_webidq_qc
from .stats_feature import test_metabolites
from .stats_gsea import gsea_pathway_sets
from .stats_ora import ora_pathway_sets
from .status import (
    is_above_uloq,
    is_below_lloq,
    is_below_lod,
    is_missing,
    is_valid,
)
from .transform import log2_with_pseudocount

__all__ = [
    "ConcentrationData",
    "QCData",
    "ThresholdData",
    "PreprocessedData",
    "read_webidq",
    "read_webidq_qc",
    "read_thresholds",
    "preprocess",
    "log2_with_pseudocount",
    "is_valid",
    "is_below_lloq",
    "is_below_lod",
    "is_above_uloq",
    "is_missing",
    # statistics + pathway layer
    "test_metabolites",
    "ora_pathway_sets",
    "gsea_pathway_sets",
    "annotate_metabolites",
    "list_pathway_sets",
    "get_pathway_sets",
    "get_pathway_sets_meta",
    "pathway_set_gmt_path",
    # QC diagnostics
    "qc_pca_plot",
    "qc_permanova",
    "qc_cv_plot",
    "qc_lod_rate_plot",
    "qc_filter_summary",
]

#' @include ConcentrationData.R QCData.R ThresholdData.R
NULL

#' PreprocessedData S4 class
#'
#' Output of the preprocessing pipeline. Holds the original sample/QC data
#' plus filter outcomes (sample-level outliers, metabolite CV/LOD-rate
#' decisions). The raw `assay` slots inside `@sample` and `@qc` are not
#' mutated; filter results live alongside them. See DESIGN.md §4.1.
#'
#' @slot sample [ConcentrationData]
#' @slot qc [QCData]
#' @slot thresholds [ThresholdData] or `NULL`
#' @slot sample_outliers logical vector, length `nrow(@sample@assay)`
#' @slot metabolite_pass logical vector, length `ncol(@sample@assay)`
#' @slot sample_filter_log data.frame, one row per sample (per-sample LOD /
#'   ULOQ rates and outlier-rule outcomes)
#' @slot metabolite_filter_log data.frame, one row per metabolite (CV per
#'   kit, CV pass, LOD rate, LOD-rate pass)
#' @slot params list of parameters used by `preprocess()`
#' @export
methods::setClass(
  "PreprocessedData",
  slots = c(
    sample = "ConcentrationData",
    qc = "QCData",
    thresholds = "ANY",
    sample_outliers = "logical",
    metabolite_pass = "logical",
    sample_filter_log = "data.frame",
    metabolite_filter_log = "data.frame",
    params = "list"
  ),
  validity = function(object) {
    errors <- character()
    if (!(is.null(object@thresholds) ||
          methods::is(object@thresholds, "ThresholdData"))) {
      errors <- c(errors,
                  "@thresholds must be NULL or a ThresholdData object")
    }
    n_samp <- nrow(object@sample@assay)
    n_metab <- ncol(object@sample@assay)
    if (length(object@sample_outliers) != n_samp) {
      errors <- c(errors,
                  "length(@sample_outliers) must equal nrow(@sample@assay)")
    }
    if (length(object@metabolite_pass) != n_metab) {
      errors <- c(errors,
                  "length(@metabolite_pass) must equal ncol(@sample@assay)")
    }
    if (nrow(object@sample_filter_log) != n_samp) {
      errors <- c(errors,
                  "nrow(@sample_filter_log) must equal nrow(@sample@assay)")
    }
    if (nrow(object@metabolite_filter_log) != n_metab) {
      errors <- c(errors,
                  "nrow(@metabolite_filter_log) must equal ncol(@sample@assay)")
    }
    if (ncol(object@sample@assay) != ncol(object@qc@assay)) {
      errors <- c(errors,
                  "@sample and @qc must share the same metabolite count")
    }
    if (length(errors)) errors else TRUE
  }
)

#' Construct a [PreprocessedData] object
#'
#' @param sample A [ConcentrationData].
#' @param qc A [QCData].
#' @param thresholds A [ThresholdData] or `NULL`.
#' @param sample_outliers Logical vector, length `nrow(@sample@assay)`.
#' @param metabolite_pass Logical vector, length `ncol(@sample@assay)`.
#' @param sample_filter_log Data.frame, one row per sample.
#' @param metabolite_filter_log Data.frame, one row per metabolite.
#' @param params List of parameters used by `preprocess()`.
#' @return A [PreprocessedData] object.
#' @export
PreprocessedData <- function(sample, qc, thresholds = NULL,
                             sample_outliers, metabolite_pass,
                             sample_filter_log, metabolite_filter_log,
                             params = list()) {
  methods::new("PreprocessedData",
               sample = sample, qc = qc, thresholds = thresholds,
               sample_outliers = sample_outliers,
               metabolite_pass = metabolite_pass,
               sample_filter_log = sample_filter_log,
               metabolite_filter_log = metabolite_filter_log,
               params = params)
}

#' @rdname PreprocessedData-class
#' @param object A [PreprocessedData] object.
#' @export
methods::setMethod(
  "show", "PreprocessedData",
  function(object) {
    cat("<PreprocessedData>\n")
    cat(sprintf("  sample: %d × %d, qc: %d × %d\n",
                nrow(object@sample@assay), ncol(object@sample@assay),
                nrow(object@qc@assay), ncol(object@qc@assay)))
    cat(sprintf("  sample outliers: %d / %d\n",
                sum(object@sample_outliers, na.rm = TRUE),
                length(object@sample_outliers)))
    cat(sprintf("  metabolite pass: %d / %d\n",
                sum(object@metabolite_pass, na.rm = TRUE),
                length(object@metabolite_pass)))
    cat(sprintf("  thresholds: %s\n",
                if (is.null(object@thresholds)) "none" else "attached"))
  }
)

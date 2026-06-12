#' QCData S4 class
#'
#' Container for pooled-QC concentration values with the parallel WebIDQ
#' status annotation. Layout matches [ConcentrationData] (qc_sample ×
#' metabolite). See DESIGN.md §2.6.
#'
#' @slot assay numeric matrix (qc_sample × metabolite)
#' @slot status character matrix of the same shape
#' @slot sample_meta data.frame with one row per QC sample (kit_id, pool_id)
#' @slot metabolite_names character vector of canonical short names
#' @export
methods::setClass(
  "QCData",
  slots = c(
    assay = "matrix",
    status = "matrix",
    sample_meta = "data.frame",
    metabolite_names = "character"
  ),
  validity = function(object) {
    errors <- character()
    if (!is.numeric(object@assay)) {
      errors <- c(errors, "@assay must be numeric")
    }
    if (!is.character(object@status)) {
      errors <- c(errors, "@status must be character")
    }
    if (!identical(dim(object@assay), dim(object@status))) {
      errors <- c(errors,
                  "@assay and @status must have identical dimensions")
    }
    if (nrow(object@sample_meta) != nrow(object@assay)) {
      errors <- c(errors, "nrow(@sample_meta) must equal nrow(@assay)")
    }
    if (length(object@metabolite_names) != ncol(object@assay)) {
      errors <- c(errors,
                  "length(@metabolite_names) must equal ncol(@assay)")
    }
    if (length(errors)) errors else TRUE
  }
)

#' Construct a [QCData] object
#'
#' @param assay Numeric matrix (qc_sample × metabolite).
#' @param status Character matrix of the same shape.
#' @param sample_meta Data.frame with one row per QC sample.
#' @param metabolite_names Character vector of canonical short names.
#' @return A [QCData] object.
#' @export
QCData <- function(assay, status, sample_meta, metabolite_names) {
  methods::new("QCData",
               assay = assay,
               status = status,
               sample_meta = sample_meta,
               metabolite_names = metabolite_names)
}

#' @rdname QCData-class
#' @param object A [QCData] object.
#' @export
methods::setMethod(
  "show", "QCData",
  function(object) {
    cat("<QCData>\n")
    cat(sprintf("  %d QC samples × %d metabolites\n",
                nrow(object@assay), ncol(object@assay)))
    cat(sprintf("  sample_meta columns: %s\n",
                paste(names(object@sample_meta), collapse = ", ")))
  }
)

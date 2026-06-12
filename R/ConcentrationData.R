#' ConcentrationData S4 class
#'
#' Container for sample-level concentration values with the parallel WebIDQ
#' status annotation. See DESIGN.md §2.1-2.3 for the input specification.
#'
#' Layout convention is samples × metabolites — matching the WebIDQ export.
#' The downstream statistics layer transposes to the feature × sample
#' convention when tests run.
#'
#' @slot assay numeric matrix (sample × metabolite) of concentrations
#' @slot status character matrix (sample × metabolite) of WebIDQ status
#'   strings, e.g. `"Valid"`, `"< LOD"`. Same shape as `assay`.
#' @slot sample_meta data.frame with one row per sample (kit_id and any
#'   WebIDQ-export metadata columns)
#' @slot metabolite_names character vector of canonical short names, length
#'   `ncol(assay)`
#' @export
methods::setClass(
  "ConcentrationData",
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

#' Construct a [ConcentrationData] object
#'
#' @param assay Numeric matrix (sample × metabolite) of concentrations.
#' @param status Character matrix of the same shape — WebIDQ status strings.
#' @param sample_meta Data.frame with one row per sample.
#' @param metabolite_names Character vector of canonical short names.
#' @return A [ConcentrationData] object.
#' @export
ConcentrationData <- function(assay, status, sample_meta, metabolite_names) {
  methods::new("ConcentrationData",
               assay = assay,
               status = status,
               sample_meta = sample_meta,
               metabolite_names = metabolite_names)
}

#' @rdname ConcentrationData-class
#' @param object A [ConcentrationData] object.
#' @export
methods::setMethod(
  "show", "ConcentrationData",
  function(object) {
    cat("<ConcentrationData>\n")
    cat(sprintf("  %d samples × %d metabolites\n",
                nrow(object@assay), ncol(object@assay)))
    cat(sprintf("  sample_meta columns: %s\n",
                paste(names(object@sample_meta), collapse = ", ")))
  }
)

#' ThresholdData S4 class
#'
#' Optional per-kit LOD / LLOQ / ULOQ holder. See DESIGN.md §2.5. The
#' thresholds are not used to recompute statuses (that information comes
#' from the WebIDQ export); they are attached to the assay row metadata so
#' downstream reports can cite them.
#'
#' @slot lod numeric matrix (kit × metabolite) of LOD values
#' @slot lloq numeric matrix (kit × metabolite) of LLOQ values (NA allowed)
#' @slot uloq numeric matrix (kit × metabolite) of ULOQ values (NA allowed)
#' @slot kit_ids character vector, length `nrow(lod)`
#' @slot metabolite_names character vector, length `ncol(lod)`
#' @export
methods::setClass(
  "ThresholdData",
  slots = c(
    lod = "matrix",
    lloq = "matrix",
    uloq = "matrix",
    kit_ids = "character",
    metabolite_names = "character"
  ),
  validity = function(object) {
    errors <- character()
    mats <- list(lod = object@lod, lloq = object@lloq, uloq = object@uloq)
    for (nm in names(mats)) {
      if (!is.numeric(mats[[nm]])) {
        errors <- c(errors, sprintf("@%s must be numeric", nm))
      }
    }
    if (!(identical(dim(object@lod), dim(object@lloq)) &&
          identical(dim(object@lod), dim(object@uloq)))) {
      errors <- c(errors,
                  "@lod, @lloq, @uloq must have identical dimensions")
    }
    if (length(object@kit_ids) != nrow(object@lod)) {
      errors <- c(errors, "length(@kit_ids) must equal nrow(@lod)")
    }
    if (length(object@metabolite_names) != ncol(object@lod)) {
      errors <- c(errors,
                  "length(@metabolite_names) must equal ncol(@lod)")
    }
    if (length(errors)) errors else TRUE
  }
)

#' Construct a [ThresholdData] object
#'
#' @param lod,lloq,uloq Numeric matrices (kit × metabolite). All must share
#'   the same dimensions. NA values are permitted in `lloq`/`uloq`.
#' @param kit_ids Character vector of kit IDs (e.g. `c("KIT01", "KIT02")`).
#' @param metabolite_names Character vector of canonical short names.
#' @return A [ThresholdData] object.
#' @export
ThresholdData <- function(lod, lloq, uloq, kit_ids, metabolite_names) {
  methods::new("ThresholdData",
               lod = lod, lloq = lloq, uloq = uloq,
               kit_ids = kit_ids, metabolite_names = metabolite_names)
}

#' @rdname ThresholdData-class
#' @param object A [ThresholdData] object.
#' @export
methods::setMethod(
  "show", "ThresholdData",
  function(object) {
    cat("<ThresholdData>\n")
    cat(sprintf("  %d kits × %d metabolites (lod/lloq/uloq)\n",
                nrow(object@lod), ncol(object@lod)))
    cat(sprintf("  kits: %s\n",
                paste(object@kit_ids, collapse = ", ")))
  }
)

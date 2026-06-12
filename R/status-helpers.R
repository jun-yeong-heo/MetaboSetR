#' Status annotation helpers
#'
#' Predicate helpers over Biocrates status annotations. See DESIGN.md §3.2
#' for the canonical status categories. `< threshold` and `< LLOQ` are
#' considered valid for filter-logic purposes (decisions.md #34); the
#' original annotation is preserved on the assay.
#'
#' NA inputs propagate to NA in every predicate except [is_missing()],
#' which returns TRUE on NA (NA is itself a missing observation).
#'
#' @param status A character vector of status annotations.
#' @return A logical vector the same length as `status`.
#' @name status_helpers
NULL

#' @rdname status_helpers
#' @export
is_valid <- function(status) {
  out <- status %in% c("Valid", "< threshold", "< LLOQ")
  out[is.na(status)] <- NA
  out
}

#' @rdname status_helpers
#' @export
is_below_lloq <- function(status) {
  status == "< LLOQ"
}

#' @rdname status_helpers
#' @export
is_below_lod <- function(status) {
  status == "< LOD"
}

#' @rdname status_helpers
#' @export
is_above_uloq <- function(status) {
  status == "> ULOQ"
}

#' @rdname status_helpers
#' @export
is_missing <- function(status) {
  is.na(status) | status == "Missing"
}

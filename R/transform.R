#' log2 transform with an adaptive per-feature pseudocount
#'
#' Applies `log2(x + c)` where the pseudocount `c` is chosen *per row*
#' (i.e. per feature) as half the smallest finite positive value in that
#' row whenever the row contains a zero, and is `0` (no offset) otherwise.
#' This is the same offset rule the package's GSEA / QC / DE layers use
#' internally before any statistical test (DESIGN.md §4.5, decisions.md
#' #25-26).
#'
#' Conventions:
#'   - `x` is treated row-wise: rows are features (metabolites), columns
#'     are samples. Pass `t(...)` if your input has features as columns.
#'   - A vector input is treated as one feature.
#'   - Inf / NaN cells pass through (log2 of Inf is Inf, of NaN is NaN).
#'   - If `pseudocount` is supplied explicitly, the same scalar is added
#'     to every cell (legacy "global" behaviour, retained for callers
#'     that want to override the per-feature default).
#'
#' @param x A numeric matrix (rows = features) or a numeric vector
#'   (one feature).
#' @param pseudocount Optional scalar. If given, used directly for every
#'   cell; otherwise the per-feature half-min rule applies.
#' @return A matrix or vector of the same shape as `x`, log2-transformed.
#' @examples
#' x <- rbind(low  = c(0, 0.1, 0.5, 1),
#'            high = c(100, 200, 0, 400))
#' log2_with_pseudocount(x)
#' @export
log2_with_pseudocount <- function(x, pseudocount = NULL) {
  if (!is.null(pseudocount)) {
    return(log2(x + pseudocount))
  }
  if (!is.matrix(x)) {
    return(.log2_one_row(x))
  }
  out <- x
  for (i in seq_len(nrow(x))) {
    out[i, ] <- .log2_one_row(x[i, ])
  }
  out
}

.log2_one_row <- function(row) {
  has_zero <- any(row == 0, na.rm = TRUE)
  if (!has_zero) return(log2(row))
  finite_pos <- row[is.finite(row) & row > 0]
  pc <- if (length(finite_pos) > 0L) min(finite_pos) / 2 else 0
  log2(row + pc)
}

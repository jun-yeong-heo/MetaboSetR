#' MetaboSetR
#'
#' Biocrates MxP Quant 1000 preprocessing, metabolite-level statistics, and
#' curated pathway-set GSEA / ORA. See `docs/DESIGN.md` for the full design.
#'
#' @keywords internal
"_PACKAGE"

# Internal helpers for the package live here. Each helper stays unexported
# unless promoted to a top-level function in DESIGN.md §7.3.

.not_implemented <- function(what) {
  stop(sprintf("%s is not yet implemented (skeleton stage).", what),
       call. = FALSE)
}

# Resolve the `group` argument used by stats functions against a
# PreprocessedData. Accepts either a single column name from
# `prep@sample@sample_meta` or a vector of length `nrow(prep@sample@assay)`
# (all samples, before outlier exclusion). Always returns a factor.
.resolve_group <- function(group, prep) {
  meta <- prep@sample@sample_meta
  n <- nrow(prep@sample@assay)
  if (is.character(group) && length(group) == 1L) {
    if (!group %in% names(meta)) {
      stop("`group` column '", group, "' not in prep@sample@sample_meta. ",
           "Available: ", paste(names(meta), collapse = ", "),
           call. = FALSE)
    }
    g <- meta[[group]]
  } else if (length(group) == n) {
    g <- group
  } else {
    stop("`group` must be a single sample_meta column name or a vector of ",
         "length nrow(prep@sample@assay).", call. = FALSE)
  }
  as.factor(g)
}

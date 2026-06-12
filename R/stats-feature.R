# Approach A: per-metabolite group testing (DESIGN.md §6.1).

# limma moderated t. Returns a data.frame ordered to match rownames(assay_log).
.test_limma <- function(assay_log, group, padjust) {
  group <- droplevels(as.factor(group))
  if (length(levels(group)) != 2L) {
    stop("limma test requires exactly two non-empty group levels; got ",
         length(levels(group)), ".", call. = FALSE)
  }
  design <- stats::model.matrix(~ group)
  fit <- limma::lmFit(assay_log, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2L, number = Inf, sort.by = "none")
  out <- data.frame(
    metabolite = rownames(tt),
    log2FC     = tt$logFC,
    AveExpr    = tt$AveExpr,
    t          = tt$t,
    P.Value    = tt$P.Value,
    stringsAsFactors = FALSE
  )
  for (m in padjust) {
    out[[paste0("padj_", m)]] <- stats::p.adjust(out$P.Value, method = m)
  }
  out
}

# Wilcoxon rank-sum, two groups only. Per-metabolite log2FC = mean(g2) -
# mean(g1) on the log-transformed values.
.test_wilcoxon <- function(assay_log, group, padjust) {
  group <- droplevels(as.factor(group))
  levels_g <- levels(group)
  if (length(levels_g) != 2L) {
    stop("Wilcoxon test requires exactly two non-empty group levels; got ",
         length(levels_g), ".", call. = FALSE)
  }
  g1 <- group == levels_g[1L]
  g2 <- group == levels_g[2L]
  n <- nrow(assay_log)
  pvals  <- rep(NA_real_, n)
  fcs    <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    x1 <- assay_log[i, g1]; x1f <- x1[is.finite(x1)]
    x2 <- assay_log[i, g2]; x2f <- x2[is.finite(x2)]
    if (length(x1f) < 1L || length(x2f) < 1L) next
    res <- tryCatch(
      suppressWarnings(stats::wilcox.test(x2f, x1f, exact = FALSE)),
      error = function(e) NULL
    )
    pvals[i] <- if (is.null(res)) NA_real_ else res$p.value
    fcs[i]   <- mean(x2f) - mean(x1f)
  }
  out <- data.frame(
    metabolite = rownames(assay_log),
    log2FC     = fcs,
    P.Value    = pvals,
    stringsAsFactors = FALSE
  )
  for (m in padjust) {
    out[[paste0("padj_", m)]] <- stats::p.adjust(out$P.Value, method = m)
  }
  out
}

#' Feature-level group testing of metabolites
#'
#' Approach A from DESIGN.md §6.1: treat each metabolite as an individual
#' feature and run a group comparison via limma moderated t (default) or
#' Wilcoxon rank-sum, with multiple-testing correction. Concentrations are
#' log2-transformed with an adaptive per-feature pseudocount before testing
#' (DESIGN.md §4.5). Outlier samples (`@sample_outliers`) are excluded.
#'
#' Metabolite-to-metabolite correlation is high within a chemical class, so
#' the BH/BY pair is reported by default (decisions.md #29).
#'
#' @param prep A [PreprocessedData] object.
#' @param group Either the name of a column in `prep@sample@sample_meta` or
#'   a vector of length `nrow(prep@sample@assay)` (all samples, before
#'   outlier exclusion).
#' @param method One of `"limma"` (default) or `"wilcoxon"`.
#' @param padjust Character vector of multiple-testing methods to report;
#'   each becomes a `padj_<method>` column. Default `c("BH", "BY")`.
#' @return A data.frame with one row per metabolite.
#' @export
test_metabolites <- function(prep, group,
                             method = c("limma", "wilcoxon"),
                             padjust = c("BH", "BY")) {
  method <- match.arg(method)
  group <- .resolve_group(group, prep)
  keep <- !prep@sample_outliers
  group <- group[keep]
  assay <- prep@sample@assay[keep, , drop = FALSE]
  assay_log <- log2_with_pseudocount(t(assay))
  rownames(assay_log) <- prep@sample@metabolite_names
  if (method == "limma") {
    .test_limma(assay_log, group, padjust)
  } else {
    .test_wilcoxon(assay_log, group, padjust)
  }
}

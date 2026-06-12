#' @include PreprocessedData.R transform.R
NULL

# Build the matrix used by PCA / PERMANOVA: log2-transformed sample × metab,
# with metabolite columns that have any non-finite values dropped (PCA
# cannot tolerate NA/Inf and we'd rather drop columns than impute silently).
.qc_log_matrix <- function(assay) {
  mat_log <- t(log2_with_pseudocount(t(assay)))
  finite_keep <- colSums(!is.finite(mat_log)) == 0L
  mat_log <- mat_log[, finite_keep, drop = FALSE]
  # Drop zero-variance columns: prcomp(scale. = TRUE) divides by SD and
  # would otherwise abort on a constant metabolite (e.g. all-zero in QC).
  sd_keep <- apply(mat_log, 2L, stats::sd) > 0
  mat_log[, sd_keep, drop = FALSE]
}

#' QC PCA score plot
#'
#' Principal-component score plot of the pooled-QC samples (or QC + biological
#' samples when `mode = "with_samples"`). Diagnostic only — not part of the
#' automated pipeline (DESIGN.md §4.4).
#'
#' @param prep A [PreprocessedData] object.
#' @param mode One of `"qc_only"` (default) or `"with_samples"`.
#' @return A `ggplot` object.
#' @export
qc_pca_plot <- function(prep, mode = c("qc_only", "with_samples")) {
  mode <- match.arg(mode)
  qc_assay <- prep@qc@assay
  qc_meta  <- prep@qc@sample_meta
  qc_meta$type <- "QC"

  if (mode == "with_samples") {
    s_assay <- prep@sample@assay
    s_meta  <- prep@sample@sample_meta
    s_meta$type <- "Sample"
    common <- intersect(names(qc_meta), names(s_meta))
    mat  <- rbind(qc_assay, s_assay)
    meta <- rbind(qc_meta[, common, drop = FALSE],
                  s_meta[, common, drop = FALSE])
  } else {
    mat <- qc_assay
    meta <- qc_meta
  }

  mat_log <- .qc_log_matrix(mat)
  if (ncol(mat_log) < 2L) {
    stop("PCA needs at least two metabolite columns with all-finite ",
         "values; only ", ncol(mat_log), " survive filtering.",
         call. = FALSE)
  }
  pca <- stats::prcomp(mat_log, center = TRUE, scale. = TRUE)

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$kit_id <- meta$kit_id
  scores$type   <- meta$type
  pct <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)

  ggplot2::ggplot(scores,
    ggplot2::aes(x = PC1, y = PC2, color = kit_id, shape = type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", pct[1L]),
                  y = sprintf("PC2 (%.1f%%)", pct[2L])) +
    ggplot2::theme_minimal()
}

#' QC PERMANOVA (kit effect)
#'
#' Quantify how much of the QC-sample variation is explained by `kit_id`
#' via `vegan::adonis2` on Euclidean distance over log-transformed values.
#'
#' @param prep A [PreprocessedData] object.
#' @return The `adonis2` result table (a data.frame subclass).
#' @export
qc_permanova <- function(prep) {
  qc_meta <- prep@qc@sample_meta
  if (!"kit_id" %in% names(qc_meta) ||
      length(unique(stats::na.omit(qc_meta$kit_id))) < 2L) {
    stop("PERMANOVA needs at least two distinct kit_id values in QC ",
         "metadata.", call. = FALSE)
  }
  mat_log <- .qc_log_matrix(prep@qc@assay)
  d <- stats::dist(mat_log)
  vegan::adonis2(d ~ kit_id, data = qc_meta)
}

#' QC CV distribution plot
#'
#' Histogram of per-metabolite CV(%) computed during preprocessing, with the
#' filter threshold drawn as a vertical line.
#'
#' @param prep A [PreprocessedData] object.
#' @return A `ggplot` object.
#' @export
qc_cv_plot <- function(prep) {
  log <- prep@metabolite_filter_log
  threshold <- prep@params$cv_threshold
  cv_cols <- setdiff(grep("^cv_", names(log), value = TRUE), "cv_pass")
  if (length(cv_cols) == 0L) {
    stop("No `cv_*` columns found in metabolite_filter_log.",
         call. = FALSE)
  }

  long <- do.call(rbind, lapply(cv_cols, function(col) {
    data.frame(kit = sub("^cv_", "", col),
               cv = log[[col]],
               stringsAsFactors = FALSE)
  }))

  ggplot2::ggplot(long, ggplot2::aes(x = cv, fill = kit)) +
    ggplot2::geom_histogram(bins = 30L, position = "identity",
                            alpha = 0.5) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        color = "red") +
    ggplot2::labs(x = "CV (%)", y = "Metabolite count",
                  title = sprintf("Metabolite CV (threshold = %g%%)",
                                  threshold)) +
    ggplot2::theme_minimal()
}

#' QC LOD-rate distribution plot
#'
#' Histogram of per-sample or per-metabolite `<LOD` rates.
#'
#' @param prep A [PreprocessedData] object.
#' @param by One of `"sample"` (default) or `"metabolite"`.
#' @return A `ggplot` object.
#' @export
qc_lod_rate_plot <- function(prep, by = c("sample", "metabolite")) {
  by <- match.arg(by)
  rates <- if (by == "sample") {
    prep@sample_filter_log$lod_rate
  } else {
    prep@metabolite_filter_log$lod_rate
  }
  if (length(rates) == 0L || all(is.na(rates))) {
    stop("No `<LOD` rate data available", if (by == "metabolite")
         " (apply_lod_rate_filter may have been FALSE)" else "",
         ".", call. = FALSE)
  }
  threshold <- if (by == "metabolite") prep@params$lod_rate_threshold else NA_real_

  p <- ggplot2::ggplot(data.frame(rate = rates),
                       ggplot2::aes(x = rate)) +
    ggplot2::geom_histogram(bins = 20L) +
    ggplot2::labs(x = sprintf("<LOD rate (per %s)", by),
                  y = "Count") +
    ggplot2::theme_minimal()
  if (!is.na(threshold)) {
    p <- p + ggplot2::geom_vline(xintercept = threshold,
                                 linetype = "dashed", color = "red")
  }
  p
}

#' QC filter pass/fail summary
#'
#' Per-axis (samples / metabolites) total / pass / fail counts.
#'
#' @param prep A [PreprocessedData] object.
#' @return A data.frame with columns `item`, `total`, `pass`, `fail`.
#' @export
qc_filter_summary <- function(prep) {
  n_samp  <- nrow(prep@sample@assay)
  n_metab <- ncol(prep@sample@assay)
  data.frame(
    item  = c("samples", "metabolites"),
    total = c(n_samp, n_metab),
    pass  = c(sum(!prep@sample_outliers, na.rm = TRUE),
              sum(prep@metabolite_pass, na.rm = TRUE)),
    fail  = c(sum(prep@sample_outliers, na.rm = TRUE),
              sum(!prep@metabolite_pass, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

# Quiet `R CMD check` "no visible binding" notes for ggplot NSE columns.
utils::globalVariables(c("PC1", "PC2", "kit_id", "type", "cv", "kit", "rate"))

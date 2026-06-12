# Internal: two-step metabolite-level filter (CV first, then LOD rate).
# See DESIGN.md §4.3 and decisions.md #11-14.

# Coefficient of variation in percent. NA-safe; returns NA when the mean
# is zero or the input is degenerate.
.cv_pct <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA_real_)
  s <- stats::sd(x, na.rm = TRUE)
  100 * s / m
}

# CV filter on the QC matrix. Per-kit (default) requires both kits to
# pass; pooled treats all QC samples as one. Returns a logical pass
# vector (length = n_metabolite) and a long log frame with per-kit CVs.
.cv_filter <- function(qc_assay, qc_sample_meta,
                       cv_threshold = 20,
                       cv_mode = c("per_kit", "pooled")) {
  cv_mode <- match.arg(cv_mode)
  n_metab <- ncol(qc_assay)
  if (cv_mode == "pooled") {
    cv_values <- apply(qc_assay, 2L, .cv_pct)
    pass <- !is.na(cv_values) & cv_values <= cv_threshold
    log <- data.frame(
      cv_pooled = cv_values,
      cv_pass = pass,
      stringsAsFactors = FALSE
    )
    return(list(pass = pass, log = log))
  }

  if (!"kit_id" %in% names(qc_sample_meta)) {
    stop("CV per_kit mode requires `kit_id` in QC sample metadata.",
         call. = FALSE)
  }
  kit_ids <- sort(unique(stats::na.omit(qc_sample_meta$kit_id)))
  if (length(kit_ids) < 1L) {
    stop("No kits found in QC sample metadata.", call. = FALSE)
  }
  cv_by_kit <- vapply(kit_ids, function(k) {
    idx <- !is.na(qc_sample_meta$kit_id) & qc_sample_meta$kit_id == k
    if (!any(idx)) {
      return(rep(NA_real_, n_metab))
    }
    apply(qc_assay[idx, , drop = FALSE], 2L, .cv_pct)
  }, numeric(n_metab))
  if (length(kit_ids) == 1L) cv_by_kit <- matrix(cv_by_kit, ncol = 1L)
  pass <- apply(cv_by_kit, 1L, function(row) {
    all(!is.na(row) & row <= cv_threshold)
  })
  log_list <- stats::setNames(
    lapply(seq_along(kit_ids), function(k) cv_by_kit[, k]),
    paste0("cv_", kit_ids)
  )
  log_list$cv_pass <- pass
  log <- as.data.frame(log_list, stringsAsFactors = FALSE,
                       check.names = FALSE)
  list(pass = pass, log = log)
}

# Per-metabolite <LOD rate computed on the supplied (already-filtered)
# sample-status matrix. Pass when rate <= threshold.
.lod_rate_filter <- function(sample_status,
                             lod_rate_threshold = 0.5) {
  n <- nrow(sample_status)
  if (n == 0L) {
    stop("Cannot compute LOD rate from a zero-row status matrix.",
         call. = FALSE)
  }
  lod_rate <- colSums(sample_status == "< LOD", na.rm = TRUE) / n
  pass <- !is.na(lod_rate) & lod_rate <= lod_rate_threshold
  list(pass = pass, lod_rate = lod_rate)
}

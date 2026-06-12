# Internal: sample-level outlier detection from per-sample status rates.
# See DESIGN.md §4.2 and decisions.md #15-17.

# Apply one of three thresholding methods to a per-sample rate vector.
# Returns a logical vector marking outliers.
.apply_outlier_method <- function(rates, method, percentile_cutoff) {
  iqr_outlier <- function(x) {
    q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
    threshold <- unname(q[2L]) + 1.5 * unname(q[2L] - q[1L])
    !is.na(x) & x > threshold
  }
  pct_outlier <- function(x) {
    threshold <- unname(stats::quantile(x, percentile_cutoff, na.rm = TRUE))
    !is.na(x) & x > threshold
  }
  switch(method,
    "IQR"        = iqr_outlier(rates),
    "percentile" = pct_outlier(rates),
    "both"       = iqr_outlier(rates) | pct_outlier(rates),
    stop("Unknown outlier method: ", method, call. = FALSE)
  )
}

# Compute per-sample <LOD / >ULOQ rates and decide which samples are
# outliers. Returns a list with the logical outlier vector and a
# per-sample log frame ready to attach to PreprocessedData.
.detect_sample_outliers <- function(sample_status,
                                    metric = c("below_lod", "above_uloq", "either"),
                                    method = c("IQR", "percentile", "both"),
                                    percentile_cutoff = 0.95) {
  metric <- match.arg(metric)
  method <- match.arg(method)
  n_metab <- ncol(sample_status)
  if (n_metab == 0L) {
    stop("sample_status has zero metabolite columns.", call. = FALSE)
  }
  lod_rate  <- rowSums(sample_status == "< LOD",  na.rm = TRUE) / n_metab
  uloq_rate <- rowSums(sample_status == "> ULOQ", na.rm = TRUE) / n_metab

  outlier_lod  <- .apply_outlier_method(lod_rate, method, percentile_cutoff)
  outlier_uloq <- .apply_outlier_method(uloq_rate, method, percentile_cutoff)

  is_outlier <- switch(metric,
    "below_lod"  = outlier_lod,
    "above_uloq" = outlier_uloq,
    "either"     = outlier_lod | outlier_uloq
  )

  list(
    outliers = is_outlier,
    log = data.frame(
      lod_rate = lod_rate,
      uloq_rate = uloq_rate,
      outlier_below_lod = outlier_lod,
      outlier_above_uloq = outlier_uloq,
      is_outlier = is_outlier,
      stringsAsFactors = FALSE
    )
  )
}

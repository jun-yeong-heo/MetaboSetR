#' @include ConcentrationData.R QCData.R ThresholdData.R PreprocessedData.R
NULL

# Apply > ULOQ winsorization. Replace cells flagged "> ULOQ" with the
# kit-specific ULOQ when one is available; leave NA-ULOQ cells untouched.
.winsorize_uloq <- function(conc, thresholds) {
  for (i in seq_len(nrow(conc@assay))) {
    kit <- conc@sample_meta$kit_id[i]
    if (is.na(kit) || !(kit %in% thresholds@kit_ids)) next
    kit_idx <- match(kit, thresholds@kit_ids)
    uloq_row <- thresholds@uloq[kit_idx, ]
    cells <- which(conc@status[i, ] == "> ULOQ" & !is.na(uloq_row))
    if (length(cells) > 0L) {
      conc@assay[i, cells] <- uloq_row[cells]
    }
  }
  conc
}

#' Run the preprocessing pipeline
#'
#' Combine sample and QC data into a [PreprocessedData] object after running
#' sample-level outlier detection and the two-step metabolite filter
#' (CV first, then LOD rate). Filter results are recorded as flags;
#' the raw `assay` matrices are not removed (DESIGN.md §4 / §5.4).
#'
#' @param conc A [ConcentrationData] for sample data.
#' @param qc A [QCData] for pooled QC.
#' @param thresholds Optional [ThresholdData]. Required if
#'   `winsorize_uloq = TRUE`.
#' @param sample_outlier_method One of `"IQR"` (default), `"percentile"`,
#'   or `"both"`. See DESIGN.md §4.2.
#' @param outlier_metric One of `"below_lod"` (default), `"above_uloq"`,
#'   or `"either"`. Picks which per-sample rate drives outlier flagging.
#' @param percentile_cutoff Quantile for the `"percentile"` method
#'   (default 0.95).
#' @param cv_threshold CV cutoff in percent (default 20).
#' @param cv_mode `"per_kit"` (AND across kits, default) or `"pooled"`.
#' @param lod_rate_threshold Maximum allowed `<LOD` rate per metabolite
#'   (default 0.5).
#' @param apply_lod_rate_filter If `FALSE`, skip the LOD-rate step
#'   (Biocrates-strict CV-only protocol).
#' @param winsorize_uloq Replace `> ULOQ` cells with the kit ULOQ. Requires
#'   `thresholds` and silently skips cells whose ULOQ is `NA`.
#' @return A [PreprocessedData] object.
#' @export
preprocess <- function(conc, qc,
                       thresholds = NULL,
                       sample_outlier_method = c("IQR", "percentile", "both"),
                       outlier_metric = c("below_lod", "above_uloq", "either"),
                       percentile_cutoff = 0.95,
                       cv_threshold = 20,
                       cv_mode = c("per_kit", "pooled"),
                       lod_rate_threshold = 0.5,
                       apply_lod_rate_filter = TRUE,
                       winsorize_uloq = FALSE) {
  sample_outlier_method <- match.arg(sample_outlier_method)
  outlier_metric <- match.arg(outlier_metric)
  cv_mode <- match.arg(cv_mode)

  if (!identical(conc@metabolite_names, qc@metabolite_names)) {
    stop("`conc` and `qc` must share the same `metabolite_names`.",
         call. = FALSE)
  }
  if (winsorize_uloq && is.null(thresholds)) {
    warning("winsorize_uloq = TRUE but thresholds = NULL; skipping ",
            "winsorization.", call. = FALSE)
    winsorize_uloq <- FALSE
  }

  # 1. Sample-level filter.
  sample_res <- .detect_sample_outliers(
    conc@status,
    metric = outlier_metric,
    method = sample_outlier_method,
    percentile_cutoff = percentile_cutoff
  )

  # 2a. CV filter (always, on raw QC values per decisions.md #12).
  cv_res <- .cv_filter(qc@assay, qc@sample_meta,
                       cv_threshold = cv_threshold,
                       cv_mode = cv_mode)

  # 2b. LOD rate filter (conditionally, on samples that passed step 1
  # — outlier samples removed first per decisions.md #18).
  metab_log <- cv_res$log
  if (apply_lod_rate_filter) {
    keep <- !sample_res$outliers
    if (sum(keep) == 0L) {
      stop("All samples flagged as outliers; cannot compute LOD-rate ",
           "filter. Loosen sample_outlier_method or set ",
           "apply_lod_rate_filter = FALSE.", call. = FALSE)
    }
    lod_res <- .lod_rate_filter(conc@status[keep, , drop = FALSE],
                                lod_rate_threshold = lod_rate_threshold)
    metabolite_pass <- cv_res$pass & lod_res$pass
    metab_log$lod_rate <- lod_res$lod_rate
    metab_log$lod_pass <- lod_res$pass
  } else {
    metabolite_pass <- cv_res$pass
    metab_log$lod_rate <- NA_real_
    metab_log$lod_pass <- NA
  }
  metab_log$pass <- metabolite_pass

  # 3. Optional ULOQ winsorization.
  if (winsorize_uloq) {
    conc <- .winsorize_uloq(conc, thresholds)
  }

  PreprocessedData(
    sample = conc,
    qc = qc,
    thresholds = thresholds,
    sample_outliers = unname(sample_res$outliers),
    metabolite_pass = unname(metabolite_pass),
    sample_filter_log = sample_res$log,
    metabolite_filter_log = metab_log,
    params = list(
      sample_outlier_method = sample_outlier_method,
      outlier_metric = outlier_metric,
      percentile_cutoff = percentile_cutoff,
      cv_threshold = cv_threshold,
      cv_mode = cv_mode,
      lod_rate_threshold = lod_rate_threshold,
      apply_lod_rate_filter = apply_lod_rate_filter,
      winsorize_uloq = winsorize_uloq
    )
  )
}

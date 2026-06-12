fxpath <- function(name) {
  system.file("extdata", "synthetic", name, package = "MetaboSetR",
              mustWork = TRUE)
}

read_synth <- function() {
  list(
    cd = read_webidq(fxpath("sample_conc.xlsx"),
                     fxpath("sample_status.xlsx")),
    qc = read_webidq_qc(fxpath("qc_conc.xlsx"),
                        fxpath("qc_status.xlsx")),
    td = read_thresholds(c(fxpath("thresholds_KIT01.xlsx"),
                           fxpath("thresholds_KIT02.xlsx")))
  )
}

# ---- sample outlier helpers -----------------------------------------------

test_that(".detect_sample_outliers flags by <LOD rate via percentile cutoff", {
  status <- matrix("Valid", nrow = 4L, ncol = 10L)
  status[1L, 1:5] <- "< LOD"   # rate 0.5
  status[2L, 1L]  <- "< LOD"   # rate 0.1
  status[3L, 1:8] <- "< LOD"   # rate 0.8 — should be the only outlier
  status[4L, 1:2] <- "< LOD"   # rate 0.2
  res <- MetaboSetR:::.detect_sample_outliers(
    status, metric = "below_lod", method = "percentile",
    percentile_cutoff = 0.75)
  expect_true(res$outliers[3L])
  expect_false(res$outliers[1L])
  expect_equal(res$log$lod_rate, c(0.5, 0.1, 0.8, 0.2))
})

test_that(".detect_sample_outliers either-mode unions both metrics", {
  status <- matrix("Valid", nrow = 3L, ncol = 10L)
  status[1L, 1:9] <- "< LOD"   # high LOD
  status[2L, 1:9] <- "> ULOQ"  # high ULOQ
  status[3L, 1L]  <- "< LOD"   # low LOD, low ULOQ
  res <- MetaboSetR:::.detect_sample_outliers(
    status, metric = "either", method = "percentile",
    percentile_cutoff = 0.66)
  expect_true(res$outliers[1L])
  expect_true(res$outliers[2L])
  expect_false(res$outliers[3L])
})

# ---- CV filter helpers -----------------------------------------------------

test_that(".cv_filter pooled mode rejects high-CV metabolites", {
  qc_assay <- cbind(
    low_cv  = c(1, 1.05, 0.95, 0.98, 1.02, 1.0),    # CV ~ 3%
    high_cv = c(100, 50, 200, 75, 150, 100)          # CV ~ 50%
  )
  qc_meta <- data.frame(kit_id = rep(c("KIT01", "KIT02"), each = 3L),
                        stringsAsFactors = FALSE)
  res <- MetaboSetR:::.cv_filter(qc_assay, qc_meta,
                                      cv_threshold = 20, cv_mode = "pooled")
  expect_true(res$pass[1L])
  expect_false(res$pass[2L])
})

test_that(".cv_filter per_kit requires every kit to pass", {
  qc_assay <- cbind(
    # kit01 CV low, kit02 CV high — overall fail
    one_kit_high = c(1, 1.05, 0.95, 100, 50, 200),
    # both kits low — pass
    both_low = c(1, 1.05, 0.95, 1, 1.05, 0.95)
  )
  qc_meta <- data.frame(kit_id = rep(c("KIT01", "KIT02"), each = 3L),
                        stringsAsFactors = FALSE)
  res <- MetaboSetR:::.cv_filter(qc_assay, qc_meta,
                                      cv_threshold = 20, cv_mode = "per_kit")
  expect_false(res$pass[1L])
  expect_true(res$pass[2L])
  expect_setequal(names(res$log),
                  c("cv_KIT01", "cv_KIT02", "cv_pass"))
})

# ---- LOD rate filter -------------------------------------------------------

test_that(".lod_rate_filter rejects metabolites past the threshold", {
  status <- matrix("Valid", nrow = 10L, ncol = 3L)
  status[, 1L]   <- "Valid"    # 0% LOD
  status[1:6, 2L] <- "< LOD"   # 60% LOD — fail at 0.5
  status[1:3, 3L] <- "< LOD"   # 30% LOD — pass
  res <- MetaboSetR:::.lod_rate_filter(status, 0.5)
  expect_true(res$pass[1L])
  expect_false(res$pass[2L])
  expect_true(res$pass[3L])
  expect_equal(res$lod_rate, c(0, 0.6, 0.3))
})

# ---- preprocess() end-to-end on synthetic data ----------------------------

test_that("preprocess returns a PreprocessedData with expected dimensions", {
  d <- read_synth()
  prep <- preprocess(d$cd, d$qc)
  expect_s4_class(prep, "PreprocessedData")
  expect_length(prep@sample_outliers, 6L)
  expect_length(prep@metabolite_pass, 30L)
  expect_equal(nrow(prep@metabolite_filter_log), 30L)
  expect_equal(nrow(prep@sample_filter_log), 6L)
  expect_true("pass" %in% names(prep@metabolite_filter_log))
})

test_that("preprocess records params for reproducibility", {
  d <- read_synth()
  prep <- preprocess(d$cd, d$qc, cv_threshold = 25,
                     lod_rate_threshold = 0.6)
  expect_equal(prep@params$cv_threshold, 25)
  expect_equal(prep@params$lod_rate_threshold, 0.6)
  expect_equal(prep@params$cv_mode, "per_kit")
  expect_equal(prep@params$apply_lod_rate_filter, TRUE)
})

test_that("apply_lod_rate_filter = FALSE skips the LOD step", {
  d <- read_synth()
  prep <- preprocess(d$cd, d$qc, apply_lod_rate_filter = FALSE)
  expect_true(all(is.na(prep@metabolite_filter_log$lod_rate)))
  expect_true(all(is.na(prep@metabolite_filter_log$lod_pass)))
  # metabolite_pass driven by CV alone in this mode
  expect_equal(prep@metabolite_pass, prep@metabolite_filter_log$cv_pass)
})

test_that("winsorize_uloq with NULL thresholds warns and is disabled", {
  d <- read_synth()
  expect_warning(
    prep <- preprocess(d$cd, d$qc, winsorize_uloq = TRUE),
    "skipping"
  )
  expect_equal(prep@params$winsorize_uloq, FALSE)
})

test_that("winsorize_uloq leaves > ULOQ cells unchanged when ULOQ is NA", {
  d <- read_synth()
  # In synthetic data: row 3 col 25 is "> ULOQ" but KIT01 col 25 ULOQ is NA.
  # winsorize must not touch that cell (no ULOQ to replace with).
  pre_value <- d$cd@assay[3L, 25L]
  prep <- preprocess(d$cd, d$qc, thresholds = d$td,
                     winsorize_uloq = TRUE)
  expect_equal(prep@sample@assay[3L, 25L], pre_value)
})

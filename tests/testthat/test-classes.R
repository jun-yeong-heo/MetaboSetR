# Class definition tests: constructor accepts valid input, validity rejects
# inconsistent dimensions, accessors round-trip metadata. See DESIGN.md §7.2.

# ---- helpers ---------------------------------------------------------------

mk_conc <- function(n_sample = 3L, n_metab = 4L) {
  ConcentrationData(
    assay = matrix(seq_len(n_sample * n_metab) * 1.0,
                   nrow = n_sample, ncol = n_metab),
    status = matrix("Valid", nrow = n_sample, ncol = n_metab),
    sample_meta = data.frame(
      sample_id = paste0("S", seq_len(n_sample)),
      kit_id = rep("KIT01", n_sample)
    ),
    metabolite_names = paste0("M", seq_len(n_metab))
  )
}

mk_qc <- function(n_qc = 6L, n_metab = 4L) {
  QCData(
    assay = matrix(seq_len(n_qc * n_metab) * 1.0,
                   nrow = n_qc, ncol = n_metab),
    status = matrix("Valid", nrow = n_qc, ncol = n_metab),
    sample_meta = data.frame(
      qc_id = paste0("QC", seq_len(n_qc)),
      kit_id = rep(c("KIT01", "KIT02"), each = n_qc / 2L)
    ),
    metabolite_names = paste0("M", seq_len(n_metab))
  )
}

mk_thresh <- function(n_kit = 2L, n_metab = 4L) {
  ThresholdData(
    lod = matrix(0.1, nrow = n_kit, ncol = n_metab),
    lloq = matrix(0.5, nrow = n_kit, ncol = n_metab),
    uloq = matrix(100, nrow = n_kit, ncol = n_metab),
    kit_ids = paste0("KIT", sprintf("%02d", seq_len(n_kit))),
    metabolite_names = paste0("M", seq_len(n_metab))
  )
}

# ---- ConcentrationData -----------------------------------------------------

test_that("ConcentrationData constructor + show work on valid input", {
  cd <- mk_conc()
  expect_s4_class(cd, "ConcentrationData")
  expect_output(show(cd), "ConcentrationData")
  expect_output(show(cd), "3 .* 4 metabolites")
})

test_that("ConcentrationData rejects mismatched assay/status dimensions", {
  expect_error(
    ConcentrationData(
      assay = matrix(1, 3, 4),
      status = matrix("Valid", 3, 5),
      sample_meta = data.frame(s = paste0("S", 1:3)),
      metabolite_names = paste0("M", 1:4)
    ),
    "identical dimensions"
  )
})

test_that("ConcentrationData rejects mismatched sample_meta nrow", {
  expect_error(
    ConcentrationData(
      assay = matrix(1, 3, 4),
      status = matrix("Valid", 3, 4),
      sample_meta = data.frame(s = paste0("S", 1:2)),
      metabolite_names = paste0("M", 1:4)
    ),
    "nrow"
  )
})

test_that("ConcentrationData rejects metabolite_names length mismatch", {
  expect_error(
    ConcentrationData(
      assay = matrix(1, 3, 4),
      status = matrix("Valid", 3, 4),
      sample_meta = data.frame(s = paste0("S", 1:3)),
      metabolite_names = paste0("M", 1:5)
    ),
    "metabolite_names"
  )
})

# ---- QCData ----------------------------------------------------------------

test_that("QCData constructor works on valid input", {
  qc <- mk_qc()
  expect_s4_class(qc, "QCData")
  expect_output(show(qc), "QCData")
})

test_that("QCData rejects shape mismatch", {
  expect_error(
    QCData(assay = matrix(1, 6, 4),
           status = matrix("Valid", 6, 3),
           sample_meta = data.frame(q = paste0("Q", 1:6)),
           metabolite_names = paste0("M", 1:4)),
    "identical dimensions"
  )
})

# ---- ThresholdData ---------------------------------------------------------

test_that("ThresholdData constructor + show work on valid input", {
  td <- mk_thresh()
  expect_s4_class(td, "ThresholdData")
  expect_output(show(td), "ThresholdData")
  expect_output(show(td), "KIT01")
})

test_that("ThresholdData rejects mismatched lod/lloq/uloq dimensions", {
  expect_error(
    ThresholdData(
      lod = matrix(0.1, 2, 4),
      lloq = matrix(0.5, 2, 5),
      uloq = matrix(100, 2, 4),
      kit_ids = c("KIT01", "KIT02"),
      metabolite_names = paste0("M", 1:4)
    ),
    "identical dimensions"
  )
})

# ---- PreprocessedData ------------------------------------------------------

test_that("PreprocessedData constructor accepts NULL thresholds", {
  cd <- mk_conc()
  qc <- mk_qc()
  pp <- PreprocessedData(
    sample = cd, qc = qc, thresholds = NULL,
    sample_outliers = rep(FALSE, 3),
    metabolite_pass = rep(TRUE, 4),
    sample_filter_log = data.frame(lod_rate = rep(0.1, 3)),
    metabolite_filter_log = data.frame(cv = rep(5, 4)),
    params = list(cv_threshold = 20)
  )
  expect_s4_class(pp, "PreprocessedData")
  expect_null(pp@thresholds)
  expect_output(show(pp), "thresholds: none")
})

test_that("PreprocessedData accepts a ThresholdData", {
  pp <- PreprocessedData(
    sample = mk_conc(), qc = mk_qc(), thresholds = mk_thresh(),
    sample_outliers = rep(FALSE, 3),
    metabolite_pass = rep(TRUE, 4),
    sample_filter_log = data.frame(lod_rate = rep(0.1, 3)),
    metabolite_filter_log = data.frame(cv = rep(5, 4))
  )
  expect_s4_class(pp@thresholds, "ThresholdData")
  expect_output(show(pp), "thresholds: attached")
})

test_that("PreprocessedData rejects vector-length mismatch", {
  expect_error(
    PreprocessedData(
      sample = mk_conc(), qc = mk_qc(),
      sample_outliers = rep(FALSE, 4),    # mismatch: 4 vs 3
      metabolite_pass = rep(TRUE, 4),
      sample_filter_log = data.frame(x = rep(0, 3)),
      metabolite_filter_log = data.frame(y = rep(0, 4))
    ),
    "sample_outliers"
  )
})

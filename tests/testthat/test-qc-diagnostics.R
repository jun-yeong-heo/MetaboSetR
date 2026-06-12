fxpath <- function(name) {
  system.file("extdata", "synthetic", name, package = "MetaboSetR",
              mustWork = TRUE)
}

read_synth_prep <- function(...) {
  cd <- read_webidq(fxpath("sample_conc.xlsx"),
                    fxpath("sample_status.xlsx"))
  qc <- read_webidq_qc(fxpath("qc_conc.xlsx"),
                       fxpath("qc_status.xlsx"))
  preprocess(cd, qc, ...)
}

# ---- qc_pca_plot -----------------------------------------------------------

test_that("qc_pca_plot returns a ggplot in qc_only mode", {
  prep <- read_synth_prep()
  p <- qc_pca_plot(prep, mode = "qc_only")
  expect_s3_class(p, "ggplot")
})

test_that("qc_pca_plot returns a ggplot in with_samples mode", {
  prep <- read_synth_prep()
  p <- qc_pca_plot(prep, mode = "with_samples")
  expect_s3_class(p, "ggplot")
})

# ---- qc_permanova ----------------------------------------------------------

test_that("qc_permanova returns the adonis2 table", {
  prep <- read_synth_prep()
  res <- qc_permanova(prep)
  expect_s3_class(res, "anova")
  # adonis2 returns an anova-like data.frame; the kit_id row + Residual + Total
  expect_true(any(rownames(res) == "kit_id" |
                    rownames(res) == "Model"))
})

# ---- qc_cv_plot ------------------------------------------------------------

test_that("qc_cv_plot returns a ggplot for per_kit log", {
  prep <- read_synth_prep()
  p <- qc_cv_plot(prep)
  expect_s3_class(p, "ggplot")
})

test_that("qc_cv_plot also works in pooled CV mode", {
  prep <- read_synth_prep(cv_mode = "pooled")
  p <- qc_cv_plot(prep)
  expect_s3_class(p, "ggplot")
})

# ---- qc_lod_rate_plot ------------------------------------------------------

test_that("qc_lod_rate_plot returns a ggplot per sample/metabolite", {
  prep <- read_synth_prep()
  p1 <- qc_lod_rate_plot(prep, by = "sample")
  p2 <- qc_lod_rate_plot(prep, by = "metabolite")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("qc_lod_rate_plot(by='metabolite') errors when filter was disabled", {
  prep <- read_synth_prep(apply_lod_rate_filter = FALSE)
  expect_error(qc_lod_rate_plot(prep, by = "metabolite"),
               "No `<LOD` rate data")
})

# ---- qc_filter_summary -----------------------------------------------------

test_that("qc_filter_summary tallies samples and metabolites", {
  prep <- read_synth_prep()
  s <- qc_filter_summary(prep)
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s), c("item", "total", "pass", "fail"))
  expect_setequal(s$item, c("samples", "metabolites"))
  # totals should equal pass + fail per row
  expect_equal(s$total, s$pass + s$fail)
  # 6 samples, 30 metabolites in the synthetic fixture
  expect_equal(s$total[s$item == "samples"], 6L)
  expect_equal(s$total[s$item == "metabolites"], 30L)
})

# End-to-end smoke test: every public function in the v0.1 API runs
# successfully on the synthetic fixtures. This is *not* a correctness test
# (handled by per-module tests); it just guards against regression of the
# full pipeline wiring.

fxpath <- function(name) {
  system.file("extdata", "synthetic", name, package = "MetaboSetR",
              mustWork = TRUE)
}

test_that("v0.1 pipeline runs end-to-end on the synthetic fixtures", {
  # ---- Import -------------------------------------------------------------
  cd <- read_webidq(fxpath("sample_conc.xlsx"),
                    fxpath("sample_status.xlsx"))
  qc <- read_webidq_qc(fxpath("qc_conc.xlsx"),
                       fxpath("qc_status.xlsx"))
  td <- read_thresholds(c(fxpath("thresholds_KIT01.xlsx"),
                          fxpath("thresholds_KIT02.xlsx")))
  expect_s4_class(cd, "ConcentrationData")
  expect_s4_class(qc, "QCData")
  expect_s4_class(td, "ThresholdData")

  # ---- Preprocess ---------------------------------------------------------
  # Permissive CV/LOD so noisy synthetic data leaves the panel alive for
  # the statistical layer below.
  prep <- preprocess(cd, qc, thresholds = td,
                     cv_threshold = 200,
                     apply_lod_rate_filter = FALSE)
  expect_s4_class(prep, "PreprocessedData")

  # ---- QC diagnostics -----------------------------------------------------
  expect_s3_class(qc_pca_plot(prep, mode = "qc_only"), "ggplot")
  expect_s3_class(qc_pca_plot(prep, mode = "with_samples"), "ggplot")
  expect_s3_class(qc_cv_plot(prep), "ggplot")
  expect_s3_class(qc_lod_rate_plot(prep, by = "sample"), "ggplot")
  expect_s3_class(qc_filter_summary(prep), "data.frame")
  expect_s3_class(qc_permanova(prep), "anova")

  # ---- Layer 2A: metabolite feature-level testing ------------------------
  feat_limma <- test_metabolites(prep, group = "kit_id", method = "limma")
  feat_wilc  <- test_metabolites(prep, group = "kit_id", method = "wilcoxon")
  expect_s3_class(feat_limma, "data.frame")
  expect_s3_class(feat_wilc, "data.frame")
  expect_equal(nrow(feat_limma), length(prep@sample@metabolite_names))
  expect_equal(nrow(feat_wilc),  length(prep@sample@metabolite_names))

  # ---- Layer 2B: curated pathway-set GSEA --------------------------------
  res_gsea <- suppressWarnings(
    gsea_pathway_sets(prep, group = "kit_id", domain = "reactome",
                      minSize = 1L, maxSize = 2000L)
  )
  expect_s3_class(res_gsea, "data.table")
  expect_s3_class(gsea_to_df(res_gsea), "data.frame")

  # ---- Layer 2C: over-representation analysis ----------------------------
  sig <- feat_limma$metabolite[feat_limma$P.Value < 0.5]
  res_ora <- suppressWarnings(
    ora_pathway_sets(sig, universe = feat_limma$metabolite,
                     domain = "reactome")
  )
  expect_s3_class(res_ora, "data.frame")

  # ---- Annotation ---------------------------------------------------------
  ann <- annotate_metabolites(prep@sample@metabolite_names)
  expect_s3_class(ann, "data.frame")

  # ---- Pathway loaders ---------------------------------------------------
  expect_s3_class(list_pathway_sets(), "data.frame")
  expect_type(get_pathway_sets("reactome"), "list")
  expect_s3_class(get_pathway_sets_meta("reactome"), "data.frame")
  expect_true(file.exists(pathway_set_gmt_path("reactome")))
})

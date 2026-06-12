fxpath <- function(name) {
  system.file("extdata", "synthetic", name, package = "MetaboSetR",
              mustWork = TRUE)
}

# Build a PreprocessedData from synthetic data with a permissive CV
# threshold so the metabolite panel survives QC.
build_prep <- function() {
  cd <- read_webidq(fxpath("sample_conc.xlsx"),
                    fxpath("sample_status.xlsx"))
  qc <- read_webidq_qc(fxpath("qc_conc.xlsx"),
                       fxpath("qc_status.xlsx"))
  preprocess(cd, qc, cv_threshold = 200,
             apply_lod_rate_filter = FALSE)
}

# ---- log2 helper -----------------------------------------------------------

test_that("log2_with_pseudocount applies adaptive pseudocount per row", {
  m <- matrix(c(1, 2, 4,    # no zero -> pure log2
                0, 2, 8),    # has zero -> pseudocount = 1 (min positive 2 / 2)
              nrow = 2L, byrow = TRUE)
  out <- log2_with_pseudocount(m)
  expect_equal(out[1L, ], log2(c(1, 2, 4)))
  expect_equal(out[2L, ], log2(c(0, 2, 8) + 1))
})

test_that("log2_with_pseudocount honors a fixed pseudocount", {
  m <- matrix(c(0, 1, 2), nrow = 1L)
  out <- log2_with_pseudocount(m, pseudocount = 0.5)
  expect_equal(out[1L, ], log2(c(0, 1, 2) + 0.5))
})

# ---- group resolver --------------------------------------------------------

test_that(".resolve_group accepts column name or vector", {
  prep <- build_prep()
  n <- nrow(prep@sample@assay)
  g1 <- MetaboSetR:::.resolve_group("kit_id", prep)
  expect_s3_class(g1, "factor")
  expect_length(g1, n)

  g2 <- MetaboSetR:::.resolve_group(rep(c("A", "B"), length.out = n), prep)
  expect_s3_class(g2, "factor")
  expect_equal(levels(g2), c("A", "B"))
})

test_that(".resolve_group errors on bad input", {
  prep <- build_prep()
  expect_error(MetaboSetR:::.resolve_group("nope", prep),
               "not in prep")
  expect_error(MetaboSetR:::.resolve_group(c("A", "B"), prep),
               "length nrow")
})

# ---- test_metabolites (limma) ---------------------------------------------

test_that("test_metabolites(method='limma') returns one row per metabolite", {
  prep <- build_prep()
  res <- test_metabolites(prep, group = "kit_id", method = "limma")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(prep@sample@metabolite_names))
  expect_setequal(names(res),
                  c("metabolite", "log2FC", "AveExpr", "t",
                    "P.Value", "padj_BH", "padj_BY"))
  expect_true(all(res$metabolite %in% prep@sample@metabolite_names))
})

test_that("test_metabolites custom padjust list controls returned columns", {
  prep <- build_prep()
  res <- test_metabolites(prep, group = "kit_id", method = "limma",
                          padjust = c("BH", "bonferroni"))
  expect_true("padj_BH" %in% names(res))
  expect_true("padj_bonferroni" %in% names(res))
  expect_false("padj_BY" %in% names(res))
})

# ---- test_metabolites (wilcoxon) ------------------------------------------

test_that("test_metabolites(method='wilcoxon') returns expected columns", {
  prep <- build_prep()
  res <- test_metabolites(prep, group = "kit_id", method = "wilcoxon")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(prep@sample@metabolite_names))
  expect_setequal(names(res),
                  c("metabolite", "log2FC", "P.Value",
                    "padj_BH", "padj_BY"))
})

test_that("wilcoxon rejects non-binary groups", {
  prep <- build_prep()
  n <- nrow(prep@sample@assay)
  three_grp <- rep(c("A", "B", "C"), length.out = n)
  expect_error(
    test_metabolites(prep, group = three_grp, method = "wilcoxon"),
    "exactly two .* group levels"
  )
})

test_that("test_metabolites drops unused factor levels", {
  prep <- build_prep()
  n <- nrow(prep@sample@assay)
  fake_group <- factor(rep(c("A", "B"), length.out = n),
                       levels = c("A", "B", "C"))
  expect_no_error(
    test_metabolites(prep, group = fake_group, method = "wilcoxon")
  )
  expect_no_error(
    test_metabolites(prep, group = fake_group, method = "limma")
  )
})

# ---- gsea_to_df (decisions.md #53) ----------------------------------------

test_that("gsea_to_df collapses list columns and is writable to TSV", {
  prep <- build_prep()
  res <- suppressWarnings(
    gsea_pathway_sets(prep, group = "kit_id", domain = "reactome",
                      minSize = 1L, maxSize = 2000L)
  )
  skip_if(nrow(res) == 0L, "fgsea returned no rows on the synthetic data")

  df <- gsea_to_df(res)
  expect_s3_class(df, "data.frame")
  expect_false(inherits(df, "data.table"))
  expect_type(df$leadingEdge, "character")
  expect_equal(nrow(df), nrow(res))
  expect_equal(setdiff(names(res), names(df)), character(0L))

  # Each row's collapsed string must equal its original list joined by `;`.
  expected <- vapply(res$leadingEdge, paste, character(1L), collapse = ";")
  expect_equal(unname(df$leadingEdge), unname(expected))

  # The resulting frame must round-trip through write.table without the
  # `unimplemented type 'list' in 'EncodeElement'` failure that motivated
  # this helper.
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  expect_silent(
    write.table(df, file = tmp, sep = "\t",
                row.names = FALSE, quote = FALSE)
  )
  back <- read.table(tmp, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE)
  expect_equal(nrow(back), nrow(df))
})

test_that("gsea_to_df is a no-op on already-flat frames", {
  flat <- data.frame(pathway = c("A", "B"), pval = c(0.01, 0.5),
                     stringsAsFactors = FALSE)
  out <- gsea_to_df(flat)
  expect_equal(out, flat)
})

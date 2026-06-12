# ---- annotate_metabolites -------------------------------------------------

test_that("annotate_metabolites returns expected layout", {
  out <- annotate_metabolites(c("Nicotine", "TMAO"))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_equal(out$metabolite, c("Nicotine", "TMAO"))
  expect_setequal(
    names(out),
    c("metabolite", "long_name", "class", "kegg_id", "hmdb_id",
      "pubchem_id", "chebi_id", "refmet_id", "lipidmaps_id", "lion_id"))
})

test_that("known metabolites are annotated with their ids", {
  out <- annotate_metabolites("Nicotine")
  expect_equal(out$kegg_id, "C00745")
  expect_equal(out$hmdb_id, "HMDB0001001")
})

test_that("unmatched names return a row of NAs", {
  out <- annotate_metabolites("not_a_metabolite")
  expect_equal(nrow(out), 1L)
  expect_equal(out$metabolite, "not_a_metabolite")
  expect_true(all(is.na(out[, setdiff(names(out), "metabolite")])))
})

test_that("input order is preserved and duplicates are kept", {
  out <- annotate_metabolites(c("TMAO", "Nicotine", "TMAO"))
  expect_equal(out$metabolite, c("TMAO", "Nicotine", "TMAO"))
  expect_equal(out$kegg_id, c("C01104", "C00745", "C01104"))
})

test_that("non-character input is rejected", {
  expect_error(annotate_metabolites(1:3), "character vector")
})

fxpath <- function(name) {
  system.file("extdata", "synthetic", name, package = "MetaboSetR",
              mustWork = TRUE)
}

build_prep <- function() {
  cd <- read_webidq(fxpath("sample_conc.xlsx"),
                    fxpath("sample_status.xlsx"))
  qc <- read_webidq_qc(fxpath("qc_conc.xlsx"),
                       fxpath("qc_status.xlsx"))
  preprocess(cd, qc, cv_threshold = 200,
             apply_lod_rate_filter = FALSE)
}

# ---- list_pathway_sets ----------------------------------------------------

test_that("list_pathway_sets reports the expected layout", {
  L <- list_pathway_sets()
  expect_s3_class(L, "data.frame")
  expect_setequal(names(L), c("domain", "set_name", "n_members"))
  # 5 immunomet + 224 reactome + 145 wikipathways + 266 smpdb
  #  + 109 lion + 2 source + 51 health
  expect_equal(nrow(L), 5L + 224L + 145L + 266L + 109L + 2L + 51L)
  expect_setequal(unique(L$domain),
                  c("immunomet", "reactome", "wikipathways", "smpdb",
                    "lion", "source", "health"))
  expect_true(all(L$n_members >= 1L))
})

# ---- get_pathway_sets -----------------------------------------------------

test_that("get_pathway_sets returns the documented shape per domain", {
  imm_met <- get_pathway_sets("immunomet")
  expect_type(imm_met, "list")
  expect_equal(length(imm_met), 5L)
  expect_true(all(grepl("^IMMUNOMET_ONEIL_", names(imm_met))))

  rx <- get_pathway_sets("reactome")
  expect_equal(length(rx), 224L)
  expect_true(all(grepl("^REACTOME_", names(rx))))
})

test_that("get_pathway_sets errors on unknown domain", {
  expect_error(get_pathway_sets("unknown"), "Unknown domain")
})

# ---- get_pathway_sets_meta ------------------------------------------------

test_that("get_pathway_sets_meta returns long-format metadata", {
  meta_imm <- get_pathway_sets_meta("immunomet")
  expect_setequal(names(meta_imm),
                  c("domain", "set_name", "member", "direction",
                    "brief_note"))
  expect_equal(nrow(meta_imm), 654L)
})

test_that("get_pathway_sets_meta errors on unknown domain", {
  expect_error(get_pathway_sets_meta("unknown"), "Unknown domain")
})

# ---- pathway_set_gmt_path -------------------------------------------------

test_that("pathway_set_gmt_path returns an existing GMT file", {
  for (d in c("immunomet", "reactome", "wikipathways", "smpdb", "lion",
              "source", "health")) {
    p <- pathway_set_gmt_path(d)
    expect_true(file.exists(p), info = d)
    expect_true(grepl("\\.gmt$", p))
  }
})

test_that("installed GMT files parse with fgsea::gmtPathways", {
  # CLAUDE.md §3 pathway-set test minimum.
  for (d in c("immunomet", "reactome", "wikipathways", "smpdb", "lion",
              "source", "health")) {
    parsed <- fgsea::gmtPathways(pathway_set_gmt_path(d))
    expect_type(parsed, "list")
    expect_gt(length(parsed), 0L)
  }
})

# ---- external-DB domain snapshots -----------------------------------------

test_that("external-DB domains expose metabolite-level sets with prefixes", {
  expect_equal(length(get_pathway_sets("reactome")), 224L)
  expect_equal(length(get_pathway_sets("wikipathways")), 145L)
  expect_equal(length(get_pathway_sets("smpdb")), 266L)
  expect_equal(length(get_pathway_sets("lion")), 109L)
  src <- get_pathway_sets("source")
  expect_setequal(names(src), c("SOURCE_PLANT", "SOURCE_MICROBE"))
  expect_true("TMAO" %in% src$SOURCE_MICROBE)
  hlth <- get_pathway_sets("health")
  expect_equal(length(hlth), 51L)
  expect_true(all(grepl("^HEALTH_", names(hlth))))
  expect_true("Ala" %in% hlth$HEALTH_CANCER)

  rx <- get_pathway_sets("reactome")
  # TCA-cycle set carries the expected small molecules (snapshot).
  tca <- rx[["REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"]]
  expect_true(all(c("Citric acid", "Fumaric acid",
                    "a-Ketoglutaric acid") %in% tca))

  expect_true("PA 16:0_18:1" %in%
              get_pathway_sets("lion")[[
                "LION_ABOVE_AVERAGE_BILAYER_THICKNESS"]])
})

# ---- gsea_pathway_sets ----------------------------------------------------

test_that("gsea_pathway_sets runs a RaMP metabolite domain on synthetic data", {
  prep <- build_prep()
  res <- suppressWarnings(
    gsea_pathway_sets(prep, group = "kit_id",
                      domain = "reactome",
                      minSize = 1L, maxSize = 2000L)
  )
  expect_s3_class(res, "data.table")
  expect_true(all(res$pathway %in% names(get_pathway_sets("reactome"))))
})

test_that("gsea_pathway_sets errors when the panel lacks domain coverage", {
  # The synthetic panel has only small molecules, so no lipid overlaps the
  # lion sets — gsea must raise the informative no-overlap error.
  prep <- build_prep()
  expect_error(
    suppressWarnings(
      gsea_pathway_sets(prep, group = "kit_id",
                        domain = "lion",
                        minSize = 1L, maxSize = 2000L)
    ),
    "overlapping"
  )
})

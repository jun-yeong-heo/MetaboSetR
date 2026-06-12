# Verify the internal sysdata.rda objects exist with the expected shapes
# (CLAUDE.md §3 pathway-set test minimums + DESIGN.md §5.2 schema).

test_that("metabolite_dict has 1234 rows and the documented columns", {
  dict <- MetaboSetR:::metabolite_dict
  expect_equal(nrow(dict), 1234L)
  expect_named(dict, c("short_name", "long_name", "class",
                       "kegg_id", "hmdb_id", "pubchem_id", "chebi_id",
                       "refmet_id", "lipidmaps_id", "lion_id"))
  # RaMP/LION id coverage (decisions.md #62, #64, #66) — snapshot the counts.
  expect_equal(sum(!is.na(dict$kegg_id)), 289L)
  expect_equal(sum(!is.na(dict$hmdb_id)), 747L)
  expect_equal(sum(!is.na(dict$lion_id)), 807L)
  expect_equal(sum(!is.na(dict$pubchem_id)), 746L)
  expect_equal(sum(!is.na(dict$chebi_id)), 663L)
  expect_equal(sum(!is.na(dict$refmet_id)), 519L)
  expect_equal(sum(!is.na(dict$lipidmaps_id)), 475L)
  # Spot check known cross-refs: L-Alanine.
  ala <- dict[dict$short_name == "Ala", ]
  expect_equal(ala$kegg_id, "C00041")
  expect_equal(ala$pubchem_id, "5950")
  expect_equal(ala$chebi_id, "CHEBI:16977")
})

test_that("status_color_map covers the five canonical statuses", {
  m <- MetaboSetR:::status_color_map
  expect_named(m, c("hex", "status_text"))
  expect_setequal(m$status_text,
                  c("Valid", "< threshold", "< LOD", "< LLOQ", "> ULOQ"))
})

test_that("pathway_sets has the expected flattened shape", {
  ps <- MetaboSetR:::pathway_sets
  expect_setequal(names(ps),
                  c("immunomet", "reactome", "wikipathways", "smpdb",
                    "lion", "source", "health"))

  expect_equal(length(ps$immunomet), 5L)
  expect_equal(sum(lengths(ps$immunomet)), 654L)
  expect_equal(length(ps$reactome), 224L)
  expect_equal(sum(lengths(ps$reactome)), 2581L)
  expect_equal(length(ps$wikipathways), 145L)
  expect_equal(sum(lengths(ps$wikipathways)), 1305L)
  expect_equal(length(ps$smpdb), 266L)
  expect_equal(sum(lengths(ps$smpdb)), 1216L)
  expect_equal(length(ps$lion), 109L)
  expect_equal(sum(lengths(ps$lion)), 5859L)
  expect_equal(length(ps$source), 2L)
  expect_equal(sum(lengths(ps$source)), 265L)
  expect_equal(length(ps$health), 51L)
  expect_equal(sum(lengths(ps$health)), 933L)
})

test_that("known pathway-set members are present (snapshot)", {
  ps <- MetaboSetR:::pathway_sets
  expect_true(all(lengths(ps$immunomet) > 0L))
  expect_true("Citric acid" %in%
              ps$reactome$REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE)
  expect_true("a-Ketoglutaric acid" %in%
              ps$reactome$REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE)
  expect_true("Gln" %in% ps$smpdb$SMPDB_GLUTAMINOLYSIS_AND_CANCER)
  expect_true("PA 16:0_18:1" %in%
              ps$lion$LION_ABOVE_AVERAGE_BILAYER_THICKNESS)
  # source origin: TMAO is microbe-derived, Trigonelline plant-derived.
  expect_true("TMAO" %in% ps$source$SOURCE_MICROBE)
  expect_true("Trigonelline" %in% ps$source$SOURCE_PLANT)
  # health: disease-association sets (e.g. amino acids in cancer).
  expect_true("Ala" %in% ps$health$HEALTH_CANCER)
})

test_that("pathway_sets_meta has the documented long-format schema", {
  meta <- MetaboSetR:::pathway_sets_meta
  expect_named(meta,
               c("domain", "set_name", "member", "direction", "brief_note"))
  expect_setequal(meta$domain,
                  c("immunomet", "reactome", "wikipathways", "smpdb",
                    "lion", "source", "health"))
  expect_false(any(is.na(meta$set_name)))
  expect_false(any(is.na(meta$member)))
  expect_equal(nrow(meta),
               654L + 2581L + 1305L + 1216L + 5859L + 265L + 933L)
})

test_that("direction values are constrained to the spec set", {
  # CLAUDE.md §2.2 forbids any non-NA direction outside {up, down, altered, -}.
  meta <- MetaboSetR:::pathway_sets_meta
  dir <- meta$direction[!is.na(meta$direction)]
  expect_true(all(dir %in% c("up", "down", "altered", "-")))
  # immunomet metabolite-level sets carry no direction column (§2.2 schema).
  expect_true(all(is.na(meta$direction[meta$domain == "immunomet"])))
  # External-DB domains are hypothesis-free: direction is always "-".
  ext <- meta$direction[meta$domain %in%
                        c("reactome", "wikipathways", "smpdb", "lion",
                          "source", "health")]
  expect_true(all(ext == "-"))
})

test_that("set names follow the domain-prefix naming convention", {
  ps <- MetaboSetR:::pathway_sets
  expect_true(all(grepl("^IMMUNOMET_ONEIL_", names(ps$immunomet))))
  expect_true(all(grepl("^REACTOME_", names(ps$reactome))))
  expect_true(all(grepl("^WIKIPATHWAYS_", names(ps$wikipathways))))
  expect_true(all(grepl("^SMPDB_", names(ps$smpdb))))
  expect_true(all(grepl("^LION_", names(ps$lion))))
  expect_true(all(grepl("^SOURCE_", names(ps$source))))
  expect_true(all(grepl("^HEALTH_", names(ps$health))))
})

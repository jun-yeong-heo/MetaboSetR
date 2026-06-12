# Status helper predicates (DESIGN.md §3.2, decisions.md #34).

test_that("is_valid covers Valid / < threshold / < LLOQ and propagates NA", {
  s <- c("Valid", "< threshold", "< LLOQ", "< LOD", "> ULOQ",
         "Missing", NA_character_)
  expect_equal(is_valid(s),
               c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, NA))
})

test_that("is_below_lloq matches only `< LLOQ`", {
  s <- c("< LLOQ", "Valid", "< LOD", NA_character_)
  expect_equal(is_below_lloq(s), c(TRUE, FALSE, FALSE, NA))
})

test_that("is_below_lod matches only `< LOD`", {
  s <- c("< LOD", "< LLOQ", "Valid", NA_character_)
  expect_equal(is_below_lod(s), c(TRUE, FALSE, FALSE, NA))
})

test_that("is_above_uloq matches only `> ULOQ`", {
  s <- c("> ULOQ", "Valid", "< LOD", NA_character_)
  expect_equal(is_above_uloq(s), c(TRUE, FALSE, FALSE, NA))
})

test_that("is_missing returns TRUE for NA and for the literal `Missing`", {
  s <- c("Missing", NA_character_, "Valid", "< LOD")
  expect_equal(is_missing(s), c(TRUE, TRUE, FALSE, FALSE))
})

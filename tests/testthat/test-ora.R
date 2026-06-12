# ---- ora_pathway_sets -----------------------------------------------------

# Build a small synthetic universe/sig against a known metabolite-level set
# so the hypergeometric arithmetic can be checked by hand.

test_that("ora_pathway_sets returns the documented layout", {
  uni <- unique(unlist(get_pathway_sets("reactome")))
  sig <- get_pathway_sets("reactome")[[1]]
  out <- ora_pathway_sets(sig, uni, "reactome")
  expect_s3_class(out, "data.frame")
  expect_setequal(
    names(out),
    c("set_name", "set_size", "overlap", "expected", "pval",
      "padj_BH", "padj_BY", "overlap_members"))
  expect_true(all(out$overlap <= out$set_size))
  # Ordered by ascending pval.
  expect_false(is.unsorted(out$pval))
})

test_that("hypergeometric p-value matches a hand computation", {
  sets <- get_pathway_sets("reactome")
  uni  <- unique(unlist(sets))
  one  <- names(sets)[1]
  # Make sig = all members of `one` plus two non-members, so the chosen set
  # is over-represented and the arithmetic is non-trivial.
  non_members <- setdiff(uni, sets[[one]])[1:2]
  sig <- c(sets[[one]], non_members)
  out <- ora_pathway_sets(sig, uni, "reactome")

  N <- length(uni)
  M <- length(intersect(sets[[one]], uni))
  n <- length(intersect(sig, uni))
  k <- length(intersect(sig, sets[[one]]))
  expected_p <- stats::phyper(k - 1L, M, N - M, n, lower.tail = FALSE)

  row <- out[out$set_name == one, ]
  expect_equal(row$overlap, k)
  expect_equal(row$set_size, M)
  expect_equal(row$expected, n * M / N)
  expect_equal(row$pval, expected_p)
})

test_that("members outside universe are excluded from set_size", {
  uni <- unique(unlist(get_pathway_sets("reactome")))
  one <- names(get_pathway_sets("reactome"))[1]
  full_members <- get_pathway_sets("reactome")[[one]]
  # Drop one member from the universe; set_size for that set must shrink.
  # (The dropped member is still in `sig`, hence the expected warning.)
  uni_drop <- setdiff(uni, full_members[1])
  out <- suppressWarnings(
    ora_pathway_sets(full_members, uni_drop, "reactome"))
  ss <- out$set_size[out$set_name == one]
  expect_equal(ss, length(full_members) - 1L)
})

test_that("sig members outside universe are dropped with a warning", {
  uni <- unique(unlist(get_pathway_sets("reactome")))
  sig <- c(get_pathway_sets("reactome")[[1]], "not_a_feature")
  expect_warning(ora_pathway_sets(sig, uni, "reactome"),
                 "not in `universe`")
})

test_that("empty universe and bad input error", {
  uni <- unique(unlist(get_pathway_sets("reactome")))
  expect_error(ora_pathway_sets(character(0), character(0), "reactome"),
               "universe.*empty")
  expect_error(ora_pathway_sets(1:3, uni, "reactome"),
               "character vectors")
  expect_error(ora_pathway_sets(character(0), uni, "bogus"),
               "Unknown domain")
})

test_that("multiple metabolite-level domains work", {
  for (d in c("immunomet", "smpdb", "health")) {
    uni <- unique(unlist(get_pathway_sets(d)))
    sig <- get_pathway_sets(d)[[1]]
    out <- ora_pathway_sets(sig, uni, d)
    expect_gt(nrow(out), 0L)
    # The seeded set should be the most enriched (all its members are sig).
    expect_equal(out$overlap[1], out$set_size[1])
  }
})

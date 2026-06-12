# Smoke test: package loads, exported stubs are present and signal not-yet-
# implemented. Intentionally minimal — full unit tests arrive with each
# implementation step (see CLAUDE.md §3, decisions.md #45).

test_that("package loads", {
  expect_true("MetaboSetR" %in% loadedNamespaces() ||
                requireNamespace("MetaboSetR", quietly = TRUE))
})

# All previously-stubbed functions are now implemented; this file is kept
# as a placeholder for future skeleton-level smoke tests.
test_that("package loads cleanly", {
  expect_true("MetaboSetR" %in% loadedNamespaces() ||
                requireNamespace("MetaboSetR", quietly = TRUE))
})

# GSEA layer: enrichment of curated metabolite-level pathway sets against a
# metabolite rank statistic. DESIGN.md §6.2, decisions.md #27/#35.

# Compute per-metabolite log2FC = mean(g2) - mean(g1) over log2-transformed
# concentrations. Returns a named numeric vector keyed by metabolite name.
# Outlier samples are excluded by `keep`.
.metabolite_log2fc <- function(prep, group) {
  keep <- !prep@sample_outliers
  if (sum(keep) < 2L) {
    stop("Fewer than two non-outlier samples available; cannot compute ",
         "metabolite log2FC.", call. = FALSE)
  }
  group <- group[keep]
  assay <- prep@sample@assay[keep, , drop = FALSE]
  metabs <- prep@sample@metabolite_names

  # log2 with adaptive pseudocount per metabolite (transpose so rows = metab).
  metab_log <- log2_with_pseudocount(t(assay))
  rownames(metab_log) <- metabs

  group <- droplevels(as.factor(group))
  if (length(levels(group)) != 2L) {
    stop("GSEA requires exactly two non-empty group levels; got ",
         length(levels(group)), ".", call. = FALSE)
  }
  g1 <- group == levels(group)[1L]
  g2 <- group == levels(group)[2L]

  m1 <- rowMeans(metab_log[, g1, drop = FALSE], na.rm = TRUE)
  m2 <- rowMeans(metab_log[, g2, drop = FALSE], na.rm = TRUE)
  fc <- m2 - m1

  # fgsea cannot handle NA / non-finite ranks; drop those metabolites.
  fc <- fc[is.finite(fc)]
  fc
}

#' GSEA over curated pathway sets
#'
#' Rank statistics from metabolite log2FC (group2 vs group1) are tested for
#' enrichment against the built-in curated metabolite-level pathway sets
#' (DESIGN.md §6.2). The rank statistic is computed on the
#' `PreprocessedData` concentrations; outlier samples are excluded.
#'
#' Direction information stored in `pathway_sets_meta` is *not* passed into
#' fgsea — it is metadata for downstream interpretation only
#' (decisions.md #40).
#'
#' @param prep A [PreprocessedData] object.
#' @param group Grouping variable as in [test_metabolites()].
#' @param domain One of the available pathway-set domains; see
#'   [list_pathway_sets()] (`"immunomet"`, `"reactome"`, `"wikipathways"`,
#'   `"smpdb"`, `"lion"`, `"source"`, `"health"`).
#' @param ... Passed to `fgsea::fgseaMultilevel()` (e.g., `minSize`,
#'   `maxSize`, `nPermSimple`).
#' @return The fgsea data.table result.
#' @export
gsea_pathway_sets <- function(prep, group, domain, ...) {
  .check_domain(domain)
  group <- .resolve_group(group, prep)
  pathways <- get_pathway_sets(domain)
  fc <- .metabolite_log2fc(prep, group)

  # Restrict each pathway to members present in the rank vector so fgsea
  # does not silently drop unknown names with confusing warnings.
  pathways <- lapply(pathways, function(m) intersect(m, names(fc)))
  pathways <- pathways[lengths(pathways) > 0L]
  if (length(pathways) == 0L) {
    stop("No pathway has any member overlapping the rank vector. ",
         "Check that `prep` covers the metabolites of domain '",
         domain, "'.", call. = FALSE)
  }

  fgsea::fgseaMultilevel(pathways = pathways, stats = fc, ...)
}

#' Convert an fgsea result to a flat data.frame for export
#'
#' [gsea_pathway_sets()] returns the raw fgsea `data.table`, whose
#' `leadingEdge` column is a list of character vectors. R's
#' `write.table()` / `write.csv()` cannot encode list columns and fail with
#' `unimplemented type 'list' in 'EncodeElement'`. This helper collapses
#' every list column with `paste(..., collapse = sep)` and returns a plain
#' data.frame so the result can be written to TSV/CSV directly. Non-list
#' columns are passed through unchanged. See DESIGN.md §6.2.3 and
#' decisions.md #53.
#'
#' @param res An fgsea-style result, typically from [gsea_pathway_sets()].
#' @param sep Single-character separator used to join list-column
#'   members. Default `";"`.
#' @return A `data.frame` with the same rows and columns as `res`; any
#'   list column (e.g. `leadingEdge`) becomes character.
#' @examples
#' \dontrun{
#'   res <- gsea_pathway_sets(prep, group = "Response", domain = "reactome")
#'   write.table(gsea_to_df(res), "reactome_gsea.tsv",
#'               sep = "\t", row.names = FALSE, quote = FALSE)
#' }
#' @export
gsea_to_df <- function(res, sep = ";") {
  out <- as.data.frame(res, stringsAsFactors = FALSE)
  is_list_col <- vapply(out, is.list, logical(1L))
  for (j in which(is_list_col)) {
    out[[j]] <- vapply(out[[j]],
                       function(x) paste(x, collapse = sep),
                       character(1L))
  }
  out
}

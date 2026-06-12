# Approach C: over-representation analysis (ORA) over curated pathway sets.
# Threshold-based complement to the rank-based GSEA in stats-gsea.R.
# DESIGN.md §6.6, decisions.md #68.

#' Over-representation analysis over curated pathway sets (Approach C)
#'
#' Tests whether a threshold-selected set of significant metabolites is
#' over-represented in each curated pathway set, using a one-sided
#' hypergeometric test. This is the threshold-based complement to the
#' rank-based [gsea_pathway_sets()]; both run over the same built-in
#' `pathway_sets`. See DESIGN.md §6.6.
#'
#' The background (`universe`) is the *measured panel*: the full set of
#' metabolites you tested. Restricting the background to what was actually
#' measured is the honest choice for a fixed Biocrates panel — a metabolite
#' you never measured can never be significant, so it does not belong in the
#' denominator (decisions.md #68). Each pathway is intersected with
#' `universe` before testing; members outside `universe` are dropped.
#'
#' `sig` and `universe` are metabolite `short_name`s — the same identifier
#' space as the pathway-set members (see [get_pathway_sets()]).
#'
#' Direction information in `pathway_sets_meta` is *not* used by the test;
#' it is interpretation metadata only, mirroring GSEA (decisions.md #40).
#'
#' @param sig Character vector of significant metabolite names.
#' @param universe Character vector of all measured/tested metabolite names
#'   (the background panel). Must be non-empty; `sig` is intersected with it.
#' @param domain One of the available pathway-set domains; see
#'   [list_pathway_sets()].
#' @param padjust Character vector of multiple-testing methods to report;
#'   each becomes a `padj_<method>` column. Default `c("BH", "BY")`.
#' @return A data.frame with one row per pathway that has at least one
#'   member in `universe`, ordered by ascending `pval`. Columns:
#'   `set_name`, `set_size` (members present in `universe`), `overlap`
#'   (significant members in the set), `expected` (overlap expected by
#'   chance), `pval`, one `padj_<method>` per `padjust`, and
#'   `overlap_members` (the overlapping members, `;`-joined).
#' @examples
#' \dontrun{
#'   res <- test_metabolites(prep, group = "treatment")
#'   ora_pathway_sets(
#'     sig      = res$metabolite[res$P.Value < 0.05],
#'     universe = res$metabolite,
#'     domain   = "reactome")
#' }
#' @export
ora_pathway_sets <- function(sig, universe, domain,
                             padjust = c("BH", "BY")) {
  .check_domain(domain)
  if (!is.character(sig) || !is.character(universe)) {
    stop("`sig` and `universe` must be character vectors of metabolite names.",
         call. = FALSE)
  }
  universe <- unique(universe[!is.na(universe)])
  if (length(universe) == 0L) {
    stop("`universe` is empty; supply the measured/tested feature names.",
         call. = FALSE)
  }
  sig <- unique(sig[!is.na(sig)])
  off <- setdiff(sig, universe)
  if (length(off) > 0L) {
    warning(length(off), " of `sig` not in `universe`; dropping them. ",
            "`sig` should be a subset of the measured panel.",
            call. = FALSE)
  }
  sig <- intersect(sig, universe)

  sets <- get_pathway_sets(domain)
  N <- length(universe)
  n <- length(sig)

  rows <- list()
  for (sn in names(sets)) {
    members <- intersect(sets[[sn]], universe)
    M <- length(members)
    if (M == 0L) next                      # nothing testable in this set
    hits <- intersect(sig, members)
    k <- length(hits)
    # P(X >= k) for the hypergeometric draw.
    pval <- stats::phyper(k - 1L, M, N - M, n, lower.tail = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      set_name        = sn,
      set_size        = M,
      overlap         = k,
      expected        = n * M / N,
      pval            = pval,
      overlap_members = paste(hits, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) {
    stop("No pathway in domain '", domain,
         "' has any member in `universe`.", call. = FALSE)
  }
  out <- do.call(rbind, rows)

  padj <- lapply(padjust, function(m) stats::p.adjust(out$pval, method = m))
  names(padj) <- paste0("padj_", padjust)
  # Keep overlap_members last, padj columns just after pval.
  out <- data.frame(
    out[, c("set_name", "set_size", "overlap", "expected", "pval")],
    padj,
    overlap_members = out$overlap_members,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$pval), ]
  rownames(out) <- NULL
  out
}

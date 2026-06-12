#' List installed pathway sets
#'
#' Returns the available `(domain, set_name)` combinations from the built-in
#' `pathway_sets` list. See DESIGN.md §6.4.
#'
#' @return A data.frame with columns `domain`, `set_name`, `n_members`.
#'   One row per set.
#' @export
list_pathway_sets <- function() {
  rows <- list()
  for (dom in names(pathway_sets)) {
    sets <- pathway_sets[[dom]]
    for (sn in names(sets)) {
      rows[[length(rows) + 1L]] <- data.frame(
        domain = dom,
        set_name = sn,
        n_members = length(sets[[sn]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Validate that `domain` is one of the registered domains.
.check_domain <- function(domain) {
  if (!domain %in% names(pathway_sets)) {
    stop("Unknown domain '", domain, "'. Available: ",
         paste(names(pathway_sets), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Get a pathway-set list ready to feed into fgsea
#'
#' Returns a named list of character vectors, one per set, suitable as the
#' `pathways` argument to `fgsea::fgseaMultilevel()`. All sets are
#' metabolite-level (members are metabolite `short_name`s).
#'
#' @param domain One of the available domains; see [list_pathway_sets()].
#'   In-house curated: `"immunomet"`. External-DB derived: `"reactome"`,
#'   `"wikipathways"`, `"smpdb"` (RaMP-DB pathways), `"lion"` (LION),
#'   `"source"` (RaMP origin: Plant/Microbe), `"health"` (RaMP
#'   disease-association ontology).
#' @return A named list of character vectors.
#' @export
get_pathway_sets <- function(domain) {
  .check_domain(domain)
  pathway_sets[[domain]]
}

#' Get long-format pathway-set metadata for a domain
#'
#' Returns the underlying long-format table joining set_name, member,
#' direction, and brief_note for a domain. Useful for downstream
#' interpretation alongside an fgsea result.
#'
#' @param domain Domain name; see [list_pathway_sets()] for the full set
#'   (`"immunomet"`, `"reactome"`, `"wikipathways"`, `"smpdb"`, `"lion"`,
#'   `"source"`, `"health"`).
#' @return A data.frame with columns `domain`, `set_name`, `member`,
#'   `direction`, `brief_note`.
#' @export
get_pathway_sets_meta <- function(domain) {
  if (!domain %in% pathway_sets_meta$domain) {
    stop("Unknown domain '", domain, "'. Available: ",
         paste(unique(pathway_sets_meta$domain), collapse = ", "),
         call. = FALSE)
  }
  rows <- pathway_sets_meta[pathway_sets_meta$domain == domain, ,
                            drop = FALSE]
  rownames(rows) <- NULL
  rows
}

#' Resolve the path of an installed pathway-set GMT file
#'
#' Convenience wrapper around `system.file()`. The naming convention is
#' `<domain>_metabolites.gmt` (e.g., `reactome_metabolites.gmt`).
#'
#' @param domain Domain name.
#' @return A character path; an empty string if the file is not installed.
#' @export
pathway_set_gmt_path <- function(domain) {
  .check_domain(domain)
  fname <- paste0(domain, "_metabolites.gmt")
  system.file("extdata", "pathway_sets", fname,
              package = "MetaboSetR")
}

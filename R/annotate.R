#' Annotate Biocrates metabolites with cross-reference database IDs
#'
#' Maps Biocrates `short_name`s to the external-database identifiers held in
#' the built-in metabolite dictionary (KEGG, HMDB, PubChem, ChEBI, RefMet,
#' LIPID MAPS, LION). One representative id per metabolite per source; see
#' DESIGN.md §9. The result always has one row per input element, in input
#' order, with `NA` for unmatched names or missing ids — ready to feed into
#' KEGG-map overlays (`pathview`) or external pathway tools.
#'
#' This exposes the id annotations the package already curates so callers no
#' longer need to reach into the internal dictionary directly.
#'
#' @param x Character vector of metabolite `short_name`s, e.g. the
#'   `colnames` of a metabolite assay.
#' @return A data.frame with one row per element of `x` and columns
#'   `metabolite` (the input name), `long_name`, `class`, and the id columns
#'   `kegg_id`, `hmdb_id`, `pubchem_id`, `chebi_id`, `refmet_id`,
#'   `lipidmaps_id`, `lion_id`. Cells are `NA` where unknown.
#' @examples
#' annotate_metabolites(c("Nicotine", "TMAO", "not_a_metabolite"))
#' @export
annotate_metabolites <- function(x) {
  if (!is.character(x)) {
    stop("`x` must be a character vector of metabolite short_names.",
         call. = FALSE)
  }
  ann_cols <- c("long_name", "class", "kegg_id", "hmdb_id", "pubchem_id",
                "chebi_id", "refmet_id", "lipidmaps_id", "lion_id")
  idx <- match(x, metabolite_dict$short_name)
  out <- data.frame(metabolite = x, stringsAsFactors = FALSE)
  for (col in ann_cols) {
    out[[col]] <- metabolite_dict[[col]][idx]
  }
  rownames(out) <- NULL
  out
}

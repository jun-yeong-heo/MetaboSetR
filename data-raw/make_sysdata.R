# data-raw/make_sysdata.R
#
# Build R/sysdata.rda from the curated source files. See DESIGN.md §5.2 and
# CLAUDE.md §6 for the contract.
#
# Run from the project root:
#   Rscript data-raw/make_sysdata.R

if (!file.exists("DESCRIPTION")) {
  stop("Run from the project root: `Rscript data-raw/make_sysdata.R`",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

read_tsv <- function(path) {
  read.delim(path, stringsAsFactors = FALSE, na.strings = "",
             quote = "\"", comment.char = "", check.names = FALSE)
}

load_set_dict <- function(master_path, member_col) {
  df <- read_tsv(master_path)
  set_names <- unique(df$set_name)
  out <- lapply(set_names, function(sn) unique(df[df$set_name == sn, member_col]))
  setNames(out, set_names)
}

load_set_long <- function(master_path, member_col, domain) {
  df <- read_tsv(master_path)
  direction <- if ("direction" %in% names(df)) df$direction else NA_character_
  brief_note <- if ("brief_note" %in% names(df)) df$brief_note else NA_character_
  data.frame(
    domain = domain,
    set_name = df$set_name,
    member = df[[member_col]],
    direction = direction,
    brief_note = brief_note,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# metabolite_dict
# ---------------------------------------------------------------------------

met_raw <- read_tsv(file.path("data-raw", "reference",
                              "biocrates_Quant1000_metabolites.tsv"))
metabolite_dict <- data.frame(
  short_name = met_raw$shortname,
  long_name = met_raw$fullname,
  class = met_raw$analyte_class,
  stringsAsFactors = FALSE
)
stopifnot(nrow(metabolite_dict) == 1234L)

# RaMP / LION ID annotations (DESIGN.md §5.2, decisions.md #62, #64).
# Empty cells become NA via read_tsv's na.strings = "".
ramp_ids <- read_tsv(file.path("data-raw", "pathway_sets", "ramp",
                               "ramp_metabolite_ids.tsv"))
lion_ids <- read_tsv(file.path("data-raw", "pathway_sets", "lion",
                               "lion_lipid_ids.tsv"))
mi <- match(metabolite_dict$short_name, ramp_ids$short_name)
li <- match(metabolite_dict$short_name, lion_ids$short_name)
stopifnot(!anyNA(mi), !anyNA(li))
# IDs are identifiers, not numbers — force character so an all-digit
# column (e.g. pubchem CID) is not read back as integer.
metabolite_dict$kegg_id <- as.character(ramp_ids$kegg_id[mi])
metabolite_dict$hmdb_id <- as.character(ramp_ids$hmdb_id[mi])
metabolite_dict$pubchem_id <- as.character(ramp_ids$pubchem_id[mi])
metabolite_dict$chebi_id <- as.character(ramp_ids$chebi_id[mi])
metabolite_dict$refmet_id <- as.character(ramp_ids$refmet_id[mi])
metabolite_dict$lipidmaps_id <- as.character(ramp_ids$lipidmaps_id[mi])
metabolite_dict$lion_id <- as.character(lion_ids$lion_id[li])

# ---------------------------------------------------------------------------
# status_color_map (DESIGN.md §2.4)
# ---------------------------------------------------------------------------

status_color_map <- data.frame(
  hex = c("#B9DE83", "#BBA7B9", "#A28BA3", "#B2D1DC", "#7FB2C5"),
  status_text = c("Valid", "< threshold", "< LOD", "< LLOQ", "> ULOQ"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# pathway_sets (nested: domain -> set_name -> members). All metabolite-level.
# ---------------------------------------------------------------------------

pathway_sets <- list(
  immunomet = load_set_dict(
    file.path("data-raw", "pathway_sets", "immunomet",
              "immunomet_metabolites_master.tsv"),
    "shortname"),
  reactome = load_set_dict(
    file.path("data-raw", "pathway_sets", "ramp", "reactome_master.tsv"),
    "shortname"),
  wikipathways = load_set_dict(
    file.path("data-raw", "pathway_sets", "ramp", "wikipathways_master.tsv"),
    "shortname"),
  smpdb = load_set_dict(
    file.path("data-raw", "pathway_sets", "ramp", "smpdb_master.tsv"),
    "shortname"),
  lion = load_set_dict(
    file.path("data-raw", "pathway_sets", "lion", "lion_master.tsv"),
    "shortname"),
  source = load_set_dict(
    file.path("data-raw", "pathway_sets", "ramp", "source_master.tsv"),
    "shortname"),
  health = load_set_dict(
    file.path("data-raw", "pathway_sets", "ramp", "health_master.tsv"),
    "shortname")
)

# Sanity checks against DESIGN.md §8.5
expected_immunomet <- c(
  IMMUNOMET_ONEIL_GLYCOLYSIS = 3L,
  IMMUNOMET_ONEIL_TCA = 10L,
  IMMUNOMET_ONEIL_FAO = 58L,
  IMMUNOMET_ONEIL_FAS = 445L,
  IMMUNOMET_ONEIL_AA_METABOLISM = 138L
)
stopifnot(identical(
  sort(lengths(pathway_sets$immunomet)),
  sort(expected_immunomet)
))
stopifnot(length(pathway_sets$reactome) == 224L)
stopifnot(length(pathway_sets$wikipathways) == 145L)
stopifnot(length(pathway_sets$smpdb) == 266L)
stopifnot(length(pathway_sets$lion) == 109L)
stopifnot(length(pathway_sets$source) == 2L)
stopifnot(sum(lengths(pathway_sets$source)) == 265L)
stopifnot(length(pathway_sets$health) == 51L)
stopifnot(sum(lengths(pathway_sets$health)) == 933L)

# ---------------------------------------------------------------------------
# pathway_sets_meta (long-format, per-member metadata)
# ---------------------------------------------------------------------------

pathway_sets_meta <- rbind(
  load_set_long(
    file.path("data-raw", "pathway_sets", "immunomet",
              "immunomet_metabolites_master.tsv"),
    "shortname", "immunomet"),
  load_set_long(
    file.path("data-raw", "pathway_sets", "ramp", "reactome_master.tsv"),
    "shortname", "reactome"),
  load_set_long(
    file.path("data-raw", "pathway_sets", "ramp", "wikipathways_master.tsv"),
    "shortname", "wikipathways"),
  load_set_long(
    file.path("data-raw", "pathway_sets", "ramp", "smpdb_master.tsv"),
    "shortname", "smpdb"),
  load_set_long(
    file.path("data-raw", "pathway_sets", "lion", "lion_master.tsv"),
    "shortname", "lion"),
  load_set_long(
    file.path("data-raw", "pathway_sets", "ramp", "source_master.tsv"),
    "shortname", "source"),
  load_set_long(
    file.path("data-raw", "pathway_sets", "ramp", "health_master.tsv"),
    "shortname", "health")
)
stopifnot(nrow(pathway_sets_meta) ==
          654L + 2581L + 1305L + 1216L + 5859L + 265L + 933L)

# ---------------------------------------------------------------------------
# Write R/sysdata.rda
# ---------------------------------------------------------------------------

usethis::use_data(
  metabolite_dict, status_color_map,
  pathway_sets, pathway_sets_meta,
  internal = TRUE, overwrite = TRUE
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

cat(sprintf("metabolite_dict:                    %d rows\n",
            nrow(metabolite_dict)))
cat(sprintf("status_color_map:                   %d rows\n",
            nrow(status_color_map)))
for (d in names(pathway_sets)) {
  cat(sprintf("pathway_sets$%-14s %d sets, %d members\n",
              paste0(d, ":"), length(pathway_sets[[d]]),
              sum(lengths(pathway_sets[[d]]))))
}
cat(sprintf(paste0("metabolite_dict id coverage: kegg=%d hmdb=%d lion=%d ",
                   "pubchem=%d chebi=%d refmet=%d lipidmaps=%d\n"),
            sum(!is.na(metabolite_dict$kegg_id)),
            sum(!is.na(metabolite_dict$hmdb_id)),
            sum(!is.na(metabolite_dict$lion_id)),
            sum(!is.na(metabolite_dict$pubchem_id)),
            sum(!is.na(metabolite_dict$chebi_id)),
            sum(!is.na(metabolite_dict$refmet_id)),
            sum(!is.na(metabolite_dict$lipidmaps_id))))
cat(sprintf("pathway_sets_meta:                  %d rows\n",
            nrow(pathway_sets_meta)))

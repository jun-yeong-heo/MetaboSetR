# Reader for the optional per-kit LOD/LLOQ/ULOQ threshold xlsx files
# (DESIGN.md §2.5). The thresholds are not used to recompute statuses;
# the WebIDQ export already carries authoritative status annotations.

# Pull `KIT\d+` out of a file's basename. Used because the threshold xlsx
# itself does not carry the kit ID anywhere (decisions.md #51).
.kit_id_from_filename <- function(path) {
  base <- basename(path)
  pos <- regexpr("KIT\\d+", base)
  if (is.na(pos) || pos < 1L) {
    stop("Could not find a `KIT\\d+` token in the threshold file name '",
         base, "'. Rename the file to include the kit ID (e.g. ",
         "`thresholds_KIT01.xlsx`).", call. = FALSE)
  }
  regmatches(base, pos)
}

# Read one per-kit threshold xlsx into a list with the four fields needed
# to assemble a ThresholdData. Vendor exports (decisions.md #51) use a
# 4 rows × (1 + N) cols layout — sheet name is not validated:
#   row 1: "Measurement time" + metabolite names
#   row 2: "LOD (calc.)" + LOD values
#   row 3: "ULOQ"        + ULOQ values
#   row 4: "LLOQ"        + LLOQ values
# Kit ID is recovered from the filename, not from any cell.
.read_one_threshold_file <- function(path) {
  raw <- readxl::read_excel(path, col_names = FALSE,
                            .name_repair = "minimal")
  raw <- as.matrix(as.data.frame(raw, stringsAsFactors = FALSE,
                                 check.names = FALSE))
  if (nrow(raw) < 4L) {
    stop("Threshold file '", path, "' has fewer than 4 rows; expected ",
         "Measurement-time / LOD / ULOQ / LLOQ rows.", call. = FALSE)
  }
  kit_id <- .kit_id_from_filename(path)
  metabolite_names <- as.character(raw[1L, -1L, drop = TRUE])
  lod  <- suppressWarnings(as.numeric(raw[2L, -1L, drop = TRUE]))
  uloq <- suppressWarnings(as.numeric(raw[3L, -1L, drop = TRUE]))
  lloq <- suppressWarnings(as.numeric(raw[4L, -1L, drop = TRUE]))
  list(kit_id = kit_id,
       lod = lod, lloq = lloq, uloq = uloq,
       metabolite_names = metabolite_names)
}

#' Read per-kit LOD/LLOQ/ULOQ threshold files
#'
#' Optional metadata reader. See DESIGN.md §2.5. Each input file covers one
#' kit; supply one path per kit and they are stacked by kit. The
#' thresholds are not used to recompute statuses (that information comes
#' from the WebIDQ export); they are attached to the assay row metadata so
#' downstream reports can cite them.
#'
#' @param threshold_files Character vector of paths to per-kit threshold
#'   xlsx files (e.g., `c("lod_kit01.xlsx", "lod_kit02.xlsx")`).
#' @param synonyms Optional named character vector mapping typed metabolite
#'   names in the threshold file to canonical short names.
#' @return A [ThresholdData] object.
#' @export
read_thresholds <- function(threshold_files, synonyms = NULL) {
  if (length(threshold_files) == 0L) {
    stop("threshold_files must contain at least one path.", call. = FALSE)
  }
  per_kit <- lapply(threshold_files, .read_one_threshold_file)

  ref_metabs <- per_kit[[1L]]$metabolite_names
  for (i in seq_along(per_kit)) {
    if (!identical(per_kit[[i]]$metabolite_names, ref_metabs)) {
      stop("Threshold files do not share the same metabolite columns.",
           call. = FALSE)
    }
  }
  resolved <- .resolve_metabolite_names(ref_metabs, synonyms)
  kit_ids <- vapply(per_kit, `[[`, character(1L), "kit_id")

  to_matrix <- function(field) {
    m <- do.call(rbind, lapply(per_kit, `[[`, field))
    rownames(m) <- kit_ids
    colnames(m) <- resolved
    m
  }

  ThresholdData(
    lod = to_matrix("lod"),
    lloq = to_matrix("lloq"),
    uloq = to_matrix("uloq"),
    kit_ids = kit_ids,
    metabolite_names = resolved
  )
}

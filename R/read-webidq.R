# Reader for WebIDQ-style concentration + status xlsx exports.
# DESIGN.md §2.1-2.4 specifies the input format; §7.3 gives the user-facing
# signatures. Helpers are unexported.

# Header columns that must precede the metabolite columns in any WebIDQ
# concentration / status export. Order is preserved as given by Biocrates.
.WEBIDQ_META_COLS <- c(
  "Sample identification", "Sample description",
  "Submission name", "Collection date",
  "Species", "Material"
)

# Read the concentration sheet of a WebIDQ-style xlsx and return the
# split (sample_meta, assay matrix, metabolite-name vector). All metabolite
# columns are coerced to numeric; meta columns stay character.
.read_webidq_conc_sheet <- function(path) {
  # Sheet name is not validated; the reader assumes the input xlsx contains
  # exactly one sheet and reads the first one. See decisions.md #50.
  n_cols <- ncol(readxl::read_excel(path, n_max = 1L))
  col_types <- c(rep("text", length(.WEBIDQ_META_COLS)),
                 rep("numeric", n_cols - length(.WEBIDQ_META_COLS)))
  # Vendor exports include a metabolite-class header row right under the
  # column headers (cells like "Acylcarnitines", "Alkaloids", ...). It
  # provokes hundreds of "Expecting numeric ... got <text>" warnings as
  # readxl coerces those text cells to NA. The row is dropped immediately
  # below by .row_meta_all_na, so the warnings are pure noise.
  raw <- suppressWarnings(readxl::read_excel(path, col_types = col_types))
  raw <- as.data.frame(raw, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(.WEBIDQ_META_COLS, names(raw))
  if (length(missing) > 0L) {
    stop("Concentration file '", path,
         "' is missing required meta columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  # Drop vendor padding rows (decisions.md #51, refined by #52). A row is
  # treated as padding iff `Sample identification` is NA, because vendor
  # exports leave other meta columns (Collection date, Material) NA on
  # legitimate rows but always fill Sample identification on real samples.
  keep <- !.is_padding_row(raw[["Sample identification"]])
  raw <- raw[keep, , drop = FALSE]
  metabolite_names <- setdiff(names(raw), .WEBIDQ_META_COLS)
  sample_meta <- raw[, .WEBIDQ_META_COLS, drop = FALSE]
  rownames(sample_meta) <- NULL
  assay <- as.matrix(raw[, metabolite_names, drop = FALSE])
  storage.mode(assay) <- "double"
  rownames(assay) <- NULL
  list(sample_meta = sample_meta,
       assay = assay,
       metabolite_names = metabolite_names)
}

# TRUE for rows that look like vendor padding (leading metabolite-class
# header row, trailing blank rows where ID columns are unset, etc.). The
# rule is "Sample identification is NA" — see decisions.md #51 / #52.
# Takes the `Sample identification` column directly so both the wide-form
# conc reader and the long-form color reader (decisions.md #55 (c)) can
# share the same predicate.
.is_padding_row <- function(sample_id) {
  is.na(sample_id)
}

# Read a parallel status-text xlsx and return the status matrix aligned to
# the metabolite columns of `conc_data`.
.read_webidq_status_text <- function(path, conc_data) {
  # See note in .read_webidq_conc_sheet — first sheet, name not validated.
  n_cols <- ncol(readxl::read_excel(path, n_max = 1L))
  # Same blanket suppression as .read_webidq_conc_sheet — text mode here
  # rarely warns, but vendor padding cells can still trip readxl edge
  # cases (e.g. mixed empty/text rows) that we already handle by dropping.
  raw <- suppressWarnings(
    readxl::read_excel(path, col_types = rep("text", n_cols))
  )
  raw <- as.data.frame(raw, stringsAsFactors = FALSE, check.names = FALSE)
  meta_missing <- setdiff(.WEBIDQ_META_COLS, names(raw))
  if (length(meta_missing) > 0L) {
    stop("Status file '", path, "' is missing required meta columns: ",
         paste(meta_missing, collapse = ", "), call. = FALSE)
  }
  # Drop the same kind of blank padding rows the conc reader drops, so
  # row counts can match (decisions.md #51 / #52).
  raw <- raw[!.is_padding_row(raw[["Sample identification"]]), , drop = FALSE]
  if (nrow(raw) != nrow(conc_data$sample_meta)) {
    stop("Status file row count (", nrow(raw),
         ") differs from concentration file (",
         nrow(conc_data$sample_meta),
         "); both files were processed with the same blank-row drop rule, ",
         "so this means the two files cover different sample sets.",
         call. = FALSE)
  }
  missing <- setdiff(conc_data$metabolite_names, names(raw))
  if (length(missing) > 0L) {
    stop("Status file is missing metabolite columns: ",
         paste(head(missing, 10L), collapse = ", "),
         if (length(missing) > 10L) ", ..." else "",
         call. = FALSE)
  }
  status <- as.matrix(raw[, conc_data$metabolite_names, drop = FALSE])
  storage.mode(status) <- "character"
  rownames(status) <- NULL
  status
}

# Extract kit ID from each `Submission name` value. DESIGN.md §2.2: regex
# `KIT\d+`. NA submission_name -> NA kit_id; non-matching string -> NA. The
# explicit !is.na(pos) guard avoids R's NA-in-logical-subscript error
# (decisions.md #51).
.extract_kit_id <- function(submission_names) {
  out <- rep(NA_character_, length(submission_names))
  pos <- regexpr("KIT\\d+", submission_names)
  has <- !is.na(pos) & pos > 0L
  if (any(has)) {
    out[has] <- regmatches(
      submission_names[has],
      regexpr("KIT\\d+", submission_names[has])
    )
  }
  out
}

# Resolve user-facing metabolite names against the canonical short_name set.
# Strategy follows DESIGN.md §2.7:
#   pass 1 — exact match against short_name; or fullname -> short_name.
#   pass 2 — whitespace collapse + leading "Total " strip, retry pass 1.
#   pass 3 — apply user-supplied `synonyms` (typed -> canonical) and retry.
# Anything still unresolved raises an error naming the offenders.
.resolve_metabolite_names <- function(input_names, synonyms = NULL) {
  canonical <- metabolite_dict$short_name
  fullname_map <- setNames(metabolite_dict$short_name,
                           metabolite_dict$long_name)

  resolve_one <- function(nm) {
    if (nm %in% canonical) return(nm)
    if (nm %in% names(fullname_map)) return(unname(fullname_map[nm]))
    nm2 <- sub("^Total\\s+", "", gsub("\\s+", " ", trimws(nm)))
    if (nm2 %in% canonical) return(nm2)
    if (nm2 %in% names(fullname_map)) return(unname(fullname_map[nm2]))
    if (!is.null(synonyms) && nm %in% names(synonyms)) {
      mapped <- unname(synonyms[[nm]])
      if (mapped %in% canonical) return(mapped)
      if (mapped %in% names(fullname_map)) {
        return(unname(fullname_map[mapped]))
      }
    }
    NA_character_
  }

  resolved <- vapply(input_names, resolve_one, character(1L))
  if (any(is.na(resolved))) {
    bad <- input_names[is.na(resolved)]
    stop("Could not resolve these metabolite names to canonical short ",
         "names:\n  ", paste(head(bad, 10L), collapse = ", "),
         if (length(bad) > 10L) ", ..." else "",
         "\nProvide them via the `synonyms = c(typed = canonical)` argument.",
         call. = FALSE)
  }
  unname(resolved)
}

# Internal helper shared by read_webidq() and read_webidq_qc().
.read_webidq_into_matrices <- function(concentration_file, status_file,
                                       status_method, synonyms) {
  conc_data <- .read_webidq_conc_sheet(concentration_file)
  status <- if (status_method == "text") {
    if (is.null(status_file)) {
      stop("status_file is required when status_method = 'text'.",
           call. = FALSE)
    }
    .read_webidq_status_text(status_file, conc_data)
  } else {
    .read_webidq_status_color(concentration_file, conc_data)
  }
  conc_data$sample_meta$kit_id <- .extract_kit_id(
    conc_data$sample_meta[["Submission name"]])
  metabolite_names <- .resolve_metabolite_names(
    conc_data$metabolite_names, synonyms)
  list(assay = conc_data$assay,
       status = status,
       sample_meta = conc_data$sample_meta,
       metabolite_names = metabolite_names)
}

#' Read a WebIDQ concentration export
#'
#' Read a Biocrates WebIDQ concentration xlsx export and pair it with the
#' parallel status annotation. Use `status_method = "text"` (the default,
#' recommended) when you also have the WebIDQ status xlsx; use
#' `status_method = "color"` to recover status from the cell fill colors of
#' the concentration file alone (DESIGN.md §2.3-2.4).
#'
#' @param concentration_file Path to the concentration xlsx (sheet
#'   `Conc_raw data`).
#' @param status_file Path to the status text xlsx. Required when
#'   `status_method = "text"`; ignored otherwise.
#' @param status_method One of `"text"` or `"color"`.
#' @param synonyms Optional named character vector mapping typed metabolite
#'   names in the input to canonical short names (DESIGN.md §2.7).
#' @return A [ConcentrationData] object.
#' @export
read_webidq <- function(concentration_file,
                        status_file = NULL,
                        status_method = c("text", "color"),
                        synonyms = NULL) {
  status_method <- match.arg(status_method)
  parts <- .read_webidq_into_matrices(concentration_file, status_file,
                                      status_method, synonyms)
  ConcentrationData(
    assay = parts$assay,
    status = parts$status,
    sample_meta = parts$sample_meta,
    metabolite_names = parts$metabolite_names
  )
}

#' Read a WebIDQ pooled-QC export
#'
#' Same shape and parsing logic as [read_webidq()], but produces a [QCData]
#' object suitable for downstream CV computation.
#'
#' @inheritParams read_webidq
#' @return A [QCData] object.
#' @export
read_webidq_qc <- function(concentration_file,
                           status_file = NULL,
                           status_method = c("text", "color"),
                           synonyms = NULL) {
  status_method <- match.arg(status_method)
  parts <- .read_webidq_into_matrices(concentration_file, status_file,
                                      status_method, synonyms)
  QCData(
    assay = parts$assay,
    status = parts$status,
    sample_meta = parts$sample_meta,
    metabolite_names = parts$metabolite_names
  )
}

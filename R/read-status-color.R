# Color-based status fallback (DESIGN.md §2.4). Used when the user does not
# have the dedicated status xlsx and we need to recover the status matrix
# from cell fill colors of the concentration xlsx itself.

# Map a 6-character hex (no `#`) to the canonical status text using the
# `status_color_map` sysdata. Returns NA for unknown / blank fills.
.hex_to_status <- function(hex) {
  if (is.na(hex) || nchar(hex) != 6L) return(NA_character_)
  ref <- toupper(sub("^#", "", status_color_map$hex))
  m <- match(toupper(hex), ref)
  if (is.na(m)) NA_character_ else status_color_map$status_text[m]
}

# Read the status matrix from cell fills of `path`. The output is aligned
# to `conc_data$metabolite_names`. Cells without a recognised fill default
# to "Valid" (DESIGN.md §2.4: in WebIDQ exports a Valid cell may carry a
# Valid-green fill or no fill at all, and we treat both as Valid).
.read_webidq_status_color <- function(path, conc_data) {
  # Sheet name is not validated; assume one sheet per file and pick the
  # first. See decisions.md #50.
  first_sheet <- readxl::excel_sheets(path)[1L]
  cells <- tidyxl::xlsx_cells(path, sheets = first_sheet)
  fmts <- tidyxl::xlsx_formats(path)
  fill_rgb <- fmts$local$fill$patternFill$fgColor$rgb
  # tidyxl returns ARGB ("FFB9DE83"); strip the alpha prefix to get RGB hex.
  hex_per_format <- substr(fill_rgb, 3L, 8L)

  # Build column index -> column name mapping from header row (row 1).
  header <- cells[cells$row == 1L, ]
  col_to_name <- setNames(header$character, as.character(header$col))

  # The conc reader drops vendor padding rows where `Sample identification`
  # is NA (decisions.md #51 / #52 / #55 (c)). The xlsx-row -> sample-row
  # mapping must honour the same rule via .is_padding_row() so cell-fills
  # line up with sample_meta.
  sid_xlsx_col <- as.integer(names(col_to_name)[
    !is.na(col_to_name) & col_to_name == "Sample identification"
  ])
  if (length(sid_xlsx_col) != 1L) {
    stop("Could not locate `Sample identification` header column in '",
         path, "' for the color-status fallback.", call. = FALSE)
  }
  sid_cells <- cells[cells$row >= 2L & cells$col == sid_xlsx_col, ]
  kept_xlsx_rows <- sid_cells$row[!.is_padding_row(sid_cells$character)]

  n_samples <- nrow(conc_data$sample_meta)
  if (length(kept_xlsx_rows) != n_samples) {
    stop("Color-status reader saw ", length(kept_xlsx_rows),
         " non-padding sample rows but the conc reader retained ",
         n_samples, "; padding-row rules are out of sync.", call. = FALSE)
  }
  xlsx_to_i <- setNames(seq_along(kept_xlsx_rows),
                        as.character(kept_xlsx_rows))

  metabolite_names <- conc_data$metabolite_names
  out <- matrix("Valid",
                nrow = n_samples,
                ncol = length(metabolite_names),
                dimnames = list(NULL, metabolite_names))

  data_rows <- cells[cells$row %in% kept_xlsx_rows, ]
  for (k in seq_len(nrow(data_rows))) {
    col_name <- col_to_name[as.character(data_rows$col[k])]
    if (is.na(col_name) || !(col_name %in% metabolite_names)) next
    fmt_id <- data_rows$local_format_id[k]
    if (is.na(fmt_id) || fmt_id < 1L || fmt_id > length(hex_per_format)) {
      next
    }
    status <- .hex_to_status(hex_per_format[fmt_id])
    if (is.na(status)) next
    i <- xlsx_to_i[as.character(data_rows$row[k])]
    if (is.na(i)) next
    j <- match(col_name, metabolite_names)
    out[i, j] <- status
  }
  out
}

#' Cross-validate text vs color status annotations
#'
#' Compares the primary text-based status (from a status xlsx) against the
#' fallback color-based status (from cell fills of the concentration xlsx).
#' Reports any cells where the two disagree. See DESIGN.md §2.4.
#'
#' @param prep A [PreprocessedData] object that retains both text and color
#'   status sources.
#' @return A data.frame of mismatched cells (sample, metabolite,
#'   text_status, color_status). Empty if no mismatch.
#' @export
validate_status <- function(prep) {
  .not_implemented("validate_status()")
}

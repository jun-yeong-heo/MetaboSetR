# data-raw/make_synthetic.R
#
# Build small synthetic WebIDQ-style fixtures used by the reader unit tests
# and by README/example code. See DESIGN.md §2 for the input specification
# and ROADMAP.md "Synthetic example data 전략".
#
# Output (all deterministic with seed = 42):
#   inst/extdata/synthetic/sample_conc.xlsx       sample × metabolite
#                                                  concentration with cell
#                                                  fill colors encoding the
#                                                  status (fallback path)
#   inst/extdata/synthetic/sample_status.xlsx     parallel status text matrix
#   inst/extdata/synthetic/qc_conc.xlsx           pooled-QC concentration
#   inst/extdata/synthetic/qc_status.xlsx         pooled-QC status text
#   inst/extdata/synthetic/thresholds_KIT01.xlsx  per-kit LOD/LLOQ/ULOQ
#   inst/extdata/synthetic/thresholds_KIT02.xlsx
#
# Run from the project root:
#   Rscript data-raw/make_synthetic.R

if (!file.exists("DESCRIPTION")) {
  stop("Run from the project root: `Rscript data-raw/make_synthetic.R`",
       call. = FALSE)
}

suppressPackageStartupMessages(library(openxlsx))

OUT_DIR <- file.path("inst", "extdata", "synthetic")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Status hex map (DESIGN.md §2.4) — used for cell-fill encoding in
# sample_conc.xlsx / qc_conc.xlsx.
STATUS_HEX <- c(
  "Valid" = "B9DE83",
  "< threshold" = "BBA7B9",
  "< LOD" = "A28BA3",
  "< LLOQ" = "B2D1DC",
  "> ULOQ" = "7FB2C5"
)

# ---------------------------------------------------------------------------
# Pick a small, diverse panel from the real metabolite catalogue
# ---------------------------------------------------------------------------

met <- read.delim(file.path("data-raw", "reference",
                            "biocrates_Quant1000_metabolites.tsv"),
                  stringsAsFactors = FALSE, na.strings = "",
                  quote = "\"", check.names = FALSE)

picked <- c(
  # 20 standard amino acids
  "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "Ala",
  "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
  "Thr", "Trp", "Tyr", "Val",
  # 10 acylcarnitines (one with parens in the name to exercise that case)
  "C0", "C2", "C3", "C4", "C5", "C6 (C4:1-DC)", "C8",
  "C10", "C12", "C16"
)
stopifnot(all(picked %in% met$shortname))
metabolite_names <- picked
n_metab <- length(metabolite_names)

# ---------------------------------------------------------------------------
# Sample / QC metadata — six WebIDQ meta columns + identifiers
# ---------------------------------------------------------------------------

n_kits <- 2L
samples_per_kit <- 3L
qc_per_kit <- 3L
n_samples <- n_kits * samples_per_kit
n_qcs <- n_kits * qc_per_kit

sample_meta <- data.frame(
  `Sample identification` = sprintf("S%03d", seq_len(n_samples)),
  `Sample description` = paste("Sample", seq_len(n_samples)),
  `Submission name` = rep(c("PROJ-X-KIT01", "PROJ-X-KIT02"),
                          each = samples_per_kit),
  `Collection date` = "2026-01-15",
  `Species` = "Human",
  `Material` = "Plasma",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

qc_meta <- data.frame(
  `Sample identification` = sprintf("QC%02d-%s",
                                    rep(seq_len(qc_per_kit), n_kits),
                                    rep(c("KIT01", "KIT02"),
                                        each = qc_per_kit)),
  `Sample description` = "Pooled QC",
  `Submission name` = rep(c("PROJ-X-KIT01", "PROJ-X-KIT02"),
                          each = qc_per_kit),
  `Collection date` = "2026-01-15",
  `Species` = "Human",
  `Material` = "Plasma",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Generate concentration values + status matrices
# ---------------------------------------------------------------------------

set.seed(42L)
base_mu <- runif(n_metab, min = 1, max = 50)

gen_assay <- function(n_row, seed_offset) {
  set.seed(42L + seed_offset)
  m <- matrix(0, nrow = n_row, ncol = n_metab,
              dimnames = list(NULL, metabolite_names))
  for (j in seq_len(n_metab)) {
    m[, j] <- exp(rnorm(n_row, mean = log(base_mu[j]), sd = 0.3))
  }
  m
}

sample_assay <- gen_assay(n_samples, 0L)
qc_assay <- gen_assay(n_qcs, 100L)

sample_status <- matrix("Valid", nrow = n_samples, ncol = n_metab,
                        dimnames = dimnames(sample_assay))
qc_status <- matrix("Valid", nrow = n_qcs, ncol = n_metab,
                    dimnames = dimnames(qc_assay))

# Sprinkle special-status cells so each category appears at least once.
# Deterministic locations.
apply_special <- function(assay, status, specials) {
  for (cell in specials) {
    i <- cell$row; j <- cell$col; s <- cell$status
    status[i, j] <- s
    if (s == "< LOD")        assay[i, j] <- 0.0
    if (s == "< LLOQ")       assay[i, j] <- 0.05
    if (s == "> ULOQ")       assay[i, j] <- 9999
    if (s == "< threshold")  assay[i, j] <- 0.2
  }
  list(assay = assay, status = status)
}

sample_specials <- list(
  list(row = 1L, col = 5L,  status = "< LOD"),
  list(row = 1L, col = 6L,  status = "< LLOQ"),
  list(row = 2L, col = 7L,  status = "< threshold"),
  list(row = 3L, col = 25L, status = "> ULOQ"),
  list(row = 4L, col = 5L,  status = "< LOD"),
  list(row = 5L, col = 8L,  status = "< LLOQ")
)
res <- apply_special(sample_assay, sample_status, sample_specials)
sample_assay <- res$assay; sample_status <- res$status

qc_specials <- list(
  list(row = 1L, col = 25L, status = "< LOD"),
  list(row = 4L, col = 26L, status = "< LOD")
)
res <- apply_special(qc_assay, qc_status, qc_specials)
qc_assay <- res$assay; qc_status <- res$status

# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

write_conc_xlsx <- function(path, assay, sample_meta, status) {
  wb <- createWorkbook()
  sheet <- "Conc_raw data"
  addWorksheet(wb, sheet)
  full <- cbind(sample_meta, assay)
  writeData(wb, sheet, full, startRow = 1L, startCol = 1L, colNames = TRUE)
  styles <- lapply(STATUS_HEX,
                   function(hex) createStyle(fgFill = paste0("#", hex)))
  meta_cols <- ncol(sample_meta)
  for (i in seq_len(nrow(assay))) {
    for (j in seq_len(ncol(assay))) {
      st <- status[i, j]
      if (st %in% names(styles)) {
        addStyle(wb, sheet, style = styles[[st]],
                 rows = i + 1L, cols = meta_cols + j,
                 gridExpand = FALSE, stack = TRUE)
      }
    }
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

write_status_xlsx <- function(path, status, sample_meta) {
  wb <- createWorkbook()
  # Deliberately use a sheet name different from the conc fixture's
  # "Conc_raw data" — vendor exports do this in practice (decisions.md #50)
  # and exercising the divergence here prevents readers from regressing
  # back to a hardcoded sheet name.
  sheet <- "Conc"
  addWorksheet(wb, sheet)
  full <- cbind(sample_meta, status)
  writeData(wb, sheet, full, startRow = 1L, startCol = 1L, colNames = TRUE)
  saveWorkbook(wb, path, overwrite = TRUE)
}

write_threshold_xlsx <- function(path, kit_id, lod, lloq, uloq, mets) {
  # Vendor exports use a 4-row layout with no kit-ID row and no trailing
  # spaces in the labels (decisions.md #51). The kit ID is recovered from
  # the filename, so `kit_id` is unused here but kept as an argument so
  # callers can keep the previous signature (and as a sanity check that
  # the filename matches).
  if (!grepl(kit_id, basename(path), fixed = TRUE)) {
    stop("write_threshold_xlsx: kit_id '", kit_id,
         "' not present in path '", path, "'.", call. = FALSE)
  }
  wb <- createWorkbook()
  sheet <- "Conc"
  addWorksheet(wb, sheet)
  writeData(wb, sheet, t(c("Measurement time", mets)),
            startRow = 1L, startCol = 1L, colNames = FALSE)
  writeData(wb, sheet, "LOD (calc.)",
            startRow = 2L, startCol = 1L, colNames = FALSE)
  writeData(wb, sheet, t(lod),
            startRow = 2L, startCol = 2L, colNames = FALSE)
  writeData(wb, sheet, "ULOQ",
            startRow = 3L, startCol = 1L, colNames = FALSE)
  writeData(wb, sheet, t(uloq),
            startRow = 3L, startCol = 2L, colNames = FALSE)
  writeData(wb, sheet, "LLOQ",
            startRow = 4L, startCol = 1L, colNames = FALSE)
  writeData(wb, sheet, t(lloq),
            startRow = 4L, startCol = 2L, colNames = FALSE)
  saveWorkbook(wb, path, overwrite = TRUE)
}

# ---------------------------------------------------------------------------
# Generate threshold values
# ---------------------------------------------------------------------------

set.seed(42L)
lod_kit01  <- runif(n_metab, min = 0.05, max = 0.5)
lloq_kit01 <- lod_kit01 * 2
uloq_kit01 <- runif(n_metab, min = 500, max = 2000)
# One NA each in lloq/uloq to mimic the lipid-metabolite gap noted in
# DESIGN.md §2.5 / ROADMAP.md open question 2.
lloq_kit01[25L] <- NA_real_
uloq_kit01[25L] <- NA_real_

lod_kit02  <- lod_kit01  * 1.05
lloq_kit02 <- lloq_kit01 * 1.05
uloq_kit02 <- uloq_kit01 * 0.95

# ---------------------------------------------------------------------------
# Write everything
# ---------------------------------------------------------------------------

write_conc_xlsx(file.path(OUT_DIR, "sample_conc.xlsx"),
                sample_assay, sample_meta, sample_status)
write_status_xlsx(file.path(OUT_DIR, "sample_status.xlsx"),
                  sample_status, sample_meta)
write_conc_xlsx(file.path(OUT_DIR, "qc_conc.xlsx"),
                qc_assay, qc_meta, qc_status)
write_status_xlsx(file.path(OUT_DIR, "qc_status.xlsx"),
                  qc_status, qc_meta)
write_threshold_xlsx(file.path(OUT_DIR, "thresholds_KIT01.xlsx"),
                     "KIT01", lod_kit01, lloq_kit01, uloq_kit01,
                     metabolite_names)
write_threshold_xlsx(file.path(OUT_DIR, "thresholds_KIT02.xlsx"),
                     "KIT02", lod_kit02, lloq_kit02, uloq_kit02,
                     metabolite_names)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

cat(sprintf("Sample concentration: %d × %d  -> sample_conc.xlsx\n",
            n_samples, n_metab))
cat(sprintf("Sample status text:   %d × %d  -> sample_status.xlsx\n",
            n_samples, n_metab))
cat(sprintf("QC concentration:     %d × %d  -> qc_conc.xlsx\n",
            n_qcs, n_metab))
cat(sprintf("QC status text:       %d × %d  -> qc_status.xlsx\n",
            n_qcs, n_metab))
cat(sprintf("Threshold (KIT01):    %d × %d  -> thresholds_KIT01.xlsx\n",
            4L, 1L + n_metab))
cat(sprintf("Threshold (KIT02):    %d × %d  -> thresholds_KIT02.xlsx\n",
            4L, 1L + n_metab))
cat("\nSample status distribution:\n")
print(table(sample_status))
cat("\nQC status distribution:\n")
print(table(qc_status))

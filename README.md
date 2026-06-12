# MetaboSetR

Biocrates MxP Quant 1000 (WebIDQ) preprocessing, metabolite-level
statistics, and curated pathway-set enrichment (GSEA / ORA). R package.

**Status**: v0.1.0. See `docs/ROADMAP.md` for the current status and
`docs/DESIGN.md` for the design rationale.

---

## What it does

Given a Biocrates Q1000 WebIDQ export (sample concentration xlsx + status
xlsx, optional pooled-QC and per-kit threshold files), MetaboSetR:

1. **Imports** WebIDQ exports into typed S4 containers; cell fill color
   serves as a status fallback when the dedicated status xlsx is absent.
2. **Preprocesses** with sample-level outlier detection (IQR / percentile),
   per-kit (or pooled) CV filtering, and `<LOD`-rate filtering.
3. **Tests** metabolites with limma moderated t or Wilcoxon (BH + BY
   adjustment) via `test_metabolites()`.
4. **Enriches** curated metabolite-level pathway sets two ways:
   - **GSEA** — `gsea_pathway_sets()` ranks metabolites by group log2FC and
     runs `fgsea::fgseaMultilevel()`.
   - **ORA** — `ora_pathway_sets()` runs a one-sided hypergeometric test of
     a threshold-selected significant metabolite set against the same sets,
     with the measured panel as the background.

   Built-in domains: 5 immunomet sets (O'Neill 2016), plus external-DB
   metabolite sets mapped via RaMP-DB — 224 reactome, 145 wikipathways,
   266 smpdb, 109 lion lipid-functional sets, a 2-set `source` domain
   (Plant/Microbe origin), and 51 `health` disease-association sets.
5. **Annotates** metabolites with cross-reference ids via
   `annotate_metabolites()` (`kegg_id`/`hmdb_id`/`pubchem_id`/`chebi_id`/
   `refmet_id`/`lipidmaps_id`/`lion_id`), for external tools such as
   `pathview`.

---

## Install

```r
# At the project root (renv auto-activates via .Rprofile):
renv::restore()      # one-time, populates the project library
devtools::install()
```

R >= 4.3, Bioconductor >= 3.18 are assumed (see CLAUDE.md §9).

---

## Quick start (using the bundled synthetic fixtures)

```r
library(MetaboSetR)

fx <- function(name) {
  system.file("extdata", "synthetic", name, package = "MetaboSetR")
}

# 1. Import
conc <- read_webidq(fx("sample_conc.xlsx"), fx("sample_status.xlsx"))
qc   <- read_webidq_qc(fx("qc_conc.xlsx"),  fx("qc_status.xlsx"))
thr  <- read_thresholds(c(fx("thresholds_KIT01.xlsx"),
                          fx("thresholds_KIT02.xlsx")))

# 2. Preprocess. Synthetic data is intentionally noisy, so loosen CV here;
# real-data defaults are cv_threshold = 20, lod_rate_threshold = 0.5.
prep <- preprocess(conc, qc, thresholds = thr,
                   cv_threshold          = 200,
                   apply_lod_rate_filter = FALSE)

# 3. QC diagnostics (optional; never auto-invoked)
qc_pca_plot(prep, mode = "with_samples")
qc_cv_plot(prep)
qc_filter_summary(prep)
qc_permanova(prep)

# 4. Metabolite-level testing
res_feat <- test_metabolites(prep, group = "kit_id", method = "limma")

# 5. Pathway-set enrichment
res_gsea <- gsea_pathway_sets(prep, group = "kit_id",
                              domain = "reactome",
                              minSize = 1L, maxSize = 2000L)

# ORA over the same sets; universe = the measured metabolite panel.
res_ora <- ora_pathway_sets(
  sig      = res_feat$metabolite[res_feat$P.Value < 0.05],
  universe = res_feat$metabolite,
  domain   = "reactome")
```

---

## Real-data workflow

The only thing that changes is which xlsx paths you pass; defaults are
tuned for actual exports.

```r
# Real WebIDQ exports — do not commit (see CLAUDE.md §5.1).
conc <- read_webidq("data/sample_conc.xlsx",   "data/sample_status.xlsx")
qc   <- read_webidq_qc("data/qc_conc.xlsx",    "data/qc_status.xlsx")
thr  <- read_thresholds(c("data/lod_kit01.xlsx",
                          "data/lod_kit02.xlsx"))

prep <- preprocess(conc, qc, thresholds = thr)
res  <- test_metabolites(prep, group = "treatment")
res_gsea <- gsea_pathway_sets(prep, group = "treatment",
                              domain = "immunomet")
```

---

## Pathway-set inspection

```r
list_pathway_sets()
#>      domain                                    set_name n_members
#> 1 immunomet               IMMUNOMET_ONEIL_GLYCOLYSIS         3
#> ...
#>  reactome  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE                 ...

get_pathway_sets("reactome")[["REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"]]
get_pathway_sets_meta("immunomet")       # long-format metadata
pathway_set_gmt_path("reactome")         # for external tools

# External-DB metabolite-level domains (small molecules via RaMP-DB,
# lipids via LION).
get_pathway_sets("lion")      # lipid functional terms
get_pathway_sets("source")    # SOURCE_PLANT / SOURCE_MICROBE
get_pathway_sets("health")    # disease-association sets (hypothesis-level)
```

### Overlaying small molecules on a KEGG map (pathview)

`annotate_metabolites()` maps Biocrates `short_name`s to their RaMP-DB
cross-reference ids — `kegg_id` for ~290 small molecules (plus
`pubchem_id`/`chebi_id`/`refmet_id`/`lipidmaps_id`/`lion_id`). The package
does not wrap `pathview`; call it yourself:

```r
ann <- annotate_metabolites(colnames(prep@sample@assay))  # short_name -> ids
fc  <- ...  # named numeric: names = kegg_id, values = log2FC
pathview::pathview(cpd.data = fc, pathway.id = "hsa00020", species = "hsa")
```

---

## Project documentation

- `docs/DESIGN.md` — design single source of truth ("what & why")
- `docs/ROADMAP.md` — version roadmap, current status, open questions
- `docs/decisions.md` — append-only decision log
- `CLAUDE.md` — work rules for the Claude Code agent

The vendor reference TSV and synthetic fixtures are under `data-raw/`
(curation source) and `inst/extdata/` (installed assets). The
`data-raw/pathway_sets/` subdirectories carry their own README per domain.

---

## Tests

```r
devtools::test()
# [ FAIL 0 ]
```

`devtools::install()` is the install-clean bar; `R CMD check` is not
required.

---

## License

Private and proprietary for now; see `LICENSE` for the full text. Intended
for public release alongside the accompanying paper. Do not redistribute
until then.

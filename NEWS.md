# MetaboSetR 0.1.0

Initial release. Derived from the internal MetaboIndicatoR package, with the
Biocrates indicator-computation layer removed; MetaboSetR keeps the WebIDQ
preprocessing pipeline and operates entirely at the metabolite level.

## Features

- Biocrates MxP Quant 1000 WebIDQ reader (concentration + status text
  primary, cell-color fallback) and threshold (LOD/LLOQ/ULOQ) reader.
- Sample / metabolite QC: CV per kit, LOD-rate, ULOQ-rate, IQR-based
  sample outlier flagging, plus QC diagnostic plots and PERMANOVA.
- Metabolite-level feature testing via `test_metabolites()`: limma
  moderated t and Wilcoxon, BH/BY adjustment.
- Curated pathway-set enrichment over metabolite-level sets:
  - `gsea_pathway_sets()` — rank-based GSEA via `fgsea::fgseaMultilevel`.
  - `ora_pathway_sets()` — threshold-based hypergeometric over-representation
    with the measured panel as background.
  - Domains: `immunomet` (O'Neill 2016), and external-DB metabolite sets
    via RaMP-DB — `reactome`, `wikipathways`, `smpdb`, `lion`, `source`,
    `health`.
- `annotate_metabolites()` — cross-reference id annotation
  (kegg/hmdb/pubchem/chebi/refmet/lipidmaps/lion) for KEGG-map overlays.
- `gsea_to_df()` and `log2_with_pseudocount()` helpers; Biocrates status
  predicates (`is_valid()`, `is_below_lod()`, …).

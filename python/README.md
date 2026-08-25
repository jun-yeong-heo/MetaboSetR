# MetaboSetR — Python port

Pure-Python port of the R package. Design is the shared
[`../docs/DESIGN.md`](../docs/DESIGN.md); migration notes and decisions are in
[`../docs/MIGRATION.md`](../docs/MIGRATION.md).

## Scope

The full pipeline is ported:

- **Preprocessing** — WebIDQ readers (text + color-fill status), optional
  threshold reader, sample-outlier detection, two-step metabolite filter
  (CV → LOD-rate), log2 transform, four data containers.
- **Approach A** — `test_metabolites` (InMoose moderated-t default, SciPy
  Wilcoxon option).
- **Approach B** — `gsea_pathway_sets` (gseapy prerank; seeded, deterministic).
- **Approach C** — `ora_pathway_sets` (SciPy hypergeometric).
- **Pathway loaders / annotation / QC diagnostics**.

## Modules
`classes` · `status` · `transform` · `read_webidq` · `read_status_color` ·
`read_thresholds` · `filter_samples` · `filter_metabolites` · `preprocess` ·
`stats_feature` · `stats_gsea` · `stats_ora` · `pathways` · `annotate` ·
`qc_diagnostics` (+ `_xlsx`, `reference`, `_padjust`).

## Dependencies
- Preprocessing only: `numpy`, `pandas`, `openpyxl`.
- Full pipeline adds: `scipy`, `inmoose` (limma), `gseapy` (GSEA),
  `scikit-bio` (PERMANOVA), `matplotlib` (plots). No `statsmodels`
  (`p.adjust` BH/BY is reimplemented in `_padjust.py`).

## Usage
```python
import numpy as np, metabosetr as m

conc = m.read_webidq("sample_conc.xlsx", "sample_status.xlsx")
qc   = m.read_webidq_qc("qc_conc.xlsx", "qc_status.xlsx")
th   = m.read_thresholds(["thresholds_KIT01.xlsx", "thresholds_KIT02.xlsx"])  # optional
prep = m.preprocess(conc, qc, thresholds=th)

group = np.array(["A", "A", "A", "B", "B", "B"])
feat  = m.test_metabolites(prep, group)                       # A: moderated-t
gsea  = m.gsea_pathway_sets(prep, group, "immunomet", seed=1) # B
ora   = m.ora_pathway_sets(feat.metabolite[feat["P.Value"] < 0.05].tolist(),
                           feat.metabolite.tolist(), "immunomet")  # C
```

## Tests
```bash
cd python && python3 -m pytest -q
```
Verified against R (base `stats`, on identical inputs, to ~1e-14):
CV / LOD-rate / IQR filters, Wilcoxon p-values & log2FC, `p.adjust` BH/BY
(including NA handling), hypergeometric ORA. Text-vs-color status readers are
cross-validated; pathway-set / metabolite-dict counts match
`make_sysdata.R`'s assertions. The limma path relies on InMoose's own
drop-in fidelity to R limma; GSEA differs from fgsea by design but is
deterministic under `seed`.

## Notes
- `reference.py` / `pathways.py` read the canonical tables and pathway
  masters from `../data-raw/`. Bundling into an installed package is deferred
  (MIGRATION.md §8 Q5).
- `_xlsx.py` tolerates WebIDQ exports that declare drawing relationships
  whose parts are absent (which crash a plain `openpyxl.load_workbook`).
- InMoose labels design columns `column0`/`column1`; the group effect is
  `column1`, and `topTable(sort_by="none")` is unsupported (results are
  reindexed to feature order instead).

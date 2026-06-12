# Immunometabolism Pathway Sets

Immunometabolism GSEA-style sets for the Biocrates Quant 1000 panel,
constructed top-down from the O'Neill 2016 framework. Sets are
**metabolite-level**.

이 디렉토리는 `data-raw/pathway_sets/immunomet/` 에 위치하며, 패키지 build 시점에 여기의 master TSV가 `R/sysdata.rda`의 `pathway_sets$immunomet`로 편입되고, GMT 파일이 `inst/extdata/pathway_sets/`에 동기화된다. 상위 공통 규칙은 `../README.md`를, 패키지 전체 설계는 `docs/DESIGN.md` §8.5.2를 참조.

---

## Files

| File | Type | Description |
|---|---|---|
| `immunomet_metabolites_master.tsv` | source of truth | Long-format metabolite-to-set master |
| `immunomet_metabolites.gmt` | build output (auto) | GMT export, metabolite level |
| `assignment_log_metabolites.tsv` | build output (auto) | Per-metabolite audit log (all 1234 metabolites) |
| `pathways.py` | source (edit this) | Pathway definitions: canonical metabolites, classes, exclusions |
| `pathway_definitions.md` | documentation | Biological rationale, references, set definitions |
| `build.py` | build script | Reproducible regeneration from sources |

---

## Curation method: top-down, reference-based

### Primary reference
O'Neill LAJ, Kishton RJ, Rathmell J (2016). A guide to immunometabolism for immunologists. *Nat Rev Immunol* 16(9):553–565. doi:10.1038/nri.2016.70

Supporting references (for individual pathway scoping):
- Pearce EL, Pearce EJ (2013). Metabolic pathways in immune cell activation and quiescence. *Immunity* 38(4):633–643.
- Buck MD, O'Sullivan D, Pearce EL (2015). T cell metabolism drives immunity. *J Exp Med* 212(9):1345–1360.
- Mogilenko DA, Sergushichev A, Artyomov MN (2023). Systems immunology approaches to metabolism. *Annu Rev Immunol* 41:317–342.

### Scope
Five of the six core O'Neill pathways are implemented. The pentose phosphate pathway is excluded because the Biocrates Quant 1000 panel does not cover its canonical metabolites (G6P, 6-phosphogluconate, ribose-5-phosphate, NADPH). Implemented sets:

- `IMMUNOMET_ONEIL_GLYCOLYSIS`
- `IMMUNOMET_ONEIL_TCA`
- `IMMUNOMET_ONEIL_FAO`
- `IMMUNOMET_ONEIL_FAS`
- `IMMUNOMET_ONEIL_AA_METABOLISM`

The `_ONEIL_` infix in the set names is deliberate: it pins membership to this specific framework so that future immunometabolism frameworks (e.g., hypothetical `IMMUNOMET_MOGILENKO_*`) can coexist without namespace collision.

### Metabolite assignment
Each of the 1234 panel metabolites is assigned to zero or more O'Neill pathways based on whether it is a canonical pathway member (by exact name) or belongs to a canonical pathway class, using the source-of-truth pathway definitions in `pathways.py`. Result: 654 unique metabolites assigned.

---

## Set sizes

| Set | Metabolites |
|---|---:|
| `IMMUNOMET_ONEIL_GLYCOLYSIS` | 3 |
| `IMMUNOMET_ONEIL_TCA` | 10 |
| `IMMUNOMET_ONEIL_FAO` | 58 |
| `IMMUNOMET_ONEIL_FAS` | 445 |
| `IMMUNOMET_ONEIL_AA_METABOLISM` | 138 |

The FAS set is intentionally broad — it includes all glycerolipid and choline-containing phospholipid classes because plasma metabolomics cannot cleanly separate de novo lipogenesis from dietary/lipoprotein-derived pools.

---

## Overlap policy

Overlap between sets is permitted when biologically meaningful (see package-level policy in `../README.md`). Artifactual overlap from chemical-class grouping is explicitly removed.

In the final metabolite-level output there are **0 overlaps**. TCA-cycle intermediates (succinate, fumarate, malate, oxaloacetate, α-ketoglutarate, 2-hydroxyglutarate, itaconate) are chemically classified under Dicarboxylic acids in the Biocrates panel, which would otherwise place them in both TCA and FAO. They are explicitly excluded from FAO (exclusion list in `pathways.py`) to preserve pathway distinction.

---

## Master TSV schema

**Metabolite level** (`immunomet_metabolites_master.tsv`):

| Column | Description |
|---|---|
| `shortname` | Metabolite short name |
| `fullname` | Metabolite full name |
| `analyte_class` | Metabolite analyte class |
| `set_name` | One of five `IMMUNOMET_ONEIL_*` identifiers |
| `match_basis` | `individual` (canonical name match) or `class:<class>` (class match) |

---

## Direction values

The immunomet master TSV has **no direction column** (loaded as NA). Immunometabolism does not have a single reference state: the key contrasts (effector vs memory T cells, M1 vs M2 macrophages, Teff vs Treg) are bidirectional, and pathway activation direction depends on which contrast the user is running. Downstream tooling should interpret directionality in the context of the specific comparison.

---

## Build procedure

From this directory:

```bash
python3 build.py
```

Regenerates the metabolite master TSV, GMT file, and assignment log from `../../reference/biocrates_Quant1000_metabolites.tsv` and the pathway definitions in `pathways.py`.

To propagate changes into the installed package, follow the top-level procedure in `../README.md` (run `data-raw/make_sysdata.R`).

---

## Methods statement for manuscripts

Immunometabolism pathway sets were constructed top-down from the O'Neill et al. (2016) *Nat Rev Immunol* framework. For each of the five implemented pathways (glycolysis, TCA cycle, fatty acid oxidation, fatty acid synthesis, amino acid metabolism; pentose phosphate pathway omitted due to panel coverage), canonical metabolite membership was defined from the primary reference and supporting reviews (Pearce and Pearce 2013; Buck et al. 2015; Mogilenko et al. 2023), then mapped to the Biocrates Quant 1000 panel at the metabolite level. Where the Biocrates chemical-class grouping placed TCA-cycle intermediates within the Dicarboxylic acids chemical class, these metabolites were explicitly excluded from the FAO set to preserve biological pathway distinction while retaining the Dicarboxylic acids class membership for true ω-oxidation products. All membership rules, exclusions, and audit logs are reproducible from the build script. Overlap between sets was permitted when biologically meaningful.

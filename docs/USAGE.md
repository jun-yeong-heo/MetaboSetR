# MetaboSetR — 사용 레퍼런스

> 이 문서는 LLM과 사용자가 분석 스크립트 작성 시 빠르게 참조하기 위한
> 사용 레퍼런스입니다. 설계 근거는 `docs/DESIGN.md`, 변경 이력은
> `docs/decisions.md`를 보세요.

## 1. 패키지 한 줄 요약

Biocrates MxP Quant 1000 (WebIDQ) 데이터를 받아서 (1) 전처리,
(2) metabolite-level feature 통계, (3) 큐레이션된 pathway-set GSEA +
pathway-set ORA를 수행하는 R 패키지.

## 2. 도메인 모델 (LLM이 알아야 할 핵심)

### 2.1 데이터 단위

- **Metabolite**: Biocrates Q1000 panel의 단일 분자. 1234종. 모든 통계·pathway set의 단위.
- **Sample**: 측정 row 한 개. 6개의 메타 컬럼(`Sample identification`,
  `Sample description`, `Submission name`, `Collection date`, `Species`,
  `Material`)과 1231개 metabolite 컬럼으로 구성.
- **QC sample**: pooled QC. 보통 KIT01·KIT02 각각 3 replicate, 합 6 samples.
- **Kit ID**: `Submission name` 컬럼에서 정규식 `KIT\d+`로 추출.

### 2.2 Status 카테고리 (5종 + Missing)

```
"Valid"        — 정상 정량
"< threshold"  — Valid로 취급
"< LLOQ"       — Valid로 취급, 단 원본 annotation은 보존
"< LOD"        — 검출한계 미만
"> ULOQ"       — 상한 초과
NA / Missing   — 방어적 처리용
```

`is_valid()` 헬퍼는 `Valid`, `< threshold`, `< LLOQ` 셋 다 TRUE 반환.

### 2.3 S4 class 4종

| 클래스 | slot | 설명 |
|---|---|---|
| `ConcentrationData` | `assay`(num matrix), `status`(chr matrix), `sample_meta`(df), `metabolite_names`(chr) | sample × metabolite |
| `QCData` | 동일 구조 | pooled QC × metabolite |
| `ThresholdData` | `lod`/`lloq`/`uloq`(num matrix), `kit_ids`, `metabolite_names` | kit × metabolite, NA 허용 |
| `PreprocessedData` | `sample`, `qc`, `thresholds`, `sample_outliers`, `metabolite_pass`, `sample_filter_log`, `metabolite_filter_log`, `params` | preprocess() 출력 — 통계 함수의 입력 |

### 2.4 핵심 설계 원칙 (LLM이 위반하지 말아야 할 것)

- **Tidyverse 사용 금지** (`dplyr`, `tidyr`, `purrr`, `magrittr`, `readr` 등). base R + Bioconductor만. (CLAUDE.md §2.1)
- **Raw concentration 그대로 사용**. filter 결과는 flag로만 attach, raw 값은 절대 mutate 안 함.
- **Log2 변환은 통계 layer 진입 시점에만**. assay 자체는 raw. (§4.5)
- **이름 매칭은 exact match only**, fuzzy 금지. (§2.7)
- **실제 sample 데이터 commit 금지** (CLAUDE.md §5.1)
- **Vendor reference TSV 수정 금지** (CLAUDE.md §5.2)

## 3. 표준 파이프라인

```r
library(MetaboSetR)

# 1. Import
conc <- read_webidq("sample_conc.xlsx", "sample_status.xlsx")
qc   <- read_webidq_qc("qc_conc.xlsx", "qc_status.xlsx")
thr  <- read_thresholds(c("lod_kit01.xlsx", "lod_kit02.xlsx"))   # optional

# 2. Preprocess
prep <- preprocess(conc, qc, thresholds = thr)

# 3. QC diagnostics (optional, 자동 호출 안 됨)
qc_pca_plot(prep, mode = "with_samples")
qc_filter_summary(prep)
qc_permanova(prep)

# 4. Metabolite-level 통계
res_feat <- test_metabolites(prep, group = "treatment", method = "limma")

# 5. Pathway-set enrichment
res_gsea <- gsea_pathway_sets(prep, group = "treatment", domain = "reactome")
res_ora  <- ora_pathway_sets(
  sig      = res_feat$metabolite[res_feat$P.Value < 0.05],
  universe = res_feat$metabolite,
  domain   = "reactome")

# 6. Pathway-set inspection
list_pathway_sets()
get_pathway_sets("immunomet")
```

## 4. 함수 레퍼런스

### 4.1 Reader

#### `read_webidq(concentration_file, status_file = NULL, status_method = c("text", "color"), synonyms = NULL)`
- **반환**: `ConcentrationData`
- `status_method = "text"` (default)이면 `status_file` 필수.
- `status_method = "color"`면 `concentration_file`의 cell fill color에서 status를 추출 (DESIGN.md §2.4).
- `synonyms`: 사용자 지정 metabolite 이름 → canonical short_name 매핑 (named character vector).

#### `read_webidq_qc(concentration_file, status_file = NULL, status_method = c("text", "color"), synonyms = NULL)`
- **반환**: `QCData`. 같은 시그니처.

#### `read_thresholds(threshold_files, synonyms = NULL)`
- **반환**: `ThresholdData`. `threshold_files`는 per-kit xlsx 경로 character vector.

### 4.2 Preprocess

#### `preprocess(conc, qc, thresholds = NULL, sample_outlier_method = c("IQR", "percentile", "both"), outlier_metric = c("below_lod", "above_uloq", "either"), percentile_cutoff = 0.95, cv_threshold = 20, cv_mode = c("per_kit", "pooled"), lod_rate_threshold = 0.5, apply_lod_rate_filter = TRUE, winsorize_uloq = FALSE)`
- **반환**: `PreprocessedData`
- 흐름: sample-level outlier → CV filter → LOD-rate filter (CV pass인 metabolite만) → optional ULOQ winsorization
- `cv_mode = "per_kit"` (default): 각 kit별 CV가 모두 ≤ threshold이어야 pass (AND)
- `cv_mode = "pooled"`: 모든 QC를 한 군으로 CV 계산
- `winsorize_uloq = TRUE`는 `thresholds` 필수 (NULL이면 warning + 비활성)

### 4.3 QC diagnostics (사용자가 명시 호출, 자동 안 됨)

| 함수 | 반환 | 메모 |
|---|---|---|
| `qc_pca_plot(prep, mode = c("qc_only", "with_samples"))` | ggplot | log2-transformed assay에서 PCA |
| `qc_permanova(prep)` | adonis2 anova table | `kit_id` 효과 정량 |
| `qc_cv_plot(prep)` | ggplot | per-metabolite CV histogram + threshold line |
| `qc_lod_rate_plot(prep, by = c("sample", "metabolite"))` | ggplot | `<LOD` rate 분포 |
| `qc_filter_summary(prep)` | data.frame | sample/metabolite total/pass/fail |

### 4.4 Statistical layer

#### Approach A — Metabolite feature-level
```r
test_metabolites(prep, group, method = c("limma", "wilcoxon"),
                 padjust = c("BH", "BY"))
```
- `prep`: `PreprocessedData`. Outlier sample은 자동 제외.
- `group`: `prep@sample@sample_meta`의 컬럼명 또는 `nrow(prep@sample@assay)` 길이의 vector
- `method = "wilcoxon"`은 정확히 2 그룹 필수
- `padjust`: 추가할 보정 방법. 컬럼명은 `padj_<method>` 형태
- **반환**: data.frame (metabolite 1행), columns:
  - limma: `metabolite`, `log2FC`, `AveExpr`, `t`, `P.Value`, `padj_*`
  - wilcoxon: `metabolite`, `log2FC`, `P.Value`, `padj_*`

#### Approach B — Curated pathway-set GSEA
```r
gsea_pathway_sets(prep, group, domain, ...)   # ... → fgsea::fgseaMultilevel
```
- `prep`: `PreprocessedData`. metabolite log2FC(group2 vs group1)가 rank stat (outlier 제외).
- `domain`: in-house `"immunomet"` + 외부-DB metabolite-level `{"reactome", "wikipathways", "smpdb", "lion", "source", "health"}`. `source` = RaMP origin (SOURCE_PLANT/SOURCE_MICROBE, 2 set); `health` = RaMP 질병-연관 (disease-MSEA, 51 set — 문헌편향·순환논리 caveat, 가설 수준 해석).
- 모든 도메인 metabolite-level (멤버 = metabolite short_name).
- **반환**: fgsea data.table (`pathway`, `pval`, `padj`, `ES`, `NES`, `size`, `leadingEdge`)
- **universe 분리**: reactome/wikipathways/smpdb는 소분자, lion은 lipid만 다룸. panel에 멤버가 전혀 안 겹치면 "no overlap" 에러. 소분자/lipid를 한 분석에 섞을 때 해석 주의.
- **KEGG-map 시각화(pathview)**: 패키지 비의존. `annotate_metabolites()`(소분자 cross-ref)로 얻은 `kegg_id`로 사용자가 `pathview::pathview(cpd.data=..., pathway.id="hsa00020")` 직접 호출.

#### Approach C — Curated pathway-set ORA (threshold-based)
```r
ora_pathway_sets(sig, universe, domain, padjust = c("BH", "BY"))
```
- GSEA(rank 기반)의 threshold 기반 짝. **GSEA와 동일한 내장 `pathway_sets`**에서 동작 (DESIGN.md §6.6).
- `sig`: 유의미 metabolite 이름 (threshold로 선별). `universe`: 측정한 전체 metabolite 이름 (background panel, 인자 필수).
- **Universe = 측정 panel**: 각 pathway 멤버를 `universe`와 교집합한 뒤 검정하므로, 측정 안 한 멤버는 분모에 안 들어감.
- 검정: pathway별 단측 hypergeometric (`stats::phyper`). `sig`가 `universe`에 없으면 경고 후 drop.
- Direction은 GSEA와 동일하게 검정에 미전달 (해석용; `get_pathway_sets_meta()`로 join).
- **반환**: data.frame (pathway별 1행, `pval` 오름차순) — `set_name`, `set_size`, `overlap`, `expected`, `pval`, `padj_*`, `overlap_members`(`;` 결합). list-column이 없어 `write.table`로 바로 저장 가능.

#### TSV/CSV 저장
```r
gsea_to_df(res, sep = ";")    # list-column 을 collapsed string 으로 변환한 data.frame
```
- `gsea_pathway_sets()` 결과의 `leadingEdge` 는 list-column 이라 `write.table`/`write.csv` 가 인코딩 못 함.
- `gsea_to_df(res)` 로 한 번 변환 후 저장. 원본 `res` 는 안 건드리므로 R 안에서는 `res$leadingEdge[[i]]` 그대로 사용 가능.

### 4.5 Pathway-set loader

```r
list_pathway_sets()                # df: domain/set_name/n_members
get_pathway_sets(domain)           # named list of character vectors
get_pathway_sets_meta(domain)      # long-format metadata df
pathway_set_gmt_path(domain)       # GMT file path
```

### 4.6 log2 transform helper

```r
log2_with_pseudocount(x, pseudocount = NULL)   # rows = features, cols = samples
```
- Per-feature half-min positive pseudocount when `pseudocount = NULL` (default).
- Pass `t(x)` if your input has features as columns (e.g. sample × metabolite).
- Same offset rule the package uses internally for GSEA / DE / QC PCA.

### 4.7 Status helper

```r
is_valid(status)        # status %in% c("Valid", "< threshold", "< LLOQ")
is_below_lloq(status)   # status == "< LLOQ" only
is_below_lod(status)    # status == "< LOD"
is_above_uloq(status)   # status == "> ULOQ"
is_missing(status)      # NA or "Missing"
```

### 4.8 Metabolite id annotation

```r
annotate_metabolites(x)   # x = character vector of short_names
```
- short_name 벡터를 내장 `metabolite_dict`의 cross-ref id에 매핑하는 공개 accessor.
- **반환**: data.frame, 입력 1개당 1행 (입력 순서 보존, 미매칭은 NA). 컬럼: `metabolite`(입력 이름), `long_name`, `class`, `kegg_id`, `hmdb_id`, `pubchem_id`, `chebi_id`, `refmet_id`, `lipidmaps_id`, `lion_id`.
- pathview 등 외부 툴용 id 추출에 사용.

## 5. Internal data (sysdata.rda — `:::`로 접근)

```r
MetaboSetR:::metabolite_dict      # 1234 × 10: short_name, long_name, class, kegg_id, hmdb_id, pubchem_id, chebi_id, refmet_id, lipidmaps_id, lion_id
MetaboSetR:::status_color_map     # 5 × 2: hex, status_text
MetaboSetR:::pathway_sets         # nested list: domain → set_name → members
MetaboSetR:::pathway_sets_meta    # 12813 × 5: domain, set_name, member, direction, brief_note
```

`metabolite_dict` id annotation: `kegg_id`/`hmdb_id`/`pubchem_id`/`chebi_id`/`refmet_id`/`lipidmaps_id` 는 RaMP-DB cross-ref (커버리지 289/747/746/663/519/475), `lion_id` 는 매칭된 LION 정규 lipid 이름 (807). 모두 character, 미매칭은 `NA`. **공개 accessor는 `annotate_metabolites()`** (§4.8) — 분석 스크립트는 `:::` 대신 이것을 쓸 것.

`pathway_sets` 구조 (모두 metabolite-level):
```
$immunomet$IMMUNOMET_ONEIL_GLYCOLYSIS = c("Glucose", "Lac", "Pyruvic acid")
$reactome$REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE = c("Citric acid", "Fumaric acid", ...)
$smpdb$SMPDB_GLUTAMINOLYSIS_AND_CANCER = c("Gln", "Glu", ...)
$lion$LION_ABOVE_AVERAGE_BILAYER_THICKNESS = c("PA 16:0_18:1", ...)
$source$SOURCE_MICROBE = c("TMAO", ...)
$health$HEALTH_CANCER = c("Ala", "Arg", "Asn", ...)
```

immunomet metabolite master에는 direction 컬럼이 없어 `pathway_sets_meta`의 immunomet direction은 NA. 외부-DB 도메인은 direction `-`.

## 6. 자주 쓰는 패턴 cookbook

### 6.1 그룹 컬럼 attach 후 분석

```r
conc <- read_webidq("sample_conc.xlsx", "sample_status.xlsx")
conc@sample_meta$treatment <- factor(c("ctrl","ctrl","ctrl","drug","drug","drug"))

prep <- preprocess(conc, qc)
# group은 prep@sample@sample_meta의 컬럼명으로 참조
res <- test_metabolites(prep, group = "treatment")
```

### 6.2 Wilcoxon은 2그룹만 — multi-group은 limma

```r
res <- test_metabolites(prep, group = "treatment", method = "limma")
# wilcoxon은 3+ group에서 에러: "Wilcoxon test requires exactly two group levels"
```

### 6.3 유의미 metabolite로 ORA

```r
res  <- test_metabolites(prep, group = "treatment")
sig  <- res$metabolite[res$P.Value < 0.05]
ora  <- ora_pathway_sets(sig, universe = res$metabolite, domain = "reactome")
```

### 6.4 GSEA 결과 direction cross-reference

```r
res  <- gsea_pathway_sets(prep, group = "treatment", domain = "immunomet")
meta <- get_pathway_sets_meta("immunomet")
# leadingEdge에 들어간 metabolite의 direction(있으면) 조회
res$leadingEdge_directions <- lapply(res$leadingEdge, function(le) {
  meta$direction[meta$member %in% le]
})
```

### 6.5 사용자 metabolite 이름이 캐논과 다를 때

```r
# 예: 사용자 데이터엔 "Aspartate"라고 적혀 있음 (canonical = "Asp")
conc <- read_webidq("sample.xlsx", "status.xlsx", synonyms = c("Aspartate" = "Asp"))
```

reader가 자동으로 풀 매핑 (`Aspartic acid` → `Asp`)도 시도하므로 대부분의
WebIDQ export는 synonyms 인자 없이도 동작.

### 6.6 ULOQ 윈저화

```r
prep <- preprocess(conc, qc, thresholds = thr, winsorize_uloq = TRUE)
# 결과: prep@sample@assay에서 status가 "> ULOQ"인 cell이 kit별 ULOQ로 대체
```

### 6.7 Pathway set GMT 외부 도구로 export

```r
gmt_path <- pathway_set_gmt_path("reactome")
# 그 후 GSEA-MSigDB, gprofiler 등 외부 툴에 그대로 input
```

## 7. 자주 나오는 에러와 대처

| 에러 메시지 | 원인 | 대처 |
|---|---|---|
| `status_file is required when status_method = 'text'` | text mode인데 status xlsx 없음 | `status_method = "color"` 또는 status_file 제공 |
| `Could not resolve these metabolite names to canonical short names` | 사용자 컬럼명이 canonical/fullname 둘 다 아님 | `synonyms = c(typed = canonical)` 인자 추가 |
| `All samples flagged as outliers` | sample_outlier_method 너무 strict 또는 데이터 자체가 noisy | `sample_outlier_method = "percentile"` 또는 `apply_lod_rate_filter = FALSE` |
| `Wilcoxon test requires exactly two group levels` | 3+ 그룹에 wilcoxon 사용 | `method = "limma"` 사용 |
| `GSEA requires exactly two group levels` | 같은 이유 | 그룹을 binary로 재구성 |
| `No pathway has any member overlapping the rank vector` | panel이 해당 도메인 metabolite를 안 가짐 (예: 소분자 패널에 lion) | 도메인-패널 매칭 확인 |
| `Unknown domain` | 존재하지 않는 도메인 | `list_pathway_sets()`로 도메인 확인 |

## 8. 흔한 설계 의문에 대한 답

- **왜 sample × metabolite layout인가?**
  WebIDQ export가 sample × metabolite 형태라 reader가 native하게 그 형태로 만든다. 통계 함수에서 필요 시 transpose한다.

- **왜 `<LLOQ`를 Valid로 취급?**
  LLOQ 미만이라도 quantification 가능. 후속 분석에서 어차피 valid로 사용. 단 원본 annotation은 보존.

- **`<LOD` 값이 0인데 통계에 그대로?**
  Raw 값 보존 원칙. log2 변환 시 adaptive pseudocount (= non-zero min / 2)로 처리.

- **N=6 같이 적은 sample로 PCA/PERMANOVA 해도?**
  동작은 하지만 통계적 신뢰도는 낮다. n_samples=6에 PERMANOVA는 모든 permutation을 enumerate하므로 informational 메시지 출력 (정상).

- **indicator 계산은 어디?**
  MetaboSetR에는 없다. Biocrates indicator 레이어는 license 사유로 제외됨 (decisions.md #1). 모든 분석은 metabolite-level.

## 9. 패키지 구조 (LLM이 코드 위치 알아야 할 때)

```
R/
  ConcentrationData.R   QCData.R   ThresholdData.R   PreprocessedData.R   (S4 classes)
  read-webidq.R              (read_webidq, read_webidq_qc)
  read-status-color.R        (color fallback, validate_status stub)
  read-thresholds.R          (read_thresholds)
  preprocess.R               (preprocess + .winsorize_uloq)
  filter-samples.R           (.detect_sample_outliers)
  filter-metabolites.R       (.cv_filter, .lod_rate_filter)
  qc-diagnostics.R           (qc_pca_plot, qc_permanova, qc_cv_plot, qc_lod_rate_plot, qc_filter_summary)
  stats-feature.R            (test_metabolites, .test_limma, .test_wilcoxon)
  stats-gsea.R               (gsea_pathway_sets, gsea_to_df, .metabolite_log2fc)
  stats-ora.R                (ora_pathway_sets)
  pathway-sets.R             (list/get_pathway_sets, get_pathway_sets_meta, pathway_set_gmt_path)
  annotate.R                 (annotate_metabolites)
  status-helpers.R           (is_valid, is_below_lloq, is_below_lod, is_above_uloq, is_missing)
  transform.R                (log2_with_pseudocount)
  utils.R                    (.resolve_group, package-level doc)
  data.R                     (sysdata documentation)
  sysdata.rda                (binary, generated by data-raw/make_sysdata.R)

inst/extdata/
  pathway_sets/*.gmt         (7 files: immunomet/reactome/wikipathways/smpdb/lion/source/health)
  synthetic/*.xlsx           (fixtures)

data-raw/
  reference/                 (vendor TSV: metabolites)
  pathway_sets/              (curation source: immunomet, ramp, lion)
  make_sysdata.R             (sysdata.rda 빌더)
  make_synthetic.R           (synthetic xlsx 빌더)

tests/testthat/

docs/
  DESIGN.md   ROADMAP.md   decisions.md   USAGE.md
```

## 10. 미구현 stub

- `validate_status(prep)` — text vs color cross-validation. PreprocessedData에 두 source 동시 retain하는 reader path 필요. 보류.

## 11. 빠른 sanity check

```r
library(MetaboSetR)
nrow(MetaboSetR:::metabolite_dict)        # 1234
nrow(list_pathway_sets())                 # 802 (5 + 224 + 145 + 266 + 109 + 2 + 51)
length(get_pathway_sets("immunomet"))     # 5
length(get_pathway_sets("reactome"))      # 224
length(get_pathway_sets("lion"))          # 109
```

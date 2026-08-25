# MetaboSetR — Design Document

> **Version**: v3.0 (MetaboSetR)
> **Purpose**: MetaboSetR R package의 설계 single source of truth.
> **Scope**: "무엇을 왜 그렇게 설계했는가" (what & how). 구현 순서와 진행 상태는 `ROADMAP.md`, 의사결정 이력은 `decisions.md` 참조.
> MetaboSetR은 내부 패키지 MetaboIndicatoR에서 indicator 레이어를 제거하여 파생되었다 (→ decisions.md #1).
> **Python 포트**: 전 파이프라인이 순수 Python으로 재이식되어 `python/metabosetr/`에 있다 (→ decisions.md #6). 본 문서는 설계 SSOT로 언어 중립이며, R 판본과 Python 포트가 동일 설계를 따른다. Python 실현 상세(파일-단위 매핑·의존성·수치 재현성)는 `MIGRATION.md` 참조.

---

## 1. 프로젝트 개요

### 1.1 목적
Biocrates MxP® Quant 1000 kit의 WebIDQ export 데이터를 input으로 받아서:
1. **전처리** (QC, filtering)
2. **Metabolite-level 통계** — 각 metabolite를 feature로 group comparison
3. **큐레이션 pathway-set enrichment** — metabolite-level set에 대한 두 가지 접근:
   - GSEA (rank-based, `fgsea`)
   - ORA (threshold-based hypergeometric)

를 수행하는 R package.

### 1.2 기본 정보
- **패키지 이름**: `MetaboSetR`
- **라이선스 / 배포**: **GPL-3.0 (copyleft)** (→ decisions.md #9). GitHub private repo로 시작, 논문 공개 시점에 함께 공개 (→ decisions.md #1, #5). 코드는 GPL-3.0, 번들 데이터 라이선스는 별도 명시 예정, 외부-DB 세트 attribution은 §8.5대로 유지
- **데이터 단위**: 모든 통계·pathway set은 metabolite-level. (indicator 레이어는 license 사유로 제외 — → decisions.md #1)
- **코딩 관례**: R에서 `summarise()` 등 tidyverse 사용 금지, base R + Bioconductor만 사용
- **Tidyverse dependency 의도적 제외**

---

## 2. Input 데이터 규격

### 2.1 입력 파일 구성
사용자가 제공해야 하는 파일:

| # | 파일 | 필수/선택 | 설명 |
|---|---|---|---|
| 1 | Sample concentration file | 필수 | WebIDQ export, 셀이 numeric concentration |
| 2 | Sample status file | 필수 (또는 fallback) | 동일 구조, 셀이 status text |
| 3 | QC concentration file | 필수 | Pooled QC 샘플의 concentration |
| 4 | QC status file | 필수 (또는 fallback) | QC의 status text |
| 5 | LOD/LLOQ/ULOQ threshold file | **선택** | Kit별 threshold, reporting용 metadata |

**중요**:
- 2와 4가 없으면 1, 3의 cell color annotation에서 status를 parsing (fallback, tidyxl 사용)
- 5는 optional. 없어도 pipeline 전체 동작. 있으면 reporting/vignette에서 활용.

(→ decisions.md 참조)

### 2.2 Concentration File 포맷
- WebIDQ export xlsx
- **Input xlsx는 시트 1개여야 함**. Sheet 이름은 검증하지 않으며, reader는 항상 첫 시트를 읽는다. Vendor가 부여한 sheet 이름이 무엇이든(`Conc_raw data`, `Conc` 등) 그대로 입력 가능.
- **Header row**: `Sample identification`, `Sample description`, `Submission name`, `Collection date`, `Species`, `Material`, 이후 metabolite short names (1231개, Quant 1000 기준)
- **Data rows**: 샘플 정보 (6개 메타 컬럼) + 측정값
- **Sample file과 QC file은 분리해서 input**
- Kit 정보는 `Submission name`에서 regex `KIT\d+`로 추출. NA submission name 은 NA kit_id 로 처리
- **Vendor 가 leading/trailing 빈 행을 export 하는 경우 reader 가 자동 drop** (`Sample identification` 이 NA 인 row)
- **Target Normalization [Median QC Level 2]이 이미 적용됨** → 재normalization 불필요

### 2.3 Status File 포맷 (Primary)
- Concentration file과 **행/열 레이아웃이 동일**. 단, sheet 이름은 vendor가 conc 파일과 다른 이름을 부여하는 경우가 있으며, reader는 sheet 이름을 검증하지 않고 첫 시트를 읽는다.
- 차이: metabolite 컬럼의 셀 element가 numeric 대신 **status text string**
- 가능한 status text:
  - `Valid`
  - `< threshold` (= valid value로 취급)
  - `< LLOQ`
  - `< LOD`
  - `> ULOQ`
- **Export에 포함되지 않는 status** (WebIDQ UI에만 있음, 제외):
  - `Accuracy out of range`, `Blank out of range`, `CV out of range`, `ISTD out of range`, `No intercept`
- Missing/NA는 실제로 거의 발생하지 않으나 방어적으로 handling

### 2.4 Status Color Annotation (Fallback)
- Sample/QC concentration file 자체에 cell fill color로 status가 인코딩되어 있음
- Color → status 매핑 테이블 (WebIDQ UI 확인값):

| Status | Hex |
|---|---|
| Valid | `#B9DE83` |
| `< threshold` | `#BBA7B9` |
| `< LOD` | `#A28BA3` |
| `< LLOQ` | `#B2D1DC` |
| `> ULOQ` | `#7FB2C5` |

- 한 셀에는 한 색만 존재
- tidyxl로 parsing. 이 매핑은 `sysdata.rda` 내장.
- **Cross-validation 옵션**: text와 color 둘 다 있으면 불일치 cell 리포트 (`validate_status()`)

### 2.5 LOD/LLOQ/ULOQ Threshold File (Optional)
- Per-kit xlsx, 4 rows × (1 + 1231) cols (vendor 실측)
- R1: `Measurement time` + metabolite short names
- R2: `LOD (calc.)` + 값 (전 metabolite)
- R3: `ULOQ` + 값 (일부 NA)
- R4: `LLOQ` + 값 (일부 NA)
- Sheet 이름은 검증하지 않음. 첫 시트만 읽음.
- **Kit ID 는 파일 내부에 없음** → 파일명에서 regex `KIT\d+` 로 추출. 파일명에 KIT 패턴이 없으면 명시 에러.
- Biocrates 측 설명:
  - LOD = `LOD (calc.)` primary, 0이면 `LOD (from OP)` fallback
  - LLOQ/ULOQ는 `from OP` 값만
  - Lipid metabolite에서 LLOQ/ULOQ NA는 WebIDQ 자체 이슈 (ROADMAP.md §Open questions 참조)
- **Status 재계산에 사용하지 않음**. Metadata 전용.

### 2.6 QC Sample 구조
- Pooled QC01/02/03 × KIT01/02 = 6 rows
- Pooled QC concentration file은 전체 1231 metabolite 커버
- Sample-level CV fail은 WebIDQ가 자체적으로 걸러줌 → input QC table에는 valid sample만 포함된다는 전제
- Vendor xlsx 가 데이터 끝/시작에 빈 행을 export 하는 경우가 있음 → reader 단계에서 `Sample identification` 이 NA 인 row 는 자동 drop

### 2.7 Metabolite 이름 규칙
- **Canonical name = short name 그대로** (괄호 포함)
- 예: `C3-DC (C4-OH)`, `C5-OH (C3-DC-M)`, `C6 (C4:1-DC)`
- WebIDQ export에서도 괄호 포함 full form으로 출력됨
- **이름 매칭 정책**: exact match only, no fuzzy
  - 1st pass: 그대로 비교
  - 2nd pass: whitespace collapse, `Total ` prefix 제거 후 비교
  - 여전히 unmatch → 명시적 에러

---

## 3. Status 처리 로직

### 3.1 Status 카테고리 (Canonical)
```
Valid          : 정상 정량
< threshold    : Valid로 취급
< LLOQ         : LLOQ 미만, 정량 신뢰도 낮으나 Valid로 취급 (원본 annotation은 보존)
< LOD          : LOD 미만, 검출 한계 이하
> ULOQ         : 상한 초과
Missing        : NA (방어적 logic용)
```

**원칙**: 원본 status 문자열은 assay에 그대로 보존하되, filter/flag 판정 시에는 `< threshold`와 `< LLOQ`를 Valid로 통합해 취급한다.

### 3.2 Helper 함수
```
is_valid(status)      → status ∈ {Valid, < threshold, < LLOQ}
is_below_lloq(s)      → status ∈ {< LLOQ}    # 원본 annotation 조회용, filter 로직엔 미사용
is_below_lod(s)       → status ∈ {< LOD}
is_above_uloq(s)      → status ∈ {> ULOQ}
is_missing(s)         → status ∈ {Missing, NA}
```

### 3.3 Zero 값 처리
- `x == 0`인 cell도 status annotation이 부여되어 있음 (대부분 `< LOD` 또는 `< LLOQ`)
- Zero는 별도 카테고리로 두지 않음. Status annotation 그대로 사용.
- 모든 measurement 값은 **imputation 없이 그대로 유지** (Inf/NaN/0 포함)

### 3.4 ULOQ 값 처리
- Default: 그대로 사용
- Option `winsorize_uloq = TRUE`: `> ULOQ` 값을 ULOQ 값으로 대체
  - Threshold file이 있어야 동작. 없으면 warning + skip.

---

## 4. 전처리 (Preprocessing) 파이프라인

### 4.1 Pipeline 순서

```
1. Import
   - sample_conc, sample_status → ConcentrationData (assay + status)
   - qc_conc, qc_status → QCData (assay + status)
   - (optional) threshold_file → ThresholdData

2. Sample-level filtering
   - <LOD rate 기반 outlier detection
   - (optional) >ULOQ rate 추가 고려

3. Metabolite-level filtering (순차 적용)
   - 3a. CV filter (per-kit AND primary / pooled option)
   - 3b. LOD rate filter (default 50%, 조정/disable 가능)

4. QC diagnostics
   - PCA (pooled QC only, pooled QC + samples)
   - PERMANOVA (kit effect)
   - CV distribution plot, LOD rate distribution plot, filter summary

5. (Log2 transform: 통계 layer 진입 시점에만, 원본 불변)
```

순서 근거: sample → metabolite

### 4.2 Sample-level Filtering
- 각 sample별 `<LOD rate` 계산
  - **분모**: 전체 metabolite 수 (Missing 포함)
  - **분자**: status가 `< LOD`인 metabolite 수
- Outlier detection method:
  - `"IQR"` (default): Q3 + 1.5×IQR 초과
  - `"percentile"`: 95th percentile 초과 (`percentile_cutoff = 0.95` 조정 가능)
  - `"both"`: OR
- `outlier_metric` 인자:
  - `"below_lod"` (default): `<LOD rate` 기반만
  - `"above_uloq"`: `>ULOQ rate` 기반만
  - `"either"`: 둘 중 하나라도 outlier이면 outlier

### 4.3 Metabolite-level Filtering

#### 4.3.1 CV filter (Step 3a)
- QC concentration + QC status에서 metabolite별 CV(%) 재계산
  - Raw QC 값 그대로 사용 (status가 `<LOD`인 cell도 포함 — Biocrates 공식 방법론)
- `cv_mode`:
  - `"per_kit"` (default): KIT01 3 replicates CV + KIT02 3 replicates CV 각각 계산 → **둘 다** threshold 이하여야 pass (AND)
  - `"pooled"`: 6 replicates 전체로 CV 계산 → threshold 이하면 pass
- Default threshold: 20%
- WebIDQ export에 들어있는 pre-computed CV row는 metadata로만 보관 (자체 재계산 값과 다를 수 있음)

#### 4.3.2 LOD rate filter (Step 3b, 조건부)
- `apply_lod_rate_filter = TRUE` (default)
- CV filter를 **통과한** metabolite 중에서, `<LOD rate > threshold`인 것 추가 제거
- Default threshold: 50%
- **분모**: 전체 샘플 수 (Missing 포함)
- `apply_lod_rate_filter = FALSE`면 이 단계 건너뜀 → Biocrates 공식 프로토콜(CV only) 엄격 준수
- Filter 순서 근거: CV first, then LOD rate

### 4.4 QC Diagnostics

모두 별도 함수, pipeline에서 자동 호출 안 함. User가 명시적 호출.

- `qc_pca_plot(prep, mode = c("qc_only", "with_samples"))` — PCA score plot (ggplot)
- `qc_permanova(prep)` — kit effect 정량 (vegan::adonis2)
- `qc_cv_plot(prep)` — metabolite별 CV 분포 histogram + threshold line
- `qc_lod_rate_plot(prep, by = c("sample", "metabolite"))` — LOD rate 분포
- `qc_filter_summary(prep)` — pass/fail 요약표 (data.frame)

### 4.5 Log2 Transform
- `log2(x + pseudocount)`
- Pseudocount: **데이터 왜곡 최소화 + 계산 오류 방지** 수준
  - metabolite-wise non-zero minimum / 2
  - 실제 0이 존재하는 경우에만 적용, 0이 없으면 pseudocount = 0
- 적용 시점: **통계 layer 진입 직전**에만. Raw 원본은 언제나 불변.
- 동일 룰을 사용자 분석 스크립트에서 재사용할 수 있도록 `log2_with_pseudocount()`로 export

---

## 5. 내장 데이터 (sysdata.rda)

### 5.1 데이터 원본
- **Metabolite reference**: `biocrates_Quant1000_metabolites.tsv` (1234 metabolites, 49 classes)
- **Pathway set source**: `data-raw/pathway_sets/<domain>/` 내 master TSV (§8.5)
- **Status color map**: §2.4 hex 매핑 하드코딩

### 5.2 내장 데이터 객체
```r
metabolite_dict      # data.frame: short_name, long_name, class,
                     #   kegg_id, hmdb_id, lion_id,
                     #   pubchem_id, chebi_id, refmet_id, lipidmaps_id (1234 rows)
                     #   *_id (kegg/hmdb/pubchem/chebi/refmet/lipidmaps): RaMP-DB
                     #     cross-refs (pathview + 소분자 set join + 외부 도구 연동)
                     #   lion_id: 매칭된 LION 정규 lipid 이름 (lipid; lion set join)
                     #   미매칭은 NA. (→ decisions.md #4)
status_color_map     # data.frame: hex → status_text
pathway_sets         # nested named list: [domain][set_name] → character(metabolite_ids)
                     #   domain ∈ {"immunomet",
                     #             "reactome", "wikipathways", "smpdb", "lion",
                     #             "source", "health"} (확장 가능)
                     #   모든 도메인 metabolite-level
pathway_sets_meta    # data.frame: domain, set_name, member, direction, brief_note
                     #   long-format, per-member metadata 유지
```

**Pathway set data source**: `data-raw/pathway_sets/<domain>/` 내 master TSV 파일들. 빌드 상세는 §8.5 참조.

---

## 6. Statistical Layer

### 6.1 Approach A: Metabolite-as-Feature
- 각 metabolite를 individual feature로 취급 (`test_metabolites()`)
- 입력: `PreprocessedData` (metabolite 농도 행렬). Outlier sample은 제외.
- Group comparison: moderated t (default) 또는 Wilcoxon. Moderated-t 엔진은 R `limma`, Python 포트는 `inmoose.limma` (empirical Bayes 방법·default 지위 불변, backend만 상이 → decisions.md #7)
- Multiple testing correction: BH + BY 동시 리포트 (default)
- **주의**: 같은 chemical class 내 metabolite 간 상관이 높을 수 있음
- Output: data.frame (metabolite 1행), 컬럼 `metabolite`, `log2FC`, `P.Value`, `padj_*` 등

### 6.2 Approach B: GSEA over curated pathway sets
- 목적: 생물학적 pathway/term 수준 해석
- `gsea_pathway_sets(prep, group, domain, ...)` — `PreprocessedData`의 metabolite 농도에서
  group log2FC rank를 계산하고 `fgsea::fgseaMultilevel`로 enrichment 검정
- Set source: `sysdata.rda`의 `pathway_sets` (모두 metabolite-level)
- 지원 도메인 (전체 목록·상세는 §8.5):
  - In-house 큐레이션: `immunomet`
  - 외부 표준 DB 파생: `reactome`·`wikipathways`·`smpdb` (소분자, RaMP-DB pathway), `lion` (lipid, LION functional term), `source` (RaMP origin — Plant/Microbe), `health` (RaMP 질병-연관, disease-MSEA) (→ decisions.md #4)
- **GSEA universe 분리**: `smpdb`(소분자)와 `lion`(lipid)은 다루는 metabolite 종류가 달라 한 rank 벡터에 섞지 않고 별도 호출 권장 (fgsea normalization 일관성). 패널에 멤버가 전혀 겹치지 않으면 명시적 no-overlap 에러
- Direction 정보(`pathway_sets_meta`의 `direction`)는 fgsea에 직접 전달하지 않고, 결과 해석 참고용으로 제공
- Output: fgsea result (data.table) — `gsea_to_df()`로 list-column을 collapsed string으로 변환 후 TSV 저장 가능

#### 6.2.1 Method 선택 논거
- GSEA backend: R 판본은 `fgsea::fgseaMultilevel` (method paper 있음 — Korotkevich et al. 2021), Python 포트는 `gseapy` prerank. 두 backend는 수치 동일하지 않으며, Python 포트는 R fgsea와의 bit-identical을 요구하지 않고 seed 고정으로 결정론만 보장한다 (→ decisions.md #8)
- `limma-voom` 불필요: count data가 아닌 continuous concentration이므로

### 6.3 결과 저장 정책
- A, B(GSEA), C(ORA)의 결과를 **하나로 합치지 않음**
- 각각 따로 반환, 따로 저장. User가 원하면 스스로 join

### 6.4 Pathway set loader API
```r
# 사용 가능한 set 목록
list_pathway_sets()
# → data.frame: domain, set_name, n_members

# fgsea에 바로 투입 가능한 named list 반환
get_pathway_sets(domain = "reactome")
# → list(REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE = c("Citric acid", ...), ...)

# Metadata 포함 long-format 반환
get_pathway_sets_meta(domain = "immunomet")
# → data.frame: domain, set_name, member, direction, brief_note

# 외부 툴을 위한 GMT 파일 경로
pathway_set_gmt_path(domain = "immunomet")
# → system.file("extdata/pathway_sets/immunomet_metabolites.gmt", package = "MetaboSetR")
```

### 6.5 External pathway map 시각화 (pathview)
차등 metabolite 를 KEGG pathway map 위에 올리는 시각화는 **패키지가 wrapping 하지 않고 사용자가 `pathview` 를 직접 호출**한다. (→ decisions.md #4)

- **입력**: `annotate_metabolites(x)` 로 얻은 `kegg_id` (소분자 KEGG compound ID) 를 key 로 한 log2FC 등 named vector. `annotate_metabolites(x)` 는 short_name 벡터를 `metabolite_dict` 의 id 컬럼에 매핑한 data.frame 을 입력 순서대로 반환 (미매칭 NA)
- **KEGG embedding 안 함**: KEGG pathway set/map 은 라이선스가 bulk 재배포를 금지하므로 패키지에 포함하지 않는다. `kegg_id` 는 RaMP-DB 의 cross-ref (HMDB/WikiPathways 경유) 에서 얻은 것이라 license-clean
- **비의존**: `pathview` 는 패키지 의존성이 아니다. 사용자가 별도 설치 후 호출
```r
# 사용자 분석 스크립트 예시 (패키지 외부)
fc_by_kegg <- ...   # named numeric: names = kegg_id, value = log2FC
pathview::pathview(cpd.data = fc_by_kegg, pathway.id = "hsa00020", species = "hsa")
```
- Lipid species 는 KEGG map 대상이 아니다 (§8.5의 lion 도메인 functional enrichment 로 대체)

### 6.6 Approach C: ORA over curated pathway sets (threshold-based)
- 목적: threshold 로 고른 유의미 metabolite 집합이 어떤 curated pathway set 에 over-represent 되는지 검정. GSEA(B, rank-based) 의 상보적 접근 (→ decisions.md #3)
- Set source: GSEA 와 동일한 `sysdata.rda` 의 `pathway_sets` (동일 도메인·identifier space)
- 입력: `sig` (유의미 metabolite 이름), `universe` (측정 panel 전체 metabolite 이름)
- **Universe = 측정 panel (Option A)**: 배경집합을 사용자가 측정한 metabolite 전체로 잡는다. 측정하지 않은 metabolite 는 유의미해질 수 없으므로 배경에서 제외하는 것이 고정 패널에서 정직. 내장 세트가 panel-driven 큐레이션이라 정합 (→ decisions.md #3)
- 검정: pathway 별 hypergeometric (`stats::phyper`, over-representation 단측). 각 pathway 멤버는 universe 와 교집합한 뒤 계산
- Multiple testing: `test_metabolites()` 와 동일하게 BH + BY 동시 리포트 (default)
- Direction 정보는 GSEA 와 동일하게 검정에 전달하지 않음 (해석용; `get_pathway_sets_meta()` 로 join)
- Output: data.frame (pathway 별 1행) — `set_name`, `set_size`, `overlap`, `expected`, `pval`, `padj_*`, `overlap_members`

---

## 7. Package Architecture

### 7.1 디렉토리 구조
```
MetaboSetR/
├── R/
│   ├── data.R                 # 내장 데이터 documentation
│   ├── read-webidq.R          # WebIDQ concentration + status reader
│   ├── read-status-color.R    # color annotation parser (tidyxl, fallback)
│   ├── read-thresholds.R      # LOD/LLOQ/ULOQ reader (optional)
│   ├── status-helpers.R       # is_valid / is_below_lod / ... helpers
│   ├── ConcentrationData.R    # S4 class
│   ├── QCData.R               # S4 class
│   ├── ThresholdData.R        # S4 class
│   ├── PreprocessedData.R     # S4 class
│   ├── preprocess.R           # pipeline (sample filter, metabolite filter)
│   ├── filter-samples.R       # outlier detection
│   ├── filter-metabolites.R   # CV + LOD rate filter
│   ├── qc-diagnostics.R       # PCA, PERMANOVA, CV/LOD plots, summary
│   ├── annotate.R             # annotate_metabolites (id cross-ref accessor)
│   ├── stats-feature.R        # Approach A: test_metabolites (limma/wilcoxon)
│   ├── stats-gsea.R           # Approach B: gsea_pathway_sets (fgsea)
│   ├── stats-ora.R            # Approach C: ora_pathway_sets (hypergeometric)
│   ├── pathway-sets.R         # loader + list/get functions
│   ├── transform.R            # log2 + pseudocount
│   └── utils.R
├── data-raw/                  # 원본 TSV → sysdata 생성 스크립트
│   ├── reference/
│   │   └── biocrates_Quant1000_metabolites.tsv
│   ├── pathway_sets/
│   │   ├── immunomet/         # O'Neill 2016, metabolite-level
│   │   ├── ramp/              # reactome/wikipathways/smpdb/source/health + id annotation
│   │   ├── lion/              # LION lipid functional terms
│   │   └── README.md          # 도메인 추가 가이드
│   └── make_sysdata.R         # builds R/sysdata.rda
├── inst/extdata/
│   ├── pathway_sets/          # *_metabolites.gmt
│   └── synthetic/             # synthetic example data
├── tests/testthat/
├── vignettes/
└── DESCRIPTION
```

### 7.2 S4 Class 계층
- `ConcentrationData`, `QCData`, `ThresholdData`, `PreprocessedData` — plain S4 class
- `PreprocessedData`는 sample/QC 농도 + 필터 결과(sample outlier, metabolite pass)를 보유. 통계 함수의 입력 단위.

### 7.3 Main Pipeline API

```r
# 1. Import
conc <- read_webidq(
  concentration_file = "sample_conc.xlsx",
  status_file        = "sample_status.xlsx"
)
qc <- read_webidq_qc(
  concentration_file = "qc_conc.xlsx",
  status_file        = "qc_status.xlsx"
)
thresh <- read_thresholds(c("lod_kit01.xlsx", "lod_kit02.xlsx"))  # optional

# 2. Preprocess
prep <- preprocess(
  conc, qc,
  thresholds            = thresh,             # optional
  sample_outlier_method = "IQR",
  outlier_metric        = "below_lod",
  cv_threshold          = 20,
  cv_mode               = "per_kit",
  lod_rate_threshold    = 0.5,
  apply_lod_rate_filter = TRUE,
  winsorize_uloq        = FALSE
)

# 3. QC diagnostics (optional, user-initiated)
qc_pca_plot(prep, mode = "qc_only")
qc_permanova(prep)
qc_filter_summary(prep)

# 4. Metabolite-level testing
res_feat <- test_metabolites(
  prep, group = "treatment",
  method = "limma",         # or "wilcoxon"
  padjust = c("BH", "BY")
)

# 5. GSEA over curated pathway sets
res_gsea <- gsea_pathway_sets(prep, group = "treatment", domain = "reactome")

# 6. ORA over curated pathway sets (threshold-based)
#    universe = 측정한 전체 metabolite (background); sig = 유의미 metabolite.
res_ora <- ora_pathway_sets(
  sig      = res_feat$metabolite[res_feat$P.Value < 0.05],
  universe = res_feat$metabolite,
  domain   = "reactome"
)
# 결과는 분리 저장, 통합하지 않음

# 7. id annotation (KEGG-map overlay 등 외부 도구용)
ann <- annotate_metabolites(colnames(prep@sample@assay))
```

### 7.4 Dependencies
아래는 R 판본의 의존성이다. Python 포트의 대응 스택(numpy/pandas/openpyxl/scipy/
inmoose/gseapy/scikit-bio/matplotlib)은 `MIGRATION.md` §3.1 참조 (→ decisions.md #6).
- `readxl` — xlsx value reading
- `tidyxl` — cell color parsing (fallback status)
- `vegan` — PERMANOVA
- `ggplot2` — diagnostic plots
- `fgsea` — GSEA backend
- `limma` — moderated t-test
- Base: `methods`, `stats`, `utils`
- **Tidyverse dependency 의도적 제외**
- **`pathview` 의도적 비의존**: KEGG map 시각화는 사용자가 직접 호출. 패키지는 `kegg_id` annotation 만 제공
- RaMP-DB/LION 매칭·빌드는 `data-raw` 내 Python (sqlite3 stdlib) 으로 수행 (CLAUDE.md §5.4) — R 런타임 신규 의존성 없음

### 7.5 Testing 전략
- Unit tests: status reader (text + color), filter 함수, metabolite statistics, ORA hypergeometric (hand-computed ground truth)
- Integration test: synthetic data 전체 pipeline end-to-end
- Pathway set: loader/schema/snapshot + GMT parse

---

## 8. Reference Files

### 8.1 Biocrates Q1000 Metabolite List
- **Source**: `data-raw/reference/biocrates_Quant1000_metabolites.tsv`
- **구조**: analyte_class, shortname, fullname
- **총 1,234 metabolites**, 49 biochemical classes
- Small molecules: LC-MS/MS. Lipids: FIA-MS/MS.
- Short name ↔ long name 매핑 테이블로 사용

### 8.2 실제 Input 파일 포맷 참조
- `OMICS-262002_LOD__LLOQ__ULOQ_KIT{01,02}.xlsx` — threshold file
- `OMICS-262002_MxPQuant1000_Human_plasma_v1_QC.xlsx` — QC concentration

### 8.5 Pathway Sets

#### 8.5.1 위치 및 파일 구조
- **Source of truth**: `data-raw/pathway_sets/<domain>/` 하위의 master TSV 파일 (모두 metabolite-level)
- **GMT 사본**: 각 `build.py` 실행 시 해당 domain 디렉토리에 생성 + `inst/extdata/pathway_sets/`에 동기화
- **sysdata 객체**: `make_sysdata.R`가 master TSV를 읽어 `pathway_sets`, `pathway_sets_meta`로 변환하여 `R/sysdata.rda`에 저장
- **예외 (RaMP-derived)**: `reactome`·`wikipathways`·`smpdb`·`source`·`health` 도메인은 공유 소스(RaMP-DB)라 `data-raw/pathway_sets/ramp/`의 단일 `build.py`가 master TSV + metabolite ID annotation을 함께 생성 (per-domain 디렉토리 컨벤션의 예외)

#### 8.5.2 Domain: `immunomet`
- **방법론**: Top-down, reference-based (O'Neill 2016 framework)
- **Primary reference**: O'Neill LAJ, Kishton RJ, Rathmell J. (2016). *Nat Rev Immunol* 16(9):553–565
- **Scope**: O'Neill의 6개 core pathway 중 5개 구현 (pentose phosphate pathway는 Biocrates Q1000 panel이 canonical metabolite를 커버하지 않아 제외)
- **Set 개수**: 5 (GLYCOLYSIS, TCA, FAO, FAS, AA_METABOLISM)
- **Level**: metabolite (멤버 = metabolite shortname)
- **Direction**: master TSV에 direction 컬럼 없음 → NA (면역세포 상태는 bidirectional)
- **Overlap 정책**: 생물학적으로 의미있는 overlap 허용. 데이터 구조상 artifact overlap은 명시적 제외 (예: TCA 중간체는 Dicarboxylic acids class에 속하지만 FAO에서 제외)
- **Set 크기 (metabolite)**:

| Set | Metabolites |
|---|---:|
| IMMUNOMET_ONEIL_GLYCOLYSIS | 3 |
| IMMUNOMET_ONEIL_TCA | 10 |
| IMMUNOMET_ONEIL_FAO | 58 |
| IMMUNOMET_ONEIL_FAS | 445 |
| IMMUNOMET_ONEIL_AA_METABOLISM | 138 |

#### 8.5.6 Domains: `reactome`, `wikipathways`, `smpdb`, `source`, `health` (RaMP-derived)
- **방법론**: 외부 표준 DB (RaMP-DB v3.0.7) 파생 — in-house 큐레이션 아님 (→ decisions.md #4)
- **Source**: RaMP-DB SQLite dump (ncats/RaMP-DB). Raw (1.95GB) 는 `.gitignore`, 파생 master TSV/GMT 만 commit. 도메인들은 `data-raw/pathway_sets/ramp/` 의 단일 build.py 가 함께 생성
- **Level**: metabolite only (소분자 — lipid 은 `lion` 도메인)
- **출처 분리 (RaMP `pathway.type` 으로 필터)**:
  - `reactome` ← type=reactome (CC-BY 4.0, 인용)
  - `wikipathways` ← type=wiki (WikiPathways CC0)
  - `smpdb` ← type=hmdb (SMPDB/Wishart 학술·개인 허용 + 인용)
  - 제외: kegg-type (KEGG 라이선스), pfocr (figure OCR 텍스트마이닝, 품질 낮음)
- **ID 매칭 / annotation**: Biocrates short_name·fullname ↔ RaMP `analytesynonym` (NOCASE) → RaMP `source` 테이블 cross-ref. `kegg_id`/`hmdb_id` 에 더해 `pubchem_id`/`chebi_id`(`CHEBI:` 표준형)/`refmet_id`/`lipidmaps_id` 부여 — metabolite 당 대표 1개. 자동 1차 + 미매칭 수동 검수, `assignment_log.tsv`
- **Set 선택 (panel-driven)**: 각 도메인에서 pathway 의 metabolite 멤버를 Biocrates 매칭 패널과 교집합 → 멤버 ≥3 인 set 만 유지 + 동일 멤버십 dedupe. 전량 사용은 중복 + 다중검정 폭발이라 회피
- **Direction**: 전부 `-` (hypothesis-free 표준 set)
- **`source` 도메인**: RaMP ontology (`HMDBOntologyType='Source'`) 에서 파생. 2 set — `SOURCE_PLANT`, `SOURCE_MICROBE`. ⚠️ source-NA 는 host 가 아니라 "기원 미주석" (2차 담즙산·IPA·equol 등 미생물 대사체도 NA)
- **`health` 도메인**: `HMDBOntologyType='Health condition'` 에서 파생한 질병-연관 metabolite set (disease-MSEA 식). ~51 set (`HEALTH_*`). ⚠️ **Caveat**: HMDB 질병 연관은 문헌 기반이라 (a) 연구량 편향, (b) 동일 질병 contrast 분석 시 순환논리, (c) "association" 이지 pathway 아님 → 결과는 가설 수준

#### 8.5.7 Domain: `lion` (external DB)
- **방법론**: 외부 lipid functional ontology (LION) 파생
- **Source**: `all-LION-lipid-associations.csv` (lipidontology.com). lipid × term 이진 행렬. Raw 다운로드는 `.gitignore`
- **Primary reference**: Molenaar MR, et al. (2019). LION/web. *GigaScience* 8(6):giz061
- **Level**: metabolite only (lipid species)
- **Set 구성**: LION term → 해당 term 에 속한 lipid species. Biocrates lipid 가 매칭된 term 만 set 으로 유지 (panel-driven ≥3+dedupe)
- **Member 식별**: Biocrates lipid 표기 ↔ LION 정규 lipid 이름 매칭. 매칭된 LION 정규 이름은 `metabolite_dict$lion_id` 에도 저장
- **Direction**: 전부 `-`

#### 8.5.8 Naming convention
- 모든 set name은 대문자 + 언더스코어, domain prefix 필수: `IMMUNOMET_*`, `REACTOME_*`, `WIKIPATHWAYS_*`, `SMPDB_*`, `LION_*`, `SOURCE_*`, `HEALTH_*`
- Reference-based domain은 reference 약어도 포함: `IMMUNOMET_ONEIL_*`
- 외부 DB 파생 도메인은 원본 pathway/term 이름을 대문자 + 언더스코어로 정규화하여 prefix 뒤에 붙임 (예: `SMPDB_CITRIC_ACID_CYCLE`, `LION_MEMBRANE_FLUIDITY`)

#### 8.5.9 Master TSV schema
- **Metabolite-level** (`<domain>_metabolites_master.tsv`):
  - `shortname`, `fullname`, `analyte_class`: metabolite 식별자
  - `set_name`: 해당 set name
  - `match_basis`: `individual` 또는 `class:<class_name>`
- **External-DB metabolite-level** (`reactome_master.tsv` 등): 위 base 에 DB-specific evidence 컬럼 추가 (`match_basis`, `kegg_id`/`hmdb_id`/`ramp_id`/`pathway_source_id` 또는 `ontology_id`, `lion_id`/`lion_term_id` 등)

#### 8.5.10 Reproducibility
- 각 domain 디렉토리의 `build.py`를 실행하면 master TSV + GMT가 재생성됨
- Build의 deterministic: 동일 reference TSV 입력 → 동일 output
- Audit trail: `assignment_log*.tsv`에 모든 포함/제외 판정 기록

---

## Appendix: 관련 문서
- `ROADMAP.md` — 버전 로드맵, 현재 상태, open questions
- `decisions.md` — Key design decisions log (append-only)
- `CLAUDE.md` — Claude Code 작업 규칙

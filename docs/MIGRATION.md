# MetaboSetR — Python 마이그레이션 문서

> **Status**: **전 계층 이식 완료** (`python/metabosetr/`). 전처리 → 통계(A/B/C) → 로더/annotate/QC까지 동작, pytest 36 케이스 통과. 남은 것은 문서 정식화(decisions.md/DESIGN.md 반영, LICENSE 교체)와 커밋.
> **Purpose**: R 패키지 MetaboSetR를 순수 Python으로 재설계하기 위한 파일-단위 매핑과 이식 순서.
> **권위**: 이 문서는 DESIGN.md에 종속된다. 설계 자체를 바꾸는 결정은 먼저 `decisions.md` → `DESIGN.md`에 반영된 뒤 여기에 내려온다 (CLAUDE.md §1).
> **우선순위**: 전처리(preprocessing) 파이프라인이 최우선이었고(완료), 이어서 통계/GSEA/ORA/QC까지 이식 완료.

## 0. 이식 현황 (2026-08-25 기준)

| 계층 | Python 모듈 | 상태 | R 대비 검증 |
|---|---|---|---|
| 전처리 (reader/filter/transform/classes) | `read_webidq`, `read_status_color`, `read_thresholds`, `filter_*`, `transform`, `classes`, `preprocess` | ✅ | CV/LOD/IQR 필터 R base 함수와 ~1e-14 일치; text↔color 교차검증 |
| A. feature (default) | `stats_feature` (InMoose moderated-t) | ✅ | InMoose drop-in (R limma near-identical) |
| A. feature (support) | `stats_feature` (SciPy Wilcoxon) | ✅ | R `wilcox.test` p-value·log2FC ~1e-16 일치 |
| B. GSEA | `stats_gsea` (gseapy prerank) | ✅ | fgsea와 상이(설계); seed 결정론 |
| C. ORA | `stats_ora` (SciPy hypergeom) | ✅ | R `phyper` 일치 |
| 다중검정 | `_padjust` (BH/BY) | ✅ | R `p.adjust` 일치 (NA 처리 포함) |
| 로더/annotate | `pathways`, `annotate`, `reference` | ✅ | `make_sysdata.R` assert 카운트 전부 일치 |
| QC diagnostics | `qc_diagnostics` (PCA/PERMANOVA/plots) | ✅ (best-effort) | 스모크 |

남은 작업: (1) decisions.md 결정 A/B/C append + DESIGN.md 반영, (2) LICENSE → GPL-3.0, (3) 커밋.

---

## 1. 목표 · 범위 · 비목표

### 1.1 목표
- MetaboSetR의 기능을 **순수 Python**으로 재구현. R 런타임 의존을 완전히 제거 (rpy2 미사용).
- 최우선: **WebIDQ reader → 전처리 필터 → PreprocessedData 산출**까지의 파이프라인.
- 후속: metabolite 통계(A), GSEA(B), ORA(C), QC diagnostics, annotation.

### 1.2 범위
- R `R/` 아래의 런타임 로직 전체.
- `data-raw/`의 sysdata 빌드 (R `make_sysdata.R` → Python 스크립트, 산출물은 pickle/parquet).

### 1.3 비목표
- **Indicator 레이어 이식 안 함** — 애초에 존재하지 않음 (decisions.md #1에서 제거됨).
- **pathview wrapping 안 함** — DESIGN §6.5대로 사용자가 직접 호출. Python에서도 `kegg_id` annotation만 제공.
- `data-raw/pathway_sets/<domain>/build.py` **재작성 안 함** — 이미 Python(stdlib). 그대로 재사용.

---

## 2. 확정된 설계 결정 (decisions.md 반영 대상)

이 세 건은 대화에서 확정됐고, 실제 이식 착수 전 `decisions.md`에 신규 항목으로 append + `DESIGN.md` 본문에 반영해야 한다.

| # | 결정 | 영향 절 |
|---|---|---|
| A | **Approach A 엔진: R `limma` → `inmoose.limma`** (Python 포트). 방법(empirical Bayes moderated-t)·default 지위 **불변**. `inmoose`가 R limma의 drop-in(near-identical)이라 수치 재현성 유지. Wilcoxon(`scipy`)은 support 옵션으로 유지. | DESIGN §6.1, §7.4 |
| B | **GSEA backend: R `fgsea` → `gseapy`(prerank) 또는 `blitzgsea`.** R fgsea와 수치 동일성은 포기(수용). 재현성은 **seed 고정으로 결정론적 보장**. API에 seed 인자 노출. | DESIGN §6.2, §6.2.1, §7.4 |
| C | **라이선스: proprietary → GPL-3.0 (copyleft).** `inmoose`가 GPL-3.0 핵심 의존이라 정합. Biocrates indicator 카탈로그는 패키지에 없음(검증 완료 — §7)이라 copyleft blocker 없음. 코드=GPL-3.0, **번들 데이터 라이선스는 별도 명시** 필요(§7). | DESIGN §1.2, LICENSE, DESCRIPTION |

---

## 3. 의존성 매핑 (R → Python)

| R 패키지/함수 | Python 대응 | 사용 계층 | 재현성 |
|---|---|---|---|
| `readxl::read_excel` | `openpyxl` (또는 `pandas.read_excel(engine="openpyxl")`) | reader | 동일 |
| `tidyxl::xlsx_cells`/`xlsx_formats` (cell fill 색) | `openpyxl` `cell.fill.fgColor.rgb` | color 상태 fallback | 동일 (오히려 단순) |
| `methods::setClass` (S4) | `@dataclass` + `__post_init__` validity | 데이터 클래스 | — |
| `stats::quantile` | `numpy.quantile` (기본 linear = R type 7) | 이상치 검출 | **동일** (아래 §6.1 주의) |
| `stats::sd` | `numpy.std(ddof=1)` | CV 계산 | **주의**: numpy 기본 ddof=0. 반드시 ddof=1 |
| `rowSums`/`colSums(na.rm=TRUE)` | `numpy.nansum` / `(arr == x).sum(axis=)` | 필터 | 동일 |
| `stats::prcomp(center, scale.)` | `sklearn.decomposition.PCA` (사전 표준화) 또는 numpy SVD | QC PCA | 근접 (부호 임의성) |
| `vegan::adonis2` | `skbio.stats.distance.permanova` | QC PERMANOVA | 근접 (seed) |
| `stats::dist` (Euclidean) | `scipy.spatial.distance.pdist` | QC | 동일 |
| `ggplot2` | `matplotlib` (또는 `plotnine`) | 진단 플롯 | N/A |
| `stats::phyper` | `scipy.stats.hypergeom.sf` | ORA | **완전 동일** |
| `stats::p.adjust` (BH/BY) | `statsmodels.stats.multitest.multipletests` (`fdr_bh`, `fdr_by`) | 다중검정 | 동일 |
| `limma::lmFit`/`eBayes`/`topTable` | `inmoose.limma` 동명 함수 | Approach A | near-identical |
| `stats::wilcox.test` | `scipy.stats.mannwhitneyu` | Approach A (옵션) | 동일 |
| `fgsea::fgseaMultilevel` | `gseapy.prerank` / `blitzgsea` | Approach B | 다름 (seed 결정론) |
| `usethis::use_data(internal=TRUE)` → `sysdata.rda` | pickle / parquet + `metabosetr/_data.py` 로더 | 내장 데이터 | 동일 |

### 3.1 런타임 의존성 (예상 최소셋)
- **필수**: `numpy`, `pandas`, `scipy`, `openpyxl`, `statsmodels`
- **통계 A**: `inmoose` (GPL-3.0)
- **GSEA B**: `gseapy` (또는 `blitzgsea`)
- **QC**: `scikit-bio` (PERMANOVA), `matplotlib` (플롯)
- R판의 "tidyverse 배제" 원칙에 대응: **무거운 프레임워크(예: PyMC) 배제, scientific-core만.**

---

## 4. 계층별 재현성 등급

| 계층 | Python 엔진 | R 대비 |
|---|---|---|
| reader (value + color), 이름 해석, threshold reader | openpyxl / pandas | **동일** |
| status helpers | numpy/pandas 벡터화 | **동일** |
| sample outlier / CV / LOD-rate 필터 | numpy/scipy | **동일** (§6.1 ddof·quantile 주의 지키면) |
| log2 + pseudocount | numpy | **동일** |
| ORA (C) | scipy hypergeom | **완전 동일** |
| feature 통계 A (default) | inmoose moderated-t | **near-identical** (drop-in) |
| feature 통계 A (support) | scipy Wilcoxon | 동일 |
| GSEA (B) | gseapy/blitzgsea | 다름, seed로 결정론 |
| QC PCA / PERMANOVA | sklearn / scikit-bio | 근접 (부호·permutation seed) |

**요지**: 전처리 계층 전체는 R과 수치 동일하게 이식 가능하다. 수치가 달라지는 곳은 GSEA 한 곳뿐.

---

## 5. S4 → Python 클래스 매핑

R 4개 S4 class를 `@dataclass`로 옮긴다. Layout 규약은 **samples × metabolites** 유지 (WebIDQ export와 동일; 통계 계층 진입 시에만 전치).

| R S4 class | Python | slot/field |
|---|---|---|
| `ConcentrationData` | `ConcentrationData` (dataclass) | `assay: np.ndarray[float]` (sample×metab), `status: np.ndarray[str]` (동형), `sample_meta: pd.DataFrame`, `metabolite_names: list[str]` |
| `QCData` | `QCData` | 동일 구조 (qc_sample×metab) |
| `ThresholdData` | `ThresholdData` | `lod/lloq/uloq: np.ndarray[float]` (kit×metab), `kit_ids: list[str]`, `metabolite_names: list[str]` |
| `PreprocessedData` | `PreprocessedData` | `sample: ConcentrationData`, `qc: QCData`, `thresholds: ThresholdData | None`, `sample_outliers: np.ndarray[bool]`, `metabolite_pass: np.ndarray[bool]`, `sample_filter_log: pd.DataFrame`, `metabolite_filter_log: pd.DataFrame`, `params: dict` |

- R의 `validity=` 함수 → dataclass `__post_init__`에서 차원/타입 검증 (동일 assertion).
- R의 `setMethod("show", ...)` → `__repr__`.
- 생성자 헬퍼(`ConcentrationData(...)`)는 dataclass 기본 생성자로 대체.

---

## 6. 전처리 이식 순서 (Phase 1 — 최우선)

각 step은 검증 가능한 성공 기준을 갖는다. R 함수와 **동일 입력 → 동일 출력**이 기본 기준.

### 6.0 사전 준비
- `R/sysdata.rda`의 `metabolite_dict`, `status_color_map`를 Python이 읽을 형식(parquet/pickle)으로 export.
  - `metabolite_dict`: 이름 해석(`.resolve_metabolite_names`)에 `short_name`/`long_name` 필요.
  - `status_color_map`: color fallback(`.hex_to_status`)에 `hex`/`status_text` 필요.
  - 임시로 R에서 1회 export하거나, `data-raw`에서 Python으로 재빌드.

### 6.1 매핑할 R 파일 → Python 모듈 (전처리)

| R 파일 | Python 모듈(제안) | 핵심 함수 |
|---|---|---|
| `R/ConcentrationData.R`, `QCData.R`, `ThresholdData.R`, `PreprocessedData.R` | `metabosetr/classes.py` | 4 dataclass |
| `R/status-helpers.R` | `metabosetr/status.py` | `is_valid`, `is_below_lloq`, `is_below_lod`, `is_above_uloq`, `is_missing` |
| `R/read-webidq.R` | `metabosetr/read_webidq.py` | `read_webidq`, `read_webidq_qc`, 내부 `_read_conc_sheet`, `_read_status_text`, `_extract_kit_id`, `_resolve_metabolite_names`, `_is_padding_row` |
| `R/read-status-color.R` | `metabosetr/read_status_color.py` | `_read_status_color`, `_hex_to_status`, (`validate_status`는 R에서도 미구현) |
| `R/read-thresholds.R` | `metabosetr/read_thresholds.py` | `read_thresholds`, `_read_one_threshold_file`, `_kit_id_from_filename` |
| `R/filter-samples.R` | `metabosetr/filter_samples.py` | `_detect_sample_outliers`, `_apply_outlier_method` |
| `R/filter-metabolites.R` | `metabosetr/filter_metabolites.py` | `_cv_filter`, `_cv_pct`, `_lod_rate_filter` |
| `R/transform.py` ← `R/transform.R` | `metabosetr/transform.py` | `log2_with_pseudocount`, `_log2_one_row` |
| `R/preprocess.R` | `metabosetr/preprocess.py` | `preprocess`, `_winsorize_uloq` |

### 6.2 이식 step (순서)

1. **classes.py** → verify: 4 dataclass 생성 + validity assertion이 R validity와 동일 케이스에서 에러.
2. **status.py** → verify: hand-computed 케이스 (Valid/< threshold/< LLOQ→valid, NA 전파 규칙 §3.2) 통과.
3. **transform.py** → verify: R `log2_with_pseudocount` docstring 예제(low/high 행)와 bit-level 일치.
4. **read_webidq.py** (text 경로) → verify: `inst/extdata/synthetic/sample_{conc,status}.xlsx` 로드 결과가 R reader와 동일 shape·값. padding-row drop, kit_id 추출, 이름 3-pass 해석 포함.
5. **read_status_color.py** (color 경로) → verify: color fallback이 text 경로와 동일 status 행렬 산출 (synthetic 파일로 cross-check).
6. **read_thresholds.py** → verify: `thresholds_KIT01/02.xlsx` → `ThresholdData` 차원·값 일치.
7. **filter_samples.py** → verify: IQR/percentile/both × below_lod/above_uloq/either 조합이 R과 동일 outlier 벡터. **quantile type 7, IQR 공식 Q3+1.5·IQR**.
8. **filter_metabolites.py** → verify: per_kit(AND)·pooled CV, LOD-rate 필터가 R과 동일 pass 벡터. **CV = 100·sd(ddof=1)/mean, mean==0→NA**.
9. **preprocess.py** → verify: synthetic 전체를 통과시켜 `PreprocessedData` 산출, R `preprocess()` 결과와 `sample_outliers`/`metabolite_pass`/log 프레임 일치. winsorize_uloq 경로 포함.

### 6.3 이식 시 수치 함정 (반드시 지킬 것)
- **`sd`**: R `sd()`는 표본표준편차(n−1). numpy는 기본 ddof=0 → **`np.std(x, ddof=1)`**. CV 재현성의 핵심.
- **`quantile`**: R 기본 type 7 = numpy `np.quantile` 기본 linear. 일치하나, scipy 다른 함수 쓸 때 주의.
- **NA 처리**: R `na.rm=TRUE` → numpy `nan`-aware (`nansum`, boolean mask). status 비교 시 NA는 False로 세는 R 규약(`== "< LOD"`이 NA에서 FALSE) 재현.
- **status helper NA 전파**: `is_valid`는 NA→NA, `is_missing`은 NA→TRUE (§3.2). 벡터화 시 명시.
- **color reader ARGB**: openpyxl `fgColor.rgb`는 `"FFB9DE83"`(ARGB). alpha 2자 strip 후 매핑 (R `substr(_,3,8)`과 동일). 미인식 fill은 "Valid" 기본값.

---

## 7. 라이선스 gating 검증 결과 (완료)

**Biocrates 저작물(indicator 카탈로그)은 패키지에 없음 → GPL-3.0 copyleft blocker 없음.**

- `decisions.md #1`: redistribution license 문제는 **"Biocrates가 list-up한 indicator list에만"** 존재. 해당 레이어 전체 제거됨. 나머지(reader, metabolite_dict, 전처리/QC, 외부-DB pathway set)는 **공개 제약 없음**.
- 파일 확인: indicator 데이터 0개. 유일한 "indicator" 문자열은 `immunomet/pathways.py:83`의 생화학 주석.
- 번들 데이터: analyte roster TSV(사실 목록), immunomet + RaMP/LION 파생 GMT(license-clean 파생물), synthetic 예제.

### 동반 정리 항목 (blocker 아님)
1. **데이터 라이선스 별도 명시**: GPL-3.0은 코드 라이선스. roster/pathway 데이터는 코드=GPL-3.0과 별개로 명시 (roster는 사실 데이터, 외부 pathway는 아래 attribution 준수).
2. **외부-DB attribution 유지**: Reactome(CC-BY), WikiPathways(CC0), SMPDB(학술), LION, RaMP 인용 의무는 GPL이 대체하지 않음. DESIGN §8.5가 이미 추적.
3. **LICENSE/DESCRIPTION 교체**: 현 "All rights reserved" → GPL-3.0 전문. decisions.md 신규 항목 → 파일 변경(코드와 별도 commit, CLAUDE.md §4).

---

## 8. 미해결 질문

1. **InMoose API 매핑 검증**: `lmFit → eBayes → topTable`이 우리 입력(continuous log2 concentration, count 아님)에서 R limma와 동일 결과를 내는지 실측 필요. `voom`은 불필요(count 전용).
2. **GSEA backend 최종 선택**: `gseapy`(표준) vs `blitzgsea`(작은 set에서 안정적 gamma 근사). immunomet set이 3–10개로 작아(DESIGN §8.5.2) 이 선택이 결과 질에 영향.
3. **데이터 라이선스 문구 확정**: §7 동반 항목 1의 구체적 라이선스명 (roster·큐레이션 세트).
4. **패키지 배포 형태**: PyPI 배포 여부, 패키지명(`metabosetr`?), 최소 Python 버전.
5. **sysdata export 방식**: R 1회 export vs `data-raw`에서 Python 완전 재빌드 (후자가 R 의존 완전 제거에 부합).

---

## Appendix: 관련 문서
- `DESIGN.md` — 설계 single source of truth (R 기준, Python도 동일 설계 준수)
- `decisions.md` — §2의 결정 A/B/C가 여기 append될 예정
- `CLAUDE.md` — 작업 규칙
</content>
</invoke>

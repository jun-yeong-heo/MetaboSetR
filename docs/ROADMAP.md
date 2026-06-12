# MetaboSetR — Roadmap

> **Last updated**: 2026-06-12
> **Purpose**: 버전 로드맵, 현재 작업 상태, open questions. DESIGN.md가 "무엇을 어떻게
> 만들지"라면 이 문서는 "언제 어떤 순서로 만들지".

---

## Current status

**v0.1.0** — 내부 패키지 MetaboIndicatoR에서 indicator 레이어를 제거하여 파생
(→ decisions.md #1). WebIDQ 전처리 파이프라인 + metabolite-level 통계 + 큐레이션
pathway-set GSEA/ORA로 구성. `devtools::install()` clean, `devtools::test()` green.

공개 범위는 private로 시작하며, 논문 공개 시점에 함께 공개 예정 (→ decisions.md #5).

---

## v0.1.0 scope (shipped)

- 모든 reader (concentration, status text, status color fallback, threshold)
- Preprocessing pipeline (sample filter, metabolite filter 2-step)
- QC diagnostics 함수군 (PCA, PERMANOVA, CV/LOD plots, filter summary)
- Metabolite-level feature testing: `test_metabolites()` (limma / wilcoxon, BH/BY)
- 큐레이션 pathway-set enrichment (metabolite-level):
  - `gsea_pathway_sets()` — rank 기반 GSEA (fgsea)
  - `ora_pathway_sets()` — threshold 기반 hypergeometric ORA (universe = 측정 panel)
  - 도메인: `immunomet`, `reactome`, `wikipathways`, `smpdb`, `lion`, `source`, `health`
- Pathway set loader API (`list_pathway_sets`, `get_pathway_sets`, …)
- Metabolite id annotation: `annotate_metabolites()` (kegg/hmdb/pubchem/chebi/refmet/
  lipidmaps/lion) — KEGG-map 오버레이용
- `gsea_to_df()`, `log2_with_pseudocount()` helpers; Biocrates status predicates
- Synthetic data + unit tests (pathway set schema + snapshot 포함)

---

## Backlog

- Vignette 본문 작성 (현재 stub) — synthetic data 기반 worked example
- `validate_status()` 구현 — text vs color cross-validation (현재 stub; reader가 cell
  color를 동시에 retain하도록 import path 변경 필요)
- 추가 pathway-set domain — 본인 분석 범위 확장 시 재평가
- 공개 준비 — 논문 공개 시점에 LICENSE 전환 + 공개 repo 정리

---

## Open questions

### Biocrates 문의 대기
1. **`< threshold` status의 정확한 의미** — Valid로 취급하는 데는 합의됐으나, 어떤
   경우에 `< threshold`가 찍히는지는 미확정.
2. **Lipid metabolite LLOQ/ULOQ NA** — WebIDQ 자체 이슈. Pipeline 영향 없음 (status
   annotation 기반).

### 해소된 항목
- ~~Text annotation 파일 제공 여부~~ → 수령 완료. Primary method로 사용.

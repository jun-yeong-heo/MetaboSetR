# MetaboSetR — Key Design Decisions Log

> **Scope**: 설계 과정에서 내린 주요 의사결정의 이력. Append-only.
> **사용 원칙**:
> - 새 결정은 아래에 항목 추가.
> - 기존 결정을 뒤집어야 하면, 기존 항목은 **수정하지 않고** 새 항목으로 추가하면서
>   기존 항목에 `(superseded by #N)` 표기.
> - "Related section"은 DESIGN.md의 해당 절과의 역참조.
>
> 이 로그는 MetaboSetR 시작 시점(#1)부터 새로 시작한다. MetaboSetR은 내부 패키지
> MetaboIndicatoR에서 indicator 레이어를 제거하여 파생되었으며, indicator 시절의
> 의사결정 history는 본 로그에 이월하지 않는다 (#1 참조).

---

### #1 — MetaboIndicatoR에서 indicator 레이어를 제거하여 MetaboSetR로 분리
**선택**: Biocrates indicator 계산 레이어 전체(`compute_indicators`, formula parser/AST,
`IndicatorProfile`, `indicators_meta`, indicator-level pathway set)를 제거하고, WebIDQ
전처리 파이프라인 + metabolite-level 통계 + 큐레이션 pathway-set enrichment만 남긴 새
패키지를 fresh git history로 분리.

**근거**: redistribution license 문제가 **Biocrates가 list-up한 indicator list에만**
존재한다. indicator list를 제외하면 나머지(reader, metabolite_dict, 전처리/QC, 외부-DB
pathway set)는 공개 제약이 없다. 공개는 논문과 함께 진행 예정이라, 그 전까지 private로
시작한다. indicator-level pathway set은 멤버가 Biocrates indicator short_name이라
license 경계 안에 들어가므로 함께 제거했고, 그 결과 indicator-only였던 cancer/aging/
neurodegen 도메인은 통째로 사라졌다(실제 분석에서 효익이 낮았음 — 사용자 확인). immunomet은
metabolite-level만 유지. **Related**: §1, ROADMAP — v0.1

### #2 — 통계·GSEA 레이어를 metabolite-level로 재배선
**선택**: `test_indicators(IndicatorProfile)` → `test_metabolites(PreprocessedData)`.
`gsea_pathway_sets(prep, group, domain)`가 `PreprocessedData`의 metabolite 농도에서
log2FC rank를 계산. indicator operand-set GSEA(구 Mode B1)와 `level` 인자는 제거 —
잔존 pathway set이 전부 metabolite-level이라 `level` 차원이 무의미.

**근거**: indicator 행렬이 사라졌으므로 통계의 입력은 metabolite 농도 행렬
(`PreprocessedData`)이 자연스러운 단위. `level` 인자는 상수가 되어 indicator/metabolite
이원성의 흔적만 남기므로 제거. **Related**: §6.1, §6.2

### #3 — ORA universe = 측정 panel (Option A)
**선택**: `ora_pathway_sets(sig, universe, domain)`의 background는 사용자가 명시하는 측정
panel 전체. 각 pathway는 universe와 교집합 후 검정.

**근거**: 고정 패널(Biocrates Q1000)에서 측정하지 않은 metabolite는 애초에 유의해질 수
없으므로 분모에서 제외하는 것이 정직하다. 내장 세트가 이미 panel-driven 큐레이션(#5)이라
정합. Direction은 검정에 미전달(해석용 metadata). **Related**: §6.6

### #4 — 큐레이션 pathway set: immunomet + 외부-DB(RaMP/LION)
**선택**: 유지 도메인 — `immunomet`(O'Neill 2016, metabolite-level), 그리고 RaMP-DB/LION
파생 metabolite-level 도메인 `reactome`/`wikipathways`/`smpdb`/`lion`/`source`/`health`.
각 외부 도메인은 panel-driven 선택(Biocrates 매칭 metabolite와 교집합, 멤버 ≥3, 동일
멤버십 dedupe), direction `-`(hypothesis-free).

**근거**: 외부 DB raw(RaMP SQLite, LION CSV)는 gitignore하고 license-clean 파생물만
commit. KEGG bulk/map은 라이선스로 미embedding하되, RaMP cross-ref 경유 `kegg_id`는
license-clean이라 `metabolite_dict`에 annotation으로 제공(pathview는 사용자가 직접 호출).
**Related**: §5.2, §6.2, §8.5

### #5 — 배포 범위: private 시작, 논문과 함께 공개
**선택**: GitHub private repo로 시작. 외부 배포는 논문 공개 시점에 함께.
`devtools::install()` clean이 의무 바, `R CMD check`는 의무 아님.

**근거**: indicator list license 경계는 #1로 제거했으나, 공개 타이밍은 논문에 종속. 그
전까지는 install-clean + testthat green을 품질 바로 둔다. **Related**: §1, CLAUDE.md §3

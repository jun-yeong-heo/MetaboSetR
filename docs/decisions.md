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

### #5 — 배포 범위: private 시작, 논문과 함께 공개 (partially superseded by #11)
**선택**: GitHub private repo로 시작. 외부 배포는 논문 공개 시점에 함께.
`devtools::install()` clean이 의무 바, `R CMD check`는 의무 아님.

**근거**: indicator list license 경계는 #1로 제거했으나, 공개 타이밍은 논문에 종속. 그
전까지는 install-clean + testthat green을 품질 바로 둔다. **Related**: §1, CLAUDE.md §3

### #6 — 순수 Python 포트로 재설계 (R 런타임 의존 제거) (partially superseded by #10)
**선택**: R 패키지를 순수 Python으로 재이식하여 `python/metabosetr/`에 둔다. 전처리
→ metabolite 통계(A) → GSEA(B) → ORA(C) → pathway 로더/annotate/QC까지 전 계층 이식.
rpy2 미사용, R 런타임 의존 완전 제거. R 판본은 repo에 유지(참조/대조용)하되, 이후
개발의 기준은 Python 포트로 이동. 파일-단위 매핑·이식 순서·수치 재현성 등급은
`MIGRATION.md`에 정리.

**근거**: 대부분의 계층이 Python 생태계에 1:1 대응이 있고(readxl→openpyxl,
tidyxl→openpyxl, vegan→scikit-bio, phyper→scipy 등), pathway set build script가 이미
Python이라 절반이 기이식 자산. 통계 계층의 재현성 리스크는 #7·#8로 해소. 재현성은
"R과 bit-identical"이 아니라 "Python 내 결정론"으로 재정의(전처리·ORA·Wilcoxon·
p.adjust는 R base와 ~1e-14 일치 검증; GSEA만 설계상 상이). **Related**: §7, MIGRATION.md

### #7 — Approach A 통계 엔진: R limma → InMoose (Python 포트)
**선택**: Approach A(`test_metabolites`)의 default 엔진을 R `limma`에서 Python
`inmoose.limma`로 교체. **방법(empirical Bayes moderated-t)과 default 지위는 불변**.
Wilcoxon(`scipy.stats.mannwhitneyu`)은 support 옵션으로 유지. DESIGN §6.1의 방법론은
그대로이고 backend만 바뀐다.

**근거**: InMoose는 R limma의 소스 충실 이식(peer-reviewed drop-in, near-identical
결과)이라, moderated-t를 default로 유지하면서도 R 의존을 끊고 feature 통계의 수치
재현성을 보전할 수 있다. 완전 베이지안(BEST/PyMC)은 MCMC·무거운 의존·재현성 seed
부담으로 과하다고 판단해 배제. **Related**: §6.1, §7.4

### #8 — GSEA backend: fgsea → gseapy (§6.2.1 논거 갱신)
**선택**: Approach B(`gsea_pathway_sets`)의 backend를 R `fgsea::fgseaMultilevel`에서
Python `gseapy` prerank로 교체. R fgsea와의 **수치 동일성은 포기**(수용). 재현성은
**seed 고정으로 결정론 보장**(API에 `seed` 인자 노출). §6.2.1의 "fgsea 사용 이유
(method paper)" 논거는 Python 포트에서 무효화되며 본 결정으로 대체.

**근거**: 순수 Python 재설계에서 R fgsea와 bit-identical을 요구하지 않기로 한 이상
(#6), 남는 요건은 Python 내 결정론뿐이고 이는 seed 고정으로 충족된다. `limma-voom`
불필요 논거(continuous concentration)는 유지. **Related**: §6.2, §6.2.1

### #9 — 라이선스: proprietary → GPL-3.0 (copyleft)
**선택**: repo 라이선스를 기존 "All rights reserved / proprietary"에서 **GPL-3.0**으로
전환. `LICENSE`를 GPL-3.0 전문으로 교체, R `DESCRIPTION`의 `License`를 `GPL-3`으로
갱신. 코드=GPL-3.0. 번들 **데이터**(reference roster·pathway 세트)의 라이선스는 코드와
별개로 추후 명시하며, 외부-DB 파생 세트의 attribution 의무(Reactome CC-BY,
WikiPathways CC0, SMPDB, LION, RaMP)는 그대로 유지(§8.5).

**근거**: copyleft 전환의 유일한 gating이었던 Biocrates indicator 카탈로그가 #1로
이미 제거되어 패키지에 없음(파일 검증 완료). GPL-3.0 채택은 핵심 의존 InMoose(GPL-3.0)
와 정합하고, 기존 R 의존(limma/fgsea, GPL)과도 posture가 동일하다. #5의 "private 시작,
논문과 함께 공개" 타이밍은 불변 — 본 결정은 공개 시 적용될 라이선스를 확정할 뿐 repo
가시성을 바꾸지 않는다. **Related**: §1.2, LICENSE, DESCRIPTION, #1, #5

### #10 — 구현 기준을 R 판본으로 환원 (#6의 "기준 이동" 무효화)
**선택**: 개발·배포의 기준 구현을 **R 패키지**로 환원한다. `python/metabosetr/`의 순수
Python 포트는 제거하지 않고 **참조·대조 구현**으로 유지(R 런타임 없는 환경에서의 실행,
수치 대조용). #6 중 "이후 개발의 기준은 Python 포트로 이동" 부분만 무효화되며, 이식
자체와 #7(InMoose)·#8(gseapy)은 Python 포트 내부의 backend 결정으로 계속 유효하다.
README는 R 사용법을 우선 제시하고 Python 포트는 참조 섹션으로 안내한다.

**근거**: 패키지 정체성과 배포 경로(DESCRIPTION/NAMESPACE/man/vignettes)가 R에 있고,
DESIGN §7.1 디렉토리 구조와 §7.3 API도 R 기준으로 서술되어 있다. 두 판본을 동시에
"기준"으로 두면 API·문서 drift 비용이 크므로, 기준은 R 하나로 두고 Python 포트는
동등 기능을 가진 대조 자산으로 유지한다. **Related**: §7.1, §7.3, §7.4, #6, #7, #8

### #11 — 공개 시점 변경: 논문 이전 public 전환
**선택**: GitHub repo를 **public으로 전환**. #5의 "논문 공개 시점에 함께 공개" **타이밍만**
무효화하고, #5의 품질 바(`devtools::install()` clean 의무, `R CMD check` 비의무)는 유지.
라이선스는 #9의 GPL-3.0 그대로.

**근거**: 공개의 유일한 gating이던 Biocrates indicator 카탈로그는 #1로 제거되어 패키지에
없고(MIGRATION.md §7 파일 검증 완료), #9로 재배포 조건(GPL-3.0)도 확정됐다. 번들 데이터의
외부-DB attribution 의무(Reactome CC-BY 4.0, WikiPathways CC0, SMPDB, LION, RaMP-DB)는
공개 여부와 무관하게 §8.5대로 유지된다. **Related**: §1.2, §8.5, #5, #9

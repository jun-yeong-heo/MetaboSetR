# CLAUDE.md — MetaboSetR 작업 규칙

이 파일은 Claude Code 세션 시작 시 자동으로 읽힙니다. 모든 작업은 아래 규칙을 따라야 합니다.

---

## 0. 기본 지침 (Default behavioral guidelines)

LLM 코딩에서 흔히 발생하는 실수를 줄이기 위한 기본 원칙. 프로젝트별 규칙(§1 이하)과 충돌하면 프로젝트 규칙이 우선한다. Trivial 작업은 판단에 맡기되, 의심이 들면 caution을 택한다.

### 0.1 Think Before Coding — 가정하지 말고, 혼란을 숨기지 말고, 트레이드오프를 드러내라
구현 전에:
- 가정은 명시적으로 말한다. 불확실하면 묻는다.
- 해석이 여러 개면 모두 제시한다. 침묵 속에 하나를 고르지 않는다.
- 더 단순한 접근이 있으면 말한다. 필요하면 사용자 요청에 push back 한다.
- 불명확하면 멈춘다. 무엇이 헷갈리는지 이름 붙이고, 묻는다.

### 0.2 Simplicity First — 문제를 푸는 최소한의 코드. 추측성 코드 금지
- 요청되지 않은 기능 추가 금지.
- 단발성 코드를 위한 추상화 금지.
- 요청되지 않은 "유연성" / "설정 가능성" 금지.
- 발생할 수 없는 시나리오에 대한 error handling 금지.
- 200줄 짠 코드가 50줄로 가능하면 다시 쓴다.

자기 점검: "시니어 엔지니어가 보면 과하다고 할까?" 그렇다면 단순화.

### 0.3 Surgical Changes — 꼭 필요한 부분만 건드린다. 자신이 만든 흔적만 정리한다
기존 코드를 수정할 때:
- 인접한 코드, 주석, 포맷을 "개선"하지 않는다.
- 망가지지 않은 것을 refactor 하지 않는다.
- 본인 취향과 다르더라도 기존 스타일에 맞춘다.
- 무관한 dead code를 발견하면 **언급은 하되 삭제하지 않는다.**

변경으로 인한 orphan은 정리한다:
- 본인 변경 때문에 사용되지 않게 된 import / 변수 / 함수는 제거한다.
- 사전에 존재하던 dead code는 요청 없이 제거하지 않는다.

테스트: 변경된 모든 줄이 사용자 요청과 직접 연결되어야 한다.

### 0.4 Goal-Driven Execution — 성공 기준을 정의하고, 검증될 때까지 루프를 돈다
작업을 검증 가능한 목표로 변환:
- "Validation 추가" → "잘못된 입력에 대한 테스트를 먼저 쓰고, 그것을 통과시켜라"
- "버그 수정" → "버그를 재현하는 테스트를 먼저 쓰고, 그것을 통과시켜라"
- "X refactor" → "변경 전과 후 모두 테스트가 통과하는지 확인하라"

Multi-step 작업이면 짧은 계획을 명시:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

강한 성공 기준은 독립적인 루프를 가능하게 한다. 약한 기준("일단 동작하게")은 지속적인 clarification을 요구한다.

### 0.5 효과 측정
다음이 줄어들면 본 지침이 작동 중이라는 신호:
- diff 내 불필요한 변경
- 과복잡 구현으로 인한 재작성
- 실수 이후의 사후 clarification (질문은 실수 *전에* 와야 한다)

---

## 1. 문서 체계와 권위

프로젝트 문서는 세 개로 분리되어 있습니다. 각 문서의 역할과 권위는 다릅니다.

- **`docs/DESIGN.md`** — 설계의 single source of truth. "무엇을 왜 그렇게 설계했는가". 모든 구현은 이 문서를 근거로 한다.
- **`docs/ROADMAP.md`** — 버전 로드맵, 구현 순서, open questions. 작업 진행에 따라 갱신되는 문서. **작업 우선순위는 여기의 Implementation order를 따른다.**
- **`docs/decisions.md`** — 의사결정 로그. Append-only. DESIGN.md의 각 섹션이 "왜 그 선택을 했는지"의 근거가 여기 쌓인다.

### 작업 중 규칙
- DESIGN.md에 명시되지 않은 결정이 필요하면 **코드를 작성하기 전에 사용자에게 먼저 물어본다.**
- DESIGN.md와 다른 판단이 더 합리적으로 보일 때도 임의로 따르지 않는다. 제안은 하되, 변경은 사용자 승인 후 다음 순서로 진행:
  1. `decisions.md`에 새 결정 추가 (필요 시 기존 결정에 `(superseded by #N)` 표기)
  2. `DESIGN.md` 본문 수정
  3. 코드 수정
- 기존 `decisions.md` 항목은 수정하지 않는다. Append-only.

## 2. 코딩 컨벤션

### 2.1 R 코드
- **Tidyverse 사용 금지** (dplyr, tidyr, purrr, magrittr, readr 등). `dplyr::summarise()` 같은 함수도 base R로 대체한다.
- Base R + Bioconductor 생태계만 사용. 허용 의존성은 DESIGN.md §7.4 참조.
- 모든 export 함수는 roxygen2 문서화. `@param`, `@return`, `@export` 필수. `@examples`는 권장.
- S4 class는 `methods::setClass()`로 정의하고, 모든 slot에 타입 명시.
- 함수/객체 이름은 snake_case. S4 class는 PascalCase.
- 파일 인코딩 UTF-8, LF line endings.

### 2.2 Pathway set 컨벤션
모든 pathway set은 metabolite-level이다 (멤버 = metabolite `short_name`).

- **Set naming**: 대문자 + 언더스코어, domain prefix 필수.
  - 예: `IMMUNOMET_ONEIL_GLYCOLYSIS`, `REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE`
  - Reference-based set은 reference 축약 prefix 포함: `IMMUNOMET_ONEIL_*` (O'Neill 2016 출처 명시)
- **Master TSV column 컨벤션** (metabolite-level):
  - 필수: `shortname`, `fullname`, `analyte_class`, `set_name`, `match_basis`
  - 외부-DB 파생 도메인(reactome 등)은 build.py가 emit하는 컬럼을 따른다.
- **Direction 값** (보유 시): `up` | `down` | `altered` | `-` 넷 중 하나. 다른 값 금지.
  - `up`: 해당 domain state에서 metabolite 값 증가
  - `down`: 감소
  - `altered`: 방향성 불명확 또는 context-dependent
  - `-`: direction 정보 없음 (비어 있다는 뜻이 아니라 "알려지지 않음"을 명시)
  - 외부-DB 도메인은 hypothesis-free이므로 direction `-`. immunomet metabolite master는
    direction 컬럼 자체가 없어 NA로 로드됨 (둘 다 허용).

## 3. 테스트 정책

본 정책은 **private 시작 단계의 격하판**이다 (→ decisions.md #5). 외부 공개가 가까워지면 strict 모드로 환원.

- **핵심 모듈은 testthat 단위 테스트 필수**: metabolite statistics, status helper, pathway-set loader/ORA 등 입출력이 명확하고 정답이 검증 가능한 모듈. Reader / QC diagnostic / plot / I/O 등은 **best-effort** (smoke test 정도면 충분).
- 로직이 명확한 모듈(통계·ORA hypergeometric 등)은 **TDD**로 진행:
  - Hand-computed ground truth 5-10개를 먼저 testthat 케이스로 작성
  - 그다음 구현
- Pathway set 관련 테스트 최소 항목:
  - Loader 함수가 예상 set 개수/멤버 개수 반환하는지
  - GMT 파일이 `fgsea::gmtPathways()`로 정상 파싱되는지
  - `sysdata` 내 list 구조 schema 검증
  - 알려진 특정 멤버가 예상 set에 실제로 들어있는지 (snapshot test)
- 검증 루프: `devtools::load_all()` → `devtools::test()` → `devtools::install()`
- **`devtools::install()`이 에러 없이 통과해야 한다**. `R CMD check` / `devtools::check()` 통과는 의무 아님.
- 커버리지 목표는 강제하지 않음. 단, metabolite statistics·pathway-set 같은 핵심 모듈은 의도적으로 높게 유지.

## 4. Git 워크플로우

- **작은 단위 커밋**. 각 작업 step이 끝날 때마다 commit. 한 commit이 여러 무관한 변경을 묶지 않는다.
- 무관한 변경은 **묶지 않을 뿐 아니라 애초에 만들지 않는다** (§0.3). 사용자 요청과 직접 연결되지 않는 drive-by refactor / 인접 코드 개선 / 포맷 정리 / 미사용처럼 보이는 사전 dead code 제거는 발견 시 별도 commit 또는 별도 사용자 승인을 받는다.
- Commit 메시지는 명령형 영문, 첫 줄 50자 이내, 본문은 왜(why)를 설명.
- DESIGN.md / ROADMAP.md / decisions.md 변경은 코드 변경과 별도 commit.
- Pathway set source 변경과 regenerated artifact(GMT, sysdata)는 가급적 같은 commit 안에 묶는다 (재현성 확인 용이).

## 5. 데이터 정책

### 5.1 실제 sample 데이터
- **실제 sample 데이터는 절대 commit 금지.**
  - `data-raw/real_samples/`는 `.gitignore`에 등록.
- Synthetic example data만 `inst/extdata/` 또는 `data/`에 둔다.

### 5.2 Reference 데이터 (vendor-provided, 불변)
- `data-raw/reference/`에만 포함:
  - `biocrates_Quant1000_metabolites.tsv` — 1234 metabolite reference (Biocrates 공식 PDF에서 파싱)
- Vendor가 새 버전을 내기 전까지 수정 금지. Schema 변경 시 DESIGN.md §8 업데이트 먼저.
- 참고: Biocrates가 list-up한 indicator catalogue는 license 문제로 본 패키지에서 제외됨
  (→ decisions.md #1).

### 5.3 Pathway set 데이터 (자체 curation, editable)
- `data-raw/pathway_sets/<domain>/`에 source 유지. Domain별 디렉토리.
- 각 domain 디렉토리에 포함되는 것:
  - Master TSV (long-format, source of truth)
  - Build script (Python — 현재 컨벤션상 data-raw 내 혼용 허용)
  - Documentation (`pathway_definitions.md`, 필요 시 `keyword_rules.md`)
  - Assignment log (audit용 TSV)
- Pathway set의 installed artifacts:
  - `inst/extdata/pathway_sets/*.gmt` — build 과정에서 자동 생성, 손으로 수정 금지
  - `R/sysdata.rda` 내 `pathway_sets`, `pathway_sets_meta` 객체 — 동일 source에서 R로 재빌드

### 5.4 Build 도구 예외 조항
- `data-raw/` 내 build script는 **Python 사용 허용** (R 본체의 tidyverse 금지 원칙은 유지).
- 이유: pathway set source의 대부분이 Python으로 큐레이션되어 있고, 이를 R로 재작성할 실익이 낮음.
- Installed package의 R 코드는 여전히 R-only.

## 6. sysdata 빌드 절차

### 6.1 생성되는 객체
`R/sysdata.rda`에 들어가는 객체:

| 객체 | 타입 | 소스 |
|---|---|---|
| `metabolite_dict` | data.frame | `data-raw/reference/biocrates_Quant1000_metabolites.tsv` + RaMP/LION id |
| `status_color_map` | data.frame | DESIGN.md §2.4 hex 매핑 하드코딩 |
| `pathway_sets` | nested named list | `data-raw/pathway_sets/<domain>/*.tsv` (domain → set → members) |
| `pathway_sets_meta` | data.frame | 동일 소스, metadata 통합 |

### 6.2 빌드 순서
1. `data-raw/pathway_sets/<domain>/build.py` 실행 → 각 domain의 TSV/GMT 재생성
2. 생성된 GMT를 `inst/extdata/pathway_sets/`로 복사 또는 심볼릭 생성
3. `data-raw/make_sysdata.R` 실행:
   - Reference TSV 로드
   - Pathway set TSV들 로드하여 `pathway_sets` 중첩 list 구성
   - `usethis::use_data(..., internal = TRUE)` 로 `R/sysdata.rda` 갱신

### 6.3 Placeholder 정책
- `pathway_sets`가 특정 도메인만 존재할 수 있음. 미존재 도메인을 요청하면 명시적 에러.

## 7. 작업 진행 방식

- **한 번에 한 step만 진행**. 사용자 리뷰 없이 다음 step으로 자동 진행 금지.
- **각 step 시작 시 성공 기준을 한 줄로 명시** (§0.4). 가능한 한 검증 가능한 형태로. 예: "parser가 hand-computed ground truth 5개를 모두 통과", "`devtools::install()` 에러 없이 통과". 기준이 모호하거나 정의하기 어려우면 그 자체를 사용자에게 먼저 묻는다 (§0.1).
- 각 step 끝에 다음을 보고:
  - 만든/수정한 파일 목록
  - 추가한 의존성 (있으면)
  - 시작 시 명시한 성공 기준에 대한 검증 결과 (`devtools::install()` + `devtools::test()` 결과 포함)
  - 다음 step 제안
- 큰 변경(파일 5개 이상 동시 수정, 새 의존성 추가, S4 class 시그니처 변경)은 착수 전 계획을 먼저 보고하고 승인 받기.

### 7.1 Pathway set 변경 절차
Pathway set 추가/수정은 다음 순서로 수행, 각 단계 별도 commit:

1. `data-raw/pathway_sets/<domain>/`에서 source (TSV + docs + build.py) 수정
2. `python3 build.py` 실행하여 TSV/GMT 재생성, 출력 확인
3. 생성된 GMT를 `inst/extdata/pathway_sets/`에 동기화
4. `data-raw/make_sysdata.R` 재실행하여 `R/sysdata.rda` 갱신
5. Loader 함수와 테스트 업데이트 (필요 시)
6. `devtools::install()` + `devtools::test()` 통과 확인

## 8. 보고/응답 언어

- 사용자와의 대화는 한국어. 코드 주석/roxygen/commit message는 영어.

## 9. 환경 가정

- R ≥ 4.3
- Bioconductor ≥ 3.18 (fgsea, limma 의존)
- 개발 도구: `devtools`, `roxygen2`, `testthat (>= 3.0)`, `usethis`
- Python ≥ 3.10 (pathway set build scripts용, stdlib만 사용 권장)

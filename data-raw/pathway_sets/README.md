# Pathway Sets — Source of Truth

이 디렉토리에는 MetaboSetR 패키지에 내장될 pathway set의 **source (수정 가능한 원본)** 가 도메인별로 들어간다. 실제로 패키지에 배포되는 `inst/extdata/pathway_sets/*.gmt`와 `R/sysdata.rda`의 `pathway_sets` / `pathway_sets_meta` 객체는 이 디렉토리의 파일로부터 build된다.

모든 pathway set은 **metabolite-level** 이다 (멤버 = metabolite `short_name`).

상세 설계는 `docs/DESIGN.md` §8.5 참조. 작업 규칙은 `CLAUDE.md` §5.3, §7.1 참조.

---

## 디렉토리 구조

```
data-raw/pathway_sets/
├── README.md                       ← 이 파일: 공통 규칙 + 추가 절차
├── immunomet/                       ← in-house 큐레이션 (O'Neill 2016, 5 set)
│   ├── README.md
│   ├── build.py                    ← metabolite-level master/GMT 재생성
│   ├── pathways.py                 ← canonical metabolite 정의 + 예외 처리
│   ├── immunomet_metabolites_master.tsv
│   ├── immunomet_metabolites.gmt
│   ├── assignment_log_metabolites.tsv
│   └── pathway_definitions.md
├── ramp/                           ← 외부 DB 파생: reactome + wikipathways + smpdb + source + health
│   ├── build.py                    ← RaMP-DB 매칭 + 도메인 set 생성 (단일 script)
│   ├── ramp_metabolite_ids.tsv     ← metabolite_dict 용 kegg/hmdb/pubchem/chebi/refmet/lipidmaps id
│   ├── {reactome,wikipathways,smpdb,source,health}_master.tsv
│   ├── {reactome,wikipathways,smpdb,source,health}_metabolites.gmt
│   ├── assignment_log.tsv
│   └── raw/                        ← RaMP SQLite (gitignored, ~2GB)
└── lion/                           ← 외부 DB 파생: lipid functional ontology
    ├── build.py                    ← LION 매칭 + term→lipid set 생성
    ├── lion_lipid_ids.tsv          ← metabolite_dict 용 lion_id annotation
    ├── lion_master.tsv
    ├── lion_metabolites.gmt
    ├── assignment_log.tsv
    └── raw/                        ← LION association CSV (gitignored)
```

`ramp`/`lion` 은 외부 DB에서 build-time에 매칭·생성되는 metabolite-level 도메인으로,
in-house 큐레이션(immunomet)과 절차가 다르다 (→ decisions.md #4).

---

## 공통 규칙

### Set naming

- 대문자 + 언더스코어만 사용
- Domain prefix 필수
- Reference 기반 domain은 reference 약어 포함

| 패턴 | 예시 |
|---|---|
| `<DOMAIN>_<SET>` | `REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE`, `LION_MEMBRANE_FLUIDITY` |
| `<DOMAIN>_<REFERENCE>_<SET>` | `IMMUNOMET_ONEIL_GLYCOLYSIS` |

### Master TSV schema (metabolite-level)

| 컬럼 | 필수 | 설명 |
|---|---|---|
| `shortname` | 필수 | metabolite shortname (metabolite reference TSV와 일치) |
| `fullname` | 필수 | metabolite fullname |
| `analyte_class` | 필수 | metabolite의 analyte_class |
| `set_name` | 필수 | 해당 set name |
| `match_basis` | 필수 | `individual` 또는 `class:<analyte_class>` |

외부-DB 파생 도메인(reactome 등)은 위 base에 DB-specific evidence 컬럼을 추가한다 (build.py가 emit).

### Direction 값 (enum)

보유 시 `up` | `down` | `altered` | `-` 넷 중 하나 (CLAUDE.md §2.2). 외부-DB 도메인은
hypothesis-free라 `-`. immunomet metabolite master는 direction 컬럼이 없어 NA로 로드된다.

### Overlap 정책

- **생물학적으로 의미있는 overlap은 허용** — 한 metabolite가 여러 set에 들어갈 수 있음
- **데이터 구조 artifact로 인한 overlap은 금지** — 화학 class 분류 때문에 우연히 걸린 경우 명시적 제외
- 명시적 제외는 build script 내부에 documented exclusion list로 관리

---

## Build 절차

새 domain을 추가하거나 기존 domain을 수정한 후:

1. **해당 domain 디렉토리에서 build script 실행**
   ```bash
   cd data-raw/pathway_sets/<domain>
   python3 build.py
   ```
   → master TSV + GMT + assignment log 재생성

2. **생성된 GMT를 `inst/extdata/pathway_sets/`로 복사**
   ```bash
   cp data-raw/pathway_sets/<domain>/*.gmt inst/extdata/pathway_sets/
   ```

3. **sysdata 재빌드**
   ```r
   source("data-raw/make_sysdata.R")
   ```

4. **테스트 통과 확인**
   ```r
   devtools::test()
   devtools::install()
   ```

5. **각 단계 별도 commit** (CLAUDE.md §7.1)

---

## 세 가지 저장 형식의 관계

| 위치 | 형식 | 용도 | 유지 방식 |
|---|---|---|---|
| `data-raw/pathway_sets/<domain>/*.tsv` | long-format TSV | Source of truth | 사람이 수정 / build.py 생성 |
| `data-raw/pathway_sets/<domain>/*.gmt` | GMT | 외부 GSEA 툴 호환용 로컬 백업 | build.py가 자동 생성 |
| `inst/extdata/pathway_sets/*.gmt` | GMT | 설치된 패키지에서 `system.file()`로 접근 | 위의 GMT를 복사 |
| `R/sysdata.rda` 내 `pathway_sets`, `pathway_sets_meta` | R list / data.frame | 패키지 런타임 빠른 접근 | `make_sysdata.R`가 master TSV에서 재생성 |

---

## 참고

- 상위 설계: `docs/DESIGN.md` §8.5 Pathway Sets
- 작업 규칙: `CLAUDE.md` §5.3, §7.1
- GSEA/ORA에서의 활용: `docs/DESIGN.md` §6.2, §6.6

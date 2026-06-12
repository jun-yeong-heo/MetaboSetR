#!/usr/bin/env python3
"""Build RaMP-derived artifacts for MetaboSetR.

Step 1 (this stage): match Biocrates metabolites to RaMP and emit a
per-metabolite ID annotation table (kegg_id/hmdb_id/ramp_id) plus an
assignment log. Set-building (reactome/wikipathways/smpdb domains) is
added in a later stage and reuses the matching functions here.

Source: data-raw/pathway_sets/ramp/raw/RaMP_SQLite_v3.0.7.sqlite (gitignored).
Run from project root:  python3 data-raw/pathway_sets/ramp/build.py
"""
import csv
import os
import re
import sqlite3
from collections import defaultdict

ROOT = os.getcwd()
HERE = os.path.join(ROOT, "data-raw/pathway_sets/ramp")
DB = os.path.join(HERE, "raw/RaMP_SQLite_v3.0.7.sqlite")
BIOC = os.path.join(ROOT, "data-raw/reference/biocrates_Quant1000_metabolites.tsv")

ID_OUT = os.path.join(HERE, "ramp_metabolite_ids.tsv")
LOG_OUT = os.path.join(HERE, "assignment_log.tsv")

MIN_MEMBERS = 3
# RaMP pathway.type -> our domain. kegg-type (license) and pfocr (figure
# OCR, low quality) are intentionally excluded (decisions.md #65).
PATHWAY_TYPE_DOMAIN = {"reactome": "reactome", "wiki": "wikipathways",
                       "hmdb": "smpdb"}

# priority order for picking a canonical rampId when a name is ambiguous;
# higher pathwayCount wins first, this only breaks remaining ties.
HMDB_STATUS_RANK = {"quantified": 5, "detected": 4, "expected": 3,
                    "predicted": 2, "no_HMDB_status": 1}

LIPID_CLASSES = {
    "Acylcarnitines", "Ceramides", "Cholesteryl esters", "Diglycerides",
    "Dihexosylceramides", "Dihydroceramides", "Fatty acids",
    "Hexosylceramides", "Lysophosphatidic acids", "Lysophosphatidylcholines",
    "Lysophosphatidylethanolamines", "Lysophosphatidylglycerols",
    "Lysophosphatidylinositols", "Lysophosphatidylserines", "Monoglycerides",
    "Phosphatidic acids", "Phosphatidylcholines", "Phosphatidylethanolamines",
    "Phosphatidylglycerols", "Phosphatidylinositols", "Phosphatidylserines",
    "Sphinganine and sphingosine phosphates", "Sphinganines and sphingosines",
    "Sphingomyelins", "Triglycerides", "Trihexosylceramides",
}


def norm(s):
    if not s:
        return ""
    s = s.strip().lower()
    s = re.sub(r"^(dl|l|d)-", "", s)
    s = s.replace("-l-", "-").replace("-d-", "-")
    for g, r in (("alpha", "a"), ("beta", "b"), ("gamma", "g"), ("omega", "o")):
        s = s.replace(g, r)
    s = re.sub(r"[^a-z0-9]", "", s)
    return s


def no_paren(s):
    return re.sub(r"\s*\([^)]*\)", "", s).strip()


def stem(s):
    """Normalize + collapse acid/conjugate-base naming (aspartic acid <->
    aspartate) so stereo/charge forms of one compound group together. Used
    only for the canonical-name filter, not for synonym lookup."""
    n = norm(s)
    for suf in ("oicacid", "icacid", "oate", "ate"):
        if n.endswith(suf) and len(n) - len(suf) >= 3:
            return n[:-len(suf)]
    return n


EXTRA_ID_TYPES = ("pubchem", "chebi", "refmet", "lipidmaps")


def _pick_hmdb(ids):
    """Prefer canonical 11-char HMDB0####### form."""
    canon = sorted(x for x in ids if re.fullmatch(r"HMDB\d{7}", x))
    return canon[0] if canon else (sorted(ids)[0] if ids else "")


def _pick_id(idtype, ids):
    """Representative id for a cross-ref column (one per metabolite)."""
    if not ids:
        return ""
    if idtype in ("pubchem", "chebi"):                  # numeric -> lowest
        best = min(ids, key=lambda x: (int(x) if x.isdigit() else 1 << 62, x))
        return f"CHEBI:{best}" if idtype == "chebi" else best
    return sorted(ids)[0]                               # refmet/lipidmaps: lexical


def load_ramp(con):
    """Return per-rampId matching record `rec` (kegg/hmdb only — drives
    candidate resolution, must stay identical to keep id assignments
    stable), an `extra` id table (pubchem/chebi/refmet/lipidmaps, for
    annotation only), and synonym indices. Compounds only."""
    rec = {}    # rampId -> {kegg:set, hmdb:set, pwc:int, rank:int, cn:str}
    extra = {}  # rampId -> {pubchem:set, chebi:set, refmet:set, lipidmaps:set}
    c = con.cursor()
    for rid, idt, sid, pwc, status in c.execute(
        "SELECT rampId, IDtype, sourceId, pathwayCount, priorityHMDBStatus "
        "FROM source WHERE geneOrCompound='compound' AND IDtype IN ('kegg','hmdb')"):
        r = rec.get(rid)
        if r is None:
            r = rec[rid] = {"kegg": set(), "hmdb": set(), "pwc": 0, "rank": 0}
        r["pwc"] = max(r["pwc"], pwc or 0)
        r["rank"] = max(r["rank"], HMDB_STATUS_RANK.get(status, 0))
        val = sid.split(":", 1)[1] if ":" in sid else sid
        r[idt].add(val)

    # extra cross-ref id types — annotation only, NOT used for matching.
    for rid, idt, sid in c.execute(
        "SELECT rampId, IDtype, sourceId FROM source WHERE "
        "geneOrCompound='compound' AND IDtype IN "
        "('pubchem','chebi','refmet','LIPIDMAPS')"):
        e = extra.setdefault(rid, {k: set() for k in EXTRA_ID_TYPES})
        e[idt.lower()].add(sid.split(":", 1)[1] if ":" in sid else sid)

    # canonical name per rampId — used to filter out synonym-noise candidates
    # (e.g. a glycan/protein that carries "alanine" as a stray synonym).
    for rid, cn in c.execute("SELECT rampId, common_name FROM analyte"):
        if rid in rec:
            rec[rid]["cn"] = cn or ""

    syn_exact = defaultdict(set)
    syn_norm = defaultdict(set)
    for syn, rid in c.execute(
        "SELECT Synonym, rampId FROM analytesynonym WHERE geneOrCompound='compound'"):
        if not syn:
            continue
        syn_exact[syn.strip().lower()].add(rid)
        syn_norm[norm(syn)].add(rid)
    return rec, extra, syn_exact, syn_norm


def candidates(sn, fn, syn_exact, syn_norm):
    """First non-empty tier of candidate rampIds; returns (rampIds, basis)."""
    tiers = [("exact:fullname", syn_exact, fn.strip().lower()),
             ("norm:fullname", syn_norm, norm(fn)),
             ("norm:fullname_noparen", syn_norm, norm(no_paren(fn)))]
    if len(sn) >= 3:  # guard against 1-2 char false synonym hits
        tiers += [("exact:shortname", syn_exact, sn.strip().lower()),
                  ("norm:shortname", syn_norm, norm(sn))]
    for basis, idx, key in tiers:
        if key and key in idx:
            return idx[key], basis
    return set(), "unmatched"


def resolve(query_stems, rids, rec):
    """Pick one rampId for the metabolite.

    Candidates whose canonical name matches the query are strongly
    preferred — this drops synonym-noise (e.g. a glycan carrying a stray
    amino-acid synonym with a huge pathwayCount). Within the name-matched
    set, pathwayCount selects the physiological stereoisomer (L-amino
    acids, D-sugars are the most pathway-connected forms).
    """
    have = [(r, rec[r]) for r in rids if r in rec]
    name_hit = [(r, x) for r, x in have if stem(x.get("cn", "")) in query_stems]
    if name_hit:
        # canonical-name verified: candidates are stereo/charge forms of the
        # same compound. Among kegg-bearing forms the lowest C-number is the
        # canonical central metabolite (L-amino acids, D-sugars). KEGG ids are
        # fixed-width C#####, so lexical min == numeric min.
        kegg_pool = [(r, x) for r, x in name_hit if x["kegg"]]
        if kegg_pool:
            kegg_pool.sort(key=lambda rx: (min(rx[1]["kegg"]), -rx[1]["pwc"], rx[0]))
            return kegg_pool[0][0]
        pool = [(r, x) for r, x in name_hit if x["hmdb"]] or name_hit
        pool.sort(key=lambda rx: (-rx[1]["pwc"], -rx[1]["rank"], rx[0]))
        return pool[0][0]
    # no canonical-name match (matched via looser synonym): names unverified,
    # so fall back to pathway connectivity rather than C-number.
    kegg_pool = [(r, x) for r, x in have if x["kegg"]]
    pool = kegg_pool or [(r, x) for r, x in have if x["hmdb"]] or have
    if not pool:
        return None
    pool.sort(key=lambda rx: (-rx[1]["pwc"], -rx[1]["rank"], rx[0]))
    return pool[0][0]


def make_set_name(domain, name):
    norm_name = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").upper()
    return f"{domain.upper()}_{norm_name}"


def build_ramp_sets(con, id_rows):
    """Map RaMP pathways (reactome/wiki/hmdb) onto the matched Biocrates
    panel. Returns domain -> {set_name: set(shortname)}, the source-id per
    set, and per-metabolite evidence."""
    rid2names = defaultdict(list)
    meta = {}
    for sn, fn, cls, ramp, kegg, hmdb, *_rest in id_rows:
        if ramp:
            rid2names[ramp].append(sn)
            meta[sn] = (fn, kegg, hmdb, ramp)
    c = con.cursor()
    pw = {}
    for prid, sid, typ, cat, name in c.execute(
        "SELECT pathwayRampId, sourceId, type, pathwayCategory, pathwayName "
        "FROM pathway"):
        d = PATHWAY_TYPE_DOMAIN.get(typ)
        if d and name:
            pw[prid] = (d, name, sid)
    sets = defaultdict(lambda: defaultdict(set))
    srcid = {}
    for rid, prid in c.execute(
        "SELECT rampId, pathwayRampId FROM analytehaspathway"):
        if rid in rid2names and prid in pw:
            d, name, sid = pw[prid]
            setn = make_set_name(d, name)
            for sn in rid2names[rid]:
                sets[d][setn].add(sn)
            srcid.setdefault((d, setn), sid)
    return sets, srcid, meta


def write_domain(domain, sets_d, srcid, meta):
    """Filter to >=MIN_MEMBERS, dedupe identical membership, write master
    TSV + GMT. Returns (n_sets, n_dropped_small, n_dropped_dup)."""
    kept = {s: m for s, m in sets_d.items() if len(m) >= MIN_MEMBERS}
    dropped_small = len(sets_d) - len(kept)
    by_members = defaultdict(list)
    for s, m in kept.items():
        by_members[frozenset(m)].append(s)
    final, dropped_dup = {}, 0
    for members, names in by_members.items():
        final[sorted(names)[0]] = members          # lexically-first representative
        dropped_dup += len(names) - 1
    master = os.path.join(HERE, f"{domain}_master.tsv")
    gmt = os.path.join(HERE, f"{domain}_metabolites.gmt")
    with open(master, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["set_name", "shortname", "fullname", "direction",
                    "pathway_source_id", "kegg_id", "hmdb_id", "ramp_id"])
        for s in sorted(final):
            for sn in sorted(final[s]):
                fn, kegg, hmdb, ramp = meta[sn]
                w.writerow([s, sn, fn, "-", srcid.get((domain, s), ""),
                            kegg, hmdb, ramp])
    with open(gmt, "w", newline="", encoding="utf-8") as f:
        for s in sorted(final):
            desc = srcid.get((domain, s), domain)
            f.write("\t".join([s, desc] + sorted(final[s])) + "\n")
    return len(final), dropped_small, dropped_dup


def build_ontology_domain(con, id_rows, ontology_type, prefix, basename):
    """Build a metabolite-level domain from a RaMP ontology category
    (`HMDBOntologyType == ontology_type`): map its terms onto the matched
    panel, keep sets with >=MIN_MEMBERS, dedupe identical membership.
    Writes <basename>_master.tsv + <basename>_metabolites.gmt. Returns
    (n_sets, dropped_small, dropped_dup)."""
    rid2names = defaultdict(list)
    meta = {}
    for sn, fn, cls, ramp, *_rest in id_rows:
        if ramp:
            rid2names[ramp].append(sn)
            meta[sn] = (fn, ramp)
    c = con.cursor()
    terms = {oid: name for oid, name in c.execute(
        "SELECT rampOntologyId, commonName FROM ontology "
        "WHERE HMDBOntologyType=?", (ontology_type,))}
    sets = defaultdict(set)
    oid_of = {}
    for rid, oid in c.execute(
        "SELECT rampCompoundId, rampOntologyId FROM analytehasontology"):
        if rid in rid2names and oid in terms:
            setn = prefix + re.sub(r"[^A-Za-z0-9]+", "_",
                                   terms[oid]).strip("_").upper()
            for sn in rid2names[rid]:
                sets[setn].add(sn)
            oid_of.setdefault(setn, oid)
    kept = {s: m for s, m in sets.items() if len(m) >= MIN_MEMBERS}
    dropped_small = len(sets) - len(kept)
    by_members = defaultdict(list)
    for s, m in kept.items():
        by_members[frozenset(m)].append(s)
    final, dropped_dup = {}, 0
    for members, names in by_members.items():
        final[sorted(names)[0]] = members          # lexically-first representative
        dropped_dup += len(names) - 1
    master = os.path.join(HERE, f"{basename}_master.tsv")
    gmt = os.path.join(HERE, f"{basename}_metabolites.gmt")
    with open(master, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["set_name", "shortname", "fullname", "direction",
                    "ramp_id", "ontology_id"])
        for s in sorted(final):
            for sn in sorted(final[s]):
                fn, ramp = meta[sn]
                w.writerow([s, sn, fn, "-", ramp, oid_of[s]])
    with open(gmt, "w", newline="", encoding="utf-8") as f:
        for s in sorted(final):
            f.write("\t".join([s, oid_of[s]] + sorted(final[s])) + "\n")
    return len(final), dropped_small, dropped_dup


def main():
    con = sqlite3.connect(DB)
    rec, extra, syn_exact, syn_norm = load_ramp(con)
    print(f"RaMP compounds with kegg/hmdb: {len(rec)}  "
          f"syn_exact={len(syn_exact)} syn_norm={len(syn_norm)}")

    rows = list(csv.DictReader(open(BIOC, newline="", encoding="utf-8"),
                               delimiter="\t"))
    id_rows, log_rows = [], []
    for r in rows:
        sn, fn, cls = r["shortname"], r["fullname"], r["analyte_class"]
        rids, basis = candidates(sn, fn, syn_exact, syn_norm)
        query_stems = {stem(fn), stem(no_paren(fn))} - {""}
        chosen = resolve(query_stems, rids, rec) if rids else None
        kegg = hmdb = ramp = pubchem = chebi = refmet = lipidmaps = ""
        if chosen:
            x = rec[chosen]
            ramp = chosen
            kegg = sorted(x["kegg"])[0] if x["kegg"] else ""
            hmdb = _pick_hmdb(x["hmdb"])
            e = extra.get(chosen, {})
            pubchem = _pick_id("pubchem", e.get("pubchem", ()))
            chebi = _pick_id("chebi", e.get("chebi", ()))
            refmet = _pick_id("refmet", e.get("refmet", ()))
            lipidmaps = _pick_id("lipidmaps", e.get("lipidmaps", ()))
        else:
            basis = "unmatched"
        id_rows.append([sn, fn, cls, ramp, kegg, hmdb,
                        pubchem, chebi, refmet, lipidmaps, basis])
        log_rows.append([sn, fn, cls, basis, len(rids),
                         ";".join(sorted(rids)[:6]),
                         "multi-candidate" if len(rids) > 1 else ""])

    with open(ID_OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["short_name", "fullname", "analyte_class",
                    "ramp_id", "kegg_id", "hmdb_id", "pubchem_id",
                    "chebi_id", "refmet_id", "lipidmaps_id", "match_basis"])
        w.writerows(id_rows)
    with open(LOG_OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["short_name", "fullname", "analyte_class", "match_basis",
                    "n_candidates", "candidate_ramp_ids", "note"])
        w.writerows(log_rows)

    # summary
    def has(i, col):
        return id_rows[i][col] != ""
    sm = [i for i, r in enumerate(id_rows) if r[2] not in LIPID_CLASSES]
    lip = [i for i, r in enumerate(id_rows) if r[2] in LIPID_CLASSES]
    kegg_col, hmdb_col = 4, 5
    print(f"\nwrote {ID_OUT} ({len(id_rows)} rows), {LOG_OUT}")
    for label, idx in (("small-mol", sm), ("lipid", lip)):
        n = len(idx)
        k = sum(1 for i in idx if has(i, kegg_col))
        h = sum(1 for i in idx if has(i, hmdb_col))
        a = sum(1 for i in idx if has(i, kegg_col) or has(i, hmdb_col))
        print(f"  {label}: {n}  any-id={a}  kegg={k}  hmdb={h}")
    for col, nm in ((6, "pubchem"), (7, "chebi"), (8, "refmet"),
                    (9, "lipidmaps")):
        print(f"  {nm:10s} populated: "
              f"{sum(1 for i in range(len(id_rows)) if has(i, col))}")
    ambig = sum(1 for r in log_rows if r[6] == "multi-candidate")
    print(f"  ambiguous (multi-candidate): {ambig}")

    # ---- pathway sets (reactome / wikipathways / smpdb) ----
    sets, srcid, meta = build_ramp_sets(con, id_rows)
    print("\n  domain sets (>=%d members, deduped):" % MIN_MEMBERS)
    for d in ("reactome", "wikipathways", "smpdb"):
        n, ds, dd = write_domain(d, sets.get(d, {}), srcid, meta)
        print(f"    {d:13s} {n:4d} sets  (dropped <{MIN_MEMBERS}: {ds}, "
              f"dup-membership: {dd})")

    # ---- RaMP ontology domains (source: origin; health: disease) ----
    for dom, otype, prefix in (("source", "Source", "SOURCE_"),
                               ("health", "Health condition", "HEALTH_")):
        n, ds, dd = build_ontology_domain(con, id_rows, otype, prefix, dom)
        print(f"    {dom:13s} {n:4d} sets  (dropped <{MIN_MEMBERS}: {ds}, "
              f"dup-membership: {dd})")
    con.close()


if __name__ == "__main__":
    main()

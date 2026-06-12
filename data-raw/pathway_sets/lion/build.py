#!/usr/bin/env python3
"""Build LION-derived artifacts for MetaboSetR.

Step 2 (this stage): match Biocrates lipids to LION canonical names and
emit a per-lipid lion_id annotation table plus an assignment log. The
LION term -> lipid-set domain is built in a later stage and reuses the
matching here.

Source: data-raw/pathway_sets/lion/raw/all-LION-lipid-associations.csv
(gitignored). Run from project root:
  python3 data-raw/pathway_sets/lion/build.py
"""
import csv
import os
import re

ROOT = os.getcwd()
HERE = os.path.join(ROOT, "data-raw/pathway_sets/lion")
LION = os.path.join(HERE, "raw/all-LION-lipid-associations.csv")
BIOC = os.path.join(ROOT, "data-raw/reference/biocrates_Quant1000_metabolites.tsv")
ID_OUT = os.path.join(HERE, "lion_lipid_ids.tsv")
LOG_OUT = os.path.join(HERE, "assignment_log.tsv")
SET_MASTER = os.path.join(HERE, "lion_master.tsv")
SET_GMT = os.path.join(HERE, "lion_metabolites.gmt")
MIN_MEMBERS = 3

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

CLASS_MAP = {
    "PC": "PC", "PE": "PE", "PG": "PG", "PI": "PI", "PS": "PS", "PA": "PA",
    "SM": "SM", "CE": "CE", "DG": "DG", "TG": "TG", "MG": "MG",
    "Cer": "Cer", "Hex-Cer": "HexCer", "Hex2Cer": "Hex2Cer",
    "Hex3Cer": "Hex3Cer",
}
# LION writes lysophospholipids as parent class with a 0:0 sn-2 chain,
# e.g. LPC 16:0 -> PC(16:0/0:0).
LYSO_PARENT = {"LPC": "PC", "LPE": "PE", "LPG": "PG", "LPI": "PI",
               "LPS": "PS", "LPA": "PA"}


def nl(s):
    return s.lower().replace(" ", "")


def candidates(short):
    """Ordered (candidate_LION_name, form) pairs for a Biocrates lipid."""
    parts = short.split(" ", 1)
    if len(parts) != 2:
        return []
    ctok, spec = parts
    m = re.search(r"\b([OP])-", spec)
    ether = (m.group(1) + "-") if m else ""
    pairs = re.findall(r"(\d+):(\d+)", spec)
    if not pairs:
        return []
    sc = sum(int(c) for c, _ in pairs)
    sd = sum(int(d) for _, d in pairs)
    if ctok in LYSO_PARENT:                     # lyso: only the 0:0 form
        p = LYSO_PARENT[ctok]
        return [(f"{p}({ether}{sc}:{sd}/0:0)", "lyso")]
    lc = CLASS_MAP.get(ctok)
    if not lc:
        return []
    summ = (f"{lc}({ether}{sc}:{sd})", "sum")
    species = (f"{lc}({spec.strip()})", "species")
    # species-resolved input (sphingoid base "/" acyl) -> prefer species form
    return [species, summ] if "/" in spec else [summ, species]


def load_lion_rows():
    return list(csv.reader(open(LION, newline="")))


def name_index(rows):
    idx = {}
    for r in rows[2:]:
        if r and r[0]:
            idx.setdefault(nl(r[0]), r[0])
    return idx


def make_set_name(term_name):
    n = re.sub(r"\[[^\]]*\]", "", term_name)        # drop [GP01] class codes
    n = re.sub(r"[^A-Za-z0-9]+", "_", n).strip("_").upper()
    return "LION_" + n


def build_lion_sets(rows, id_rows):
    """term column -> set(Biocrates shortname) via matched lion_id rows."""
    from collections import defaultdict
    lid2names = defaultdict(list)
    full = {}
    for sn, fn, cls, lion_id, form in id_rows:
        if lion_id:
            lid2names[lion_id].append(sn)
            full[sn] = fn
    matched = set(lid2names)
    term_ids, term_names = rows[0][1:], rows[1][1:]
    sets = defaultdict(set)
    for row in rows[2:]:
        if not row or row[0] not in matched:
            continue
        names = lid2names[row[0]]
        for j, cell in enumerate(row[1:]):
            if cell.strip().lower() == "x":
                sets[j].update(names)
    return sets, term_ids, term_names, full


def write_lion_sets(sets, term_ids, term_names, full):
    """Filter >=MIN_MEMBERS, dedupe identical membership, write master+GMT."""
    from collections import defaultdict
    kept = {j: m for j, m in sets.items() if len(m) >= MIN_MEMBERS}
    dropped_small = len(sets) - len(kept)
    by_members = defaultdict(list)
    for j, m in kept.items():
        by_members[frozenset(m)].append(j)
    final, used, dropped_dup = {}, set(), 0
    for members, cols in by_members.items():
        j = sorted(cols, key=lambda c: term_ids[c])[0]   # lowest term id
        dropped_dup += len(cols) - 1
        name = make_set_name(term_names[j])
        if name in used:
            name = f"{name}_{term_ids[j].replace(':', '')}"
        used.add(name)
        final[name] = (members, term_ids[j])
    with open(SET_MASTER, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["set_name", "shortname", "fullname", "direction",
                    "lion_term_id"])
        for s in sorted(final):
            members, tid = final[s]
            for sn in sorted(members):
                w.writerow([s, sn, full[sn], "-", tid])
    with open(SET_GMT, "w", newline="", encoding="utf-8") as f:
        for s in sorted(final):
            members, tid = final[s]
            f.write("\t".join([s, tid] + sorted(members)) + "\n")
    return len(final), dropped_small, dropped_dup


def main():
    rows = load_lion_rows()
    idx = name_index(rows)
    print(f"LION lipids (norm keys): {len(idx)}")

    bio = list(csv.DictReader(open(BIOC, newline=""), delimiter="\t"))
    id_rows, log_rows = [], []
    matched = 0
    for r in bio:
        sn, fn, cls = r["shortname"], r["fullname"], r["analyte_class"]
        lion_id, form = "", ""
        if cls in LIPID_CLASSES:
            for cand, frm in candidates(sn):
                if nl(cand) in idx:
                    lion_id, form = idx[nl(cand)], frm
                    break
        if lion_id:
            matched += 1
        id_rows.append([sn, fn, cls, lion_id, form])
        if cls in LIPID_CLASSES:
            log_rows.append([sn, fn, cls, lion_id or "<unmatched>", form or "-",
                             ";".join(c for c, _ in candidates(sn))])

    with open(ID_OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["short_name", "fullname", "analyte_class",
                    "lion_id", "match_form"])
        w.writerows(id_rows)
    with open(LOG_OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["short_name", "fullname", "analyte_class",
                    "lion_id", "match_form", "candidates_tried"])
        w.writerows(log_rows)

    n_lip = sum(1 for r in bio if r["analyte_class"] in LIPID_CLASSES)
    print(f"wrote {ID_OUT} ({len(id_rows)} rows), {LOG_OUT}")
    print(f"lipids matched to LION: {matched}/{n_lip}")

    # ---- LION term -> lipid-set domain ----
    sets, term_ids, term_names, full = build_lion_sets(rows, id_rows)
    n, ds, dd = write_lion_sets(sets, term_ids, term_names, full)
    print(f"lion sets (>={MIN_MEMBERS} members, deduped): {n}  "
          f"(dropped <{MIN_MEMBERS}: {ds}, dup-membership: {dd})")


if __name__ == "__main__":
    main()

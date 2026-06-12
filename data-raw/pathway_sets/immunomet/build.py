#!/usr/bin/env python3
"""
Build immunometabolism pathway sets based on O'Neill 2016 framework.

Outputs (in this directory):
  - immunomet_metabolites_master.tsv  metabolite-level master (from metabolites)
  - immunomet_metabolites.gmt         metabolite-level GMT export
  - assignment_log_metabolites.tsv    audit log for metabolite assignment

Design notes:
  * Strategy C (metabolites): A metabolite is assigned to a pathway if its
    shortname/fullname matches PATHWAY_METABOLITES for that pathway, OR if
    its analyte_class matches PATHWAY_CLASSES for that pathway.
    Overlap allowed (one metabolite may be in multiple sets).
"""
import csv
import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathways import (
    PATHWAY_METABOLITES, PATHWAY_CLASSES, PATHWAY_EXCLUDE,
    PATHWAY_ORDER, PATHWAY_DESCRIPTIONS,
)

# Resolve paths relative to this script so the build is portable.
ROOT = Path(__file__).resolve().parents[3]
METABOLITES_TSV = str(ROOT / 'data-raw' / 'reference' / 'biocrates_Quant1000_metabolites.tsv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# Helpers
# =============================================================================

def load_metabolites():
    with open(METABOLITES_TSV) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def build_panel_index(metabs):
    """Return {shortname: metab_row, fullname: metab_row, and class lookup}."""
    by_name = {}
    by_class = {}
    for m in metabs:
        by_name[m['shortname']] = m
        by_name[m['fullname']] = m
        by_class.setdefault(m['analyte_class'], []).append(m)
    return by_name, by_class


# =============================================================================
# Metabolite-level assignment (Strategy C)
# =============================================================================

def assign_metabolite_pathways(metab, by_name):
    """Return list of pathway names this metabolite belongs to."""
    pathways = []
    for pw in PATHWAY_ORDER:
        individual = PATHWAY_METABOLITES[pw]
        classes = PATHWAY_CLASSES[pw]
        excluded = PATHWAY_EXCLUDE[pw]

        # Skip if metabolite is on the exclusion list for this pathway
        if metab['shortname'] in excluded or metab['fullname'] in excluded:
            continue

        if metab['shortname'] in individual or metab['fullname'] in individual:
            pathways.append(pw)
        elif metab['analyte_class'] in classes:
            pathways.append(pw)
    return pathways


# =============================================================================
# Main
# =============================================================================

def main():
    metabs = load_metabolites()
    by_name, by_class = build_panel_index(metabs)

    # ---------- Metabolite-level ----------
    metab_rows = []        # long-format master rows
    metab_log_rows = []    # one row per metabolite, with pathways joined
    for m in metabs:
        pws = assign_metabolite_pathways(m, by_name)
        metab_log_rows.append({
            'shortname': m['shortname'],
            'fullname': m['fullname'],
            'analyte_class': m['analyte_class'],
            'n_pathways': len(pws),
            'pathways': '; '.join(pws) if pws else '',
        })
        for pw in pws:
            metab_rows.append({
                'shortname': m['shortname'],
                'fullname': m['fullname'],
                'analyte_class': m['analyte_class'],
                'set_name': pw,
                'match_basis': 'individual' if (m['shortname'] in PATHWAY_METABOLITES[pw]
                                                or m['fullname'] in PATHWAY_METABOLITES[pw])
                                else f"class:{m['analyte_class']}",
            })

    # ---------- Write outputs ----------
    write_tsv(os.path.join(OUT_DIR, 'immunomet_metabolites_master.tsv'),
              ['shortname', 'fullname', 'analyte_class', 'set_name', 'match_basis'],
              metab_rows)
    write_tsv(os.path.join(OUT_DIR, 'assignment_log_metabolites.tsv'),
              ['shortname', 'fullname', 'analyte_class', 'n_pathways', 'pathways'],
              metab_log_rows)

    # GMT export
    write_gmt(os.path.join(OUT_DIR, 'immunomet_metabolites.gmt'),
              metab_rows, 'shortname')

    # ---------- Summary ----------
    print(f"Source metabolites: {len(metabs)}")
    print()
    print("=== METABOLITE-LEVEL ===")
    from collections import Counter
    metab_set_counts = Counter(r['set_name'] for r in metab_rows)
    for pw in PATHWAY_ORDER:
        print(f"  {metab_set_counts[pw]:4d}  {pw}")
    # overlap stats
    unique_metabs_in_any_set = {r['shortname'] for r in metab_rows}
    print(f"  ----")
    print(f"  Unique metabolites assigned to at least one pathway: {len(unique_metabs_in_any_set)}")
    print(f"  Total (with multi-set membership): {len(metab_rows)}")
    overlap_metabs = [k for k, v in Counter(r['shortname'] for r in metab_rows).items() if v > 1]
    print(f"  Metabolites in more than one pathway: {len(overlap_metabs)}")


def write_tsv(path, fieldnames, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                          quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def write_gmt(path, long_rows, member_key):
    """Write GMT: one line per set, tab-separated set-name, description, members."""
    from collections import defaultdict
    sets = defaultdict(list)
    for r in long_rows:
        sets[r['set_name']].append(r[member_key])
    with open(path, 'w') as f:
        for pw in PATHWAY_ORDER:
            members = sets.get(pw, [])
            if not members:
                continue
            desc = PATHWAY_DESCRIPTIONS.get(pw, '').replace('\t', ' ').replace('\n', ' ')
            line = [pw, desc] + members
            f.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()

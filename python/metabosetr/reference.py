"""Reference data used by readers, annotation, and pathway loaders.

The R package bundles these in ``R/sysdata.rda``; here they are built at
import time from the committed source TSVs (the same inputs
``data-raw/make_sysdata.R`` uses), so no ``sysdata.rda`` export is required.

* :func:`metabolite_dict` -- canonical name table joined with RaMP / LION
  cross-reference IDs (reference TSV + ``ramp_metabolite_ids.tsv`` +
  ``lion_lipid_ids.tsv``), matched by ``short_name``.
* :data:`STATUS_COLOR_MAP` -- WebIDQ cell-fill hex -> status text
  (DESIGN.md §2.4), hard-coded exactly as the R sysdata hard-codes it.

Bundling these into the installed package (vs. reading the repo TSVs) is
deferred -- see MIGRATION.md §8 Q5.
"""

from __future__ import annotations

import functools
from pathlib import Path

import pandas as pd

# DESIGN.md §2.4 -- WebIDQ cell-fill hex (uppercase, no '#') -> status text.
STATUS_COLOR_MAP: dict[str, str] = {
    "B9DE83": "Valid",
    "BBA7B9": "< threshold",
    "A28BA3": "< LOD",
    "B2D1DC": "< LLOQ",
    "7FB2C5": "> ULOQ",
}

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DATA_RAW = _REPO_ROOT / "data-raw"
_REFERENCE_TSV = _DATA_RAW / "reference" / "biocrates_Quant1000_metabolites.tsv"
_RAMP_IDS_TSV = _DATA_RAW / "pathway_sets" / "ramp" / "ramp_metabolite_ids.tsv"
_LION_IDS_TSV = _DATA_RAW / "pathway_sets" / "lion" / "lion_lipid_ids.tsv"

# Installed GMTs (mirror ``pathway_set_gmt_path``).
EXTDATA_PATHWAY_SETS = _REPO_ROOT / "inst" / "extdata" / "pathway_sets"


def _read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV the way make_sysdata.R does: empty cells -> NA, all str."""
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_values=[""],
        quoting=3,  # csv.QUOTE_NONE; masters have no quoting
    )


@functools.lru_cache(maxsize=1)
def metabolite_dict() -> pd.DataFrame:
    """Canonical name table joined with RaMP / LION IDs (1234 rows).

    Columns: ``short_name``, ``long_name``, ``class``, ``kegg_id``,
    ``hmdb_id``, ``pubchem_id``, ``chebi_id``, ``refmet_id``,
    ``lipidmaps_id``, ``lion_id``. Unmatched IDs are ``None``.
    """
    met = _read_tsv(_REFERENCE_TSV)
    out = pd.DataFrame(
        {
            "short_name": met["shortname"],
            "long_name": met["fullname"],
            "class": met["analyte_class"],
        }
    )

    ramp = _read_tsv(_RAMP_IDS_TSV).set_index("short_name")
    lion = _read_tsv(_LION_IDS_TSV).set_index("short_name")
    ramp_cols = [
        "kegg_id",
        "hmdb_id",
        "pubchem_id",
        "chebi_id",
        "refmet_id",
        "lipidmaps_id",
    ]
    for col in ramp_cols:
        out[col] = out["short_name"].map(ramp[col])
    out["lion_id"] = out["short_name"].map(lion["lion_id"])

    # Represent missing uniformly as None (NA from na_values or unmatched).
    return out.where(pd.notna(out), None)

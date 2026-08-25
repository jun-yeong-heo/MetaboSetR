"""Cross-reference ID annotation (port of ``annotate.R``).

Maps Biocrates ``short_name``s to the external-DB IDs in
:func:`metabolite_dict` (DESIGN.md §6.5 / §9). Always one row per input
element, in input order, with ``None`` for unmatched names or missing IDs --
ready to key a ``pathview`` KEGG-map overlay.
"""

from __future__ import annotations

import pandas as pd

from .reference import metabolite_dict

_ANN_COLS = [
    "long_name",
    "class",
    "kegg_id",
    "hmdb_id",
    "pubchem_id",
    "chebi_id",
    "refmet_id",
    "lipidmaps_id",
    "lion_id",
]


def annotate_metabolites(x) -> pd.DataFrame:
    """Annotate metabolite ``short_name``s with cross-reference IDs.

    ``x`` is a sequence of metabolite short names. Returns a DataFrame with
    columns ``metabolite`` (the input), ``long_name``, ``class``, and the ID
    columns; ``None`` where unknown.
    """
    x = list(x)
    md = metabolite_dict().set_index("short_name")
    out = pd.DataFrame({"metabolite": x})
    for col in _ANN_COLS:
        lookup = md[col].to_dict()  # short_name -> id (or None)
        values = [lookup.get(name) for name in x]
        # normalise residual NaN to None, and assign as an object Series so
        # pandas does not coerce None back to NaN on column assignment.
        values = [
            None if (v is None or (isinstance(v, float) and pd.isna(v))) else v
            for v in values
        ]
        out[col] = pd.Series(values, index=out.index, dtype=object)
    return out

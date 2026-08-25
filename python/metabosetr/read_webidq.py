"""WebIDQ concentration + status readers (port of ``read-webidq.R``).

DESIGN.md §2.1-2.4 for the input format; §7.3 for the user-facing
signatures. Helpers are underscore-prefixed.
"""

from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd

from . import _xlsx
from .classes import ConcentrationData, QCData
from .reference import metabolite_dict
from .read_status_color import _read_status_color

# Header columns that must precede the metabolite columns, in Biocrates order.
META_COLS = [
    "Sample identification",
    "Sample description",
    "Submission name",
    "Collection date",
    "Species",
    "Material",
]

_KIT_RE = re.compile(r"KIT\d+")


def _is_padding_row(sample_id) -> bool:
    """True for vendor padding rows: ``Sample identification`` is NA.

    Matches ``.is_padding_row`` in R -- the class-header row and trailing
    blank rows leave this cell empty while real samples always fill it.
    """
    if sample_id is None:
        return True
    if isinstance(sample_id, float) and np.isnan(sample_id):
        return True
    return False


def _rows_to_frame(rows) -> pd.DataFrame:
    """Build a DataFrame from (header, data...) rows, header = first row."""
    if not rows:
        raise ValueError("Empty sheet.")
    header = [("" if v is None else str(v)) for v in rows[0]]
    width = len(header)
    data = []
    for r in rows[1:]:
        r = list(r) + [None] * (width - len(r))
        data.append(r[:width])
    return pd.DataFrame(data, columns=header)


def _read_conc_sheet(path: str) -> dict:
    """Read the concentration sheet: (sample_meta, assay, metabolite_names)."""
    # Sheet name is not validated; first sheet only (decisions.md #50).
    _, rows = _xlsx.first_sheet_rows(path)
    raw = _rows_to_frame(rows)

    missing = [c for c in META_COLS if c not in raw.columns]
    if missing:
        raise ValueError(
            f"Concentration file '{path}' is missing required meta columns: "
            + ", ".join(missing)
        )

    keep = ~raw["Sample identification"].map(_is_padding_row)
    raw = raw[keep].reset_index(drop=True)

    metabolite_names = [c for c in raw.columns if c not in META_COLS]
    sample_meta = raw[META_COLS].reset_index(drop=True).copy()
    assay = (
        raw[metabolite_names]
        .apply(pd.to_numeric, errors="coerce")
        .to_numpy(dtype=float)
    )
    return {
        "sample_meta": sample_meta,
        "assay": assay,
        "metabolite_names": metabolite_names,
    }


def _read_status_text(path: str, conc_data: dict) -> np.ndarray:
    """Read a parallel status-text sheet aligned to ``conc_data`` columns."""
    _, rows = _xlsx.first_sheet_rows(path)
    raw = _rows_to_frame(rows)

    meta_missing = [c for c in META_COLS if c not in raw.columns]
    if meta_missing:
        raise ValueError(
            f"Status file '{path}' is missing required meta columns: "
            + ", ".join(meta_missing)
        )

    raw = raw[~raw["Sample identification"].map(_is_padding_row)].reset_index(drop=True)
    n_conc = len(conc_data["sample_meta"])
    if len(raw) != n_conc:
        raise ValueError(
            f"Status file row count ({len(raw)}) differs from concentration "
            f"file ({n_conc}); both files were processed with the same "
            "blank-row drop rule, so this means the two files cover different "
            "sample sets."
        )

    missing = [c for c in conc_data["metabolite_names"] if c not in raw.columns]
    if missing:
        raise ValueError(
            "Status file is missing metabolite columns: "
            + ", ".join(missing[:10])
            + (", ..." if len(missing) > 10 else "")
        )

    status = raw[conc_data["metabolite_names"]].to_numpy(dtype=object)
    return status


def _extract_kit_id(submission_names) -> list:
    """Extract ``KIT\\d+`` from each Submission name; NA/no-match -> None."""
    out = []
    for name in submission_names:
        if name is None or (isinstance(name, float) and np.isnan(name)):
            out.append(None)
            continue
        m = _KIT_RE.search(str(name))
        out.append(m.group(0) if m else None)
    return out


def _resolve_metabolite_names(input_names, synonyms: Optional[dict] = None) -> list:
    """Resolve names to canonical short names (DESIGN.md §2.7, 3 passes)."""
    md = metabolite_dict()
    canonical = set(md["short_name"])
    fullname_map = dict(zip(md["long_name"], md["short_name"]))

    def resolve_one(nm: str):
        if nm in canonical:
            return nm
        if nm in fullname_map:
            return fullname_map[nm]
        nm2 = re.sub(r"^Total\s+", "", re.sub(r"\s+", " ", nm.strip()))
        if nm2 in canonical:
            return nm2
        if nm2 in fullname_map:
            return fullname_map[nm2]
        if synonyms and nm in synonyms:
            mapped = synonyms[nm]
            if mapped in canonical:
                return mapped
            if mapped in fullname_map:
                return fullname_map[mapped]
        return None

    resolved = [resolve_one(nm) for nm in input_names]
    bad = [nm for nm, r in zip(input_names, resolved) if r is None]
    if bad:
        raise ValueError(
            "Could not resolve these metabolite names to canonical short names:\n  "
            + ", ".join(bad[:10])
            + (", ..." if len(bad) > 10 else "")
            + "\nProvide them via the `synonyms={typed: canonical}` argument."
        )
    return resolved


def _read_into_matrices(
    concentration_file: str,
    status_file: Optional[str],
    status_method: str,
    synonyms: Optional[dict],
) -> dict:
    conc_data = _read_conc_sheet(concentration_file)
    if status_method == "text":
        if status_file is None:
            raise ValueError("status_file is required when status_method='text'.")
        status = _read_status_text(status_file, conc_data)
    else:
        status = _read_status_color(concentration_file, conc_data)

    conc_data["sample_meta"]["kit_id"] = _extract_kit_id(
        conc_data["sample_meta"]["Submission name"].tolist()
    )
    metabolite_names = _resolve_metabolite_names(
        conc_data["metabolite_names"], synonyms
    )
    return {
        "assay": conc_data["assay"],
        "status": status,
        "sample_meta": conc_data["sample_meta"],
        "metabolite_names": metabolite_names,
    }


def read_webidq(
    concentration_file: str,
    status_file: Optional[str] = None,
    status_method: str = "text",
    synonyms: Optional[dict] = None,
) -> ConcentrationData:
    """Read a WebIDQ concentration export into a :class:`ConcentrationData`.

    Use ``status_method='text'`` (default) with the status xlsx, or
    ``'color'`` to recover status from the concentration file's cell fills.
    """
    if status_method not in ("text", "color"):
        raise ValueError("status_method must be 'text' or 'color'.")
    parts = _read_into_matrices(concentration_file, status_file, status_method, synonyms)
    return ConcentrationData(
        assay=parts["assay"],
        status=parts["status"],
        sample_meta=parts["sample_meta"],
        metabolite_names=parts["metabolite_names"],
    )


def read_webidq_qc(
    concentration_file: str,
    status_file: Optional[str] = None,
    status_method: str = "text",
    synonyms: Optional[dict] = None,
) -> QCData:
    """Read a WebIDQ pooled-QC export into a :class:`QCData`."""
    if status_method not in ("text", "color"):
        raise ValueError("status_method must be 'text' or 'color'.")
    parts = _read_into_matrices(concentration_file, status_file, status_method, synonyms)
    return QCData(
        assay=parts["assay"],
        status=parts["status"],
        sample_meta=parts["sample_meta"],
        metabolite_names=parts["metabolite_names"],
    )

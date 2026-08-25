"""Optional per-kit LOD/LLOQ/ULOQ reader (port of ``read-thresholds.R``).

DESIGN.md §2.5. Thresholds are metadata only -- they are not used to
recompute statuses. Kit ID comes from the filename (``KIT\\d+``); the file
itself does not carry it.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np

from . import _xlsx
from .classes import ThresholdData
from .read_webidq import _resolve_metabolite_names

_KIT_RE = re.compile(r"KIT\d+")


def _kit_id_from_filename(path: str) -> str:
    base = Path(path).name
    m = _KIT_RE.search(base)
    if not m:
        raise ValueError(
            f"Could not find a `KIT\\d+` token in the threshold file name "
            f"'{base}'. Rename the file to include the kit ID "
            "(e.g. `thresholds_KIT01.xlsx`)."
        )
    return m.group(0)


def _to_float_row(cells) -> np.ndarray:
    out = np.empty(len(cells), dtype=float)
    for i, v in enumerate(cells):
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            out[i] = np.nan
    return out


def _read_one_threshold_file(path: str) -> dict:
    """Read one per-kit threshold sheet (4 x (1+N) layout)."""
    _, rows = _xlsx.first_sheet_rows(path)
    if len(rows) < 4:
        raise ValueError(
            f"Threshold file '{path}' has fewer than 4 rows; expected "
            "Measurement-time / LOD / ULOQ / LLOQ rows."
        )
    metabolite_names = [str(v) for v in rows[0][1:]]
    lod = _to_float_row(rows[1][1:])
    uloq = _to_float_row(rows[2][1:])
    lloq = _to_float_row(rows[3][1:])
    return {
        "kit_id": _kit_id_from_filename(path),
        "lod": lod,
        "uloq": uloq,
        "lloq": lloq,
        "metabolite_names": metabolite_names,
    }


def read_thresholds(threshold_files, synonyms: Optional[dict] = None) -> ThresholdData:
    """Read per-kit LOD/LLOQ/ULOQ files into a :class:`ThresholdData`.

    One file per kit; supply one path per kit and they are stacked by kit.
    """
    if isinstance(threshold_files, str):
        threshold_files = [threshold_files]
    if len(threshold_files) == 0:
        raise ValueError("threshold_files must contain at least one path.")

    per_kit = [_read_one_threshold_file(p) for p in threshold_files]

    ref_metabs = per_kit[0]["metabolite_names"]
    for pk in per_kit:
        if pk["metabolite_names"] != ref_metabs:
            raise ValueError("Threshold files do not share the same metabolite columns.")

    resolved = _resolve_metabolite_names(ref_metabs, synonyms)
    kit_ids = [pk["kit_id"] for pk in per_kit]

    def to_matrix(field: str) -> np.ndarray:
        return np.vstack([pk[field] for pk in per_kit])

    return ThresholdData(
        lod=to_matrix("lod"),
        lloq=to_matrix("lloq"),
        uloq=to_matrix("uloq"),
        kit_ids=kit_ids,
        metabolite_names=resolved,
    )

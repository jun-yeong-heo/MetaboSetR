"""Color-based status fallback (port of ``read-status-color.R``).

DESIGN.md §2.4. Recovers the status matrix from cell fill colors of the
concentration xlsx when no dedicated status file is available. ``tidyxl``'s
format table is replaced by ``openpyxl``'s per-cell fill, which exposes the
fill color directly.
"""

from __future__ import annotations

import numpy as np

from . import _xlsx
from .reference import STATUS_COLOR_MAP


def _is_padding_row(sample_id) -> bool:
    if sample_id is None:
        return True
    return isinstance(sample_id, float) and np.isnan(sample_id)


def _hex_to_status(hex6):
    """Map a 6-char hex (no '#', no alpha) to status text; None if unknown."""
    if hex6 is None or len(hex6) != 6:
        return None
    return STATUS_COLOR_MAP.get(hex6.upper())


def _cell_status(cell):
    """Status from a cell's fill fgColor, or None for unrecognised/blank."""
    rgb = getattr(cell.fill.fgColor, "rgb", None)
    if not isinstance(rgb, str):
        return None  # theme/indexed color -> no recognised status
    if len(rgb) == 8:  # ARGB, strip alpha -> RGB
        hex6 = rgb[2:8]
    elif len(rgb) == 6:
        hex6 = rgb
    else:
        return None
    return _hex_to_status(hex6)


def _read_status_color(path: str, conc_data: dict) -> np.ndarray:
    """Read the status matrix from cell fills, aligned to ``conc_data``.

    Cells without a recognised fill default to ``"Valid"`` (DESIGN.md §2.4).
    """
    wb = _xlsx.load_workbook(path)
    ws = wb[wb.sheetnames[0]]

    header = {}
    for c in range(1, ws.max_column + 1):
        v = ws.cell(1, c).value
        if v is not None:
            header[str(v)] = c
    if "Sample identification" not in header:
        raise ValueError(
            f"Could not locate `Sample identification` header column in "
            f"'{path}' for the color-status fallback."
        )
    sid_col = header["Sample identification"]

    kept_rows = [
        r
        for r in range(2, ws.max_row + 1)
        if not _is_padding_row(ws.cell(r, sid_col).value)
    ]
    n_samples = len(conc_data["sample_meta"])
    if len(kept_rows) != n_samples:
        raise ValueError(
            f"Color-status reader saw {len(kept_rows)} non-padding sample rows "
            f"but the conc reader retained {n_samples}; padding-row rules are "
            "out of sync."
        )

    metabolite_names = conc_data["metabolite_names"]
    out = np.full((n_samples, len(metabolite_names)), "Valid", dtype=object)
    for i, r in enumerate(kept_rows):
        for j, name in enumerate(metabolite_names):
            c = header.get(name)
            if c is None:
                continue
            status = _cell_status(ws.cell(r, c))
            if status is not None:
                out[i, j] = status
    return out

"""log2 transform with an adaptive per-feature pseudocount.

Port of ``transform.R`` (DESIGN.md §4.5). ``log2(x + c)`` where the
pseudocount ``c`` is chosen *per row* (per feature) as half the smallest
finite positive value in that row whenever the row contains a zero, and is
``0`` otherwise.

Conventions (identical to the R version):

* ``x`` is treated row-wise: rows are features, columns are samples. Pass
  ``x.T`` if your input has features as columns.
* A 1-D input is treated as one feature.
* Inf / NaN cells pass through.
* If ``pseudocount`` is given, that scalar is added to every cell.
"""

from __future__ import annotations

import numpy as np


def _log2_one_row(row: np.ndarray) -> np.ndarray:
    row = np.asarray(row, dtype=float)
    has_zero = np.any(row == 0)  # NaN == 0 is False, matching R na.rm behaviour
    with np.errstate(divide="ignore", invalid="ignore"):
        if not has_zero:
            return np.log2(row)
        finite_pos = row[np.isfinite(row) & (row > 0)]
        pc = finite_pos.min() / 2 if finite_pos.size > 0 else 0.0
        return np.log2(row + pc)


def log2_with_pseudocount(x, pseudocount=None):
    """log2-transform ``x`` with the per-feature half-min pseudocount rule.

    Returns an array the same shape as ``x``.
    """
    if pseudocount is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log2(np.asarray(x, dtype=float) + pseudocount)
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 1:
        return _log2_one_row(arr)
    out = np.empty_like(arr, dtype=float)
    for i in range(arr.shape[0]):
        out[i, :] = _log2_one_row(arr[i, :])
    return out

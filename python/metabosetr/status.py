"""Status annotation helpers (port of ``status-helpers.R``).

DESIGN.md §3.2. ``< threshold`` and ``< LLOQ`` count as valid for
filter-logic purposes; the original annotation is preserved on the assay.

NA handling mirrors R logical-with-NA semantics: NA inputs propagate to NA
(``None``) in every predicate except :func:`is_missing`, which returns
``True`` on NA. Inputs are treated element-wise; ``None`` and ``NaN`` are
both taken as missing. Results are returned as object arrays holding
``True`` / ``False`` / ``None`` to represent the three-valued logic.
"""

from __future__ import annotations

import numpy as np

_VALID = {"Valid", "< threshold", "< LLOQ"}


def _is_na(value) -> bool:
    if value is None:
        return True
    return isinstance(value, float) and np.isnan(value)


def _predicate(status, fn, na_result):
    arr = np.asarray(status, dtype=object).ravel()
    out = np.empty(arr.shape, dtype=object)
    for i, v in enumerate(arr):
        out[i] = na_result if _is_na(v) else fn(v)
    return out.reshape(np.asarray(status, dtype=object).shape)


def is_valid(status):
    """status in {Valid, < threshold, < LLOQ}; NA -> None."""
    return _predicate(status, lambda v: v in _VALID, None)


def is_below_lloq(status):
    """status == '< LLOQ'; NA -> None."""
    return _predicate(status, lambda v: v == "< LLOQ", None)


def is_below_lod(status):
    """status == '< LOD'; NA -> None."""
    return _predicate(status, lambda v: v == "< LOD", None)


def is_above_uloq(status):
    """status == '> ULOQ'; NA -> None."""
    return _predicate(status, lambda v: v == "> ULOQ", None)


def is_missing(status):
    """NA or 'Missing' -> True."""
    return _predicate(status, lambda v: v == "Missing", True)

import numpy as np

from metabosetr import status
from metabosetr.transform import log2_with_pseudocount


def test_is_valid_categories_and_na():
    s = np.array(["Valid", "< threshold", "< LLOQ", "< LOD", "> ULOQ", None], dtype=object)
    out = status.is_valid(s)
    assert list(out) == [True, True, True, False, False, None]


def test_is_below_lod_and_uloq():
    s = np.array(["< LOD", "> ULOQ", "Valid", None], dtype=object)
    assert list(status.is_below_lod(s)) == [True, False, False, None]
    assert list(status.is_above_uloq(s)) == [False, True, False, None]


def test_is_missing_na_is_true():
    s = np.array(["Valid", "Missing", None], dtype=object)
    assert list(status.is_missing(s)) == [False, True, True]


def test_log2_per_feature_pseudocount_matches_rule():
    # Mirrors the R docstring example (transform.R).
    x = np.array([[0.0, 0.1, 0.5, 1.0], [100.0, 200.0, 0.0, 400.0]])
    out = log2_with_pseudocount(x)
    # low row: zero present -> pc = min positive / 2 = 0.05
    assert np.allclose(out[0], np.log2(np.array([0.05, 0.15, 0.55, 1.05])))
    # high row: zero present -> pc = 100 / 2 = 50
    assert np.allclose(out[1], np.log2(np.array([150.0, 250.0, 50.0, 450.0])))


def test_log2_no_zero_is_plain_log2():
    row = np.array([1.0, 2.0, 4.0])
    assert np.allclose(log2_with_pseudocount(row), np.log2(row))


def test_log2_explicit_pseudocount_overrides():
    row = np.array([0.0, 1.0, 3.0])
    assert np.allclose(log2_with_pseudocount(row, pseudocount=1.0), np.log2(row + 1.0))

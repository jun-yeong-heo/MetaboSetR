import numpy as np
import pandas as pd

from metabosetr.filter_metabolites import _cv_filter, _cv_pct, _lod_rate_filter
from metabosetr.filter_samples import _apply_outlier_method, _detect_sample_outliers


def test_cv_pct_hand_computed():
    # [10,12,14]: mean=12, sd(ddof=1)=2 -> CV = 100*2/12
    assert np.isclose(_cv_pct(np.array([10.0, 12.0, 14.0])), 100 * 2 / 12)


def test_cv_pct_zero_mean_and_singleton_are_nan():
    assert np.isnan(_cv_pct(np.array([-1.0, 0.0, 1.0])))  # mean 0
    assert np.isnan(_cv_pct(np.array([5.0])))  # sd(ddof=1) undefined


def test_cv_filter_per_kit_is_and_across_kits():
    # metabolite 0 passes both kits; metabolite 1 fails kit2.
    assay = np.array(
        [
            [10.0, 10.0],
            [12.0, 12.0],  # kit A
            [14.0, 14.0],
            [10.0, 1.0],
            [12.0, 50.0],  # kit B: metab1 huge CV
            [14.0, 100.0],
        ]
    )
    meta = pd.DataFrame({"kit_id": ["KIT01"] * 3 + ["KIT02"] * 3})
    res = _cv_filter(assay, meta, cv_threshold=20, cv_mode="per_kit")
    assert list(res["pass"]) == [True, False]
    assert "cv_KIT01" in res["log"].columns and "cv_KIT02" in res["log"].columns


def test_lod_rate_filter_hand_computed():
    status = np.array(
        [["< LOD", "Valid"], ["< LOD", "Valid"], ["Valid", "Valid"], ["Valid", "Valid"]],
        dtype=object,
    )
    res = _lod_rate_filter(status, lod_rate_threshold=0.5)
    assert np.allclose(res["lod_rate"], [0.5, 0.0])
    assert list(res["pass"]) == [True, True]  # 0.5 <= 0.5 passes


def test_iqr_outlier_threshold_formula():
    # type-7 quantiles of [1,2,3,4,100]: q1=2, q3=4, IQR=2, thr=4+1.5*2=7
    rates = np.array([1.0, 2.0, 3.0, 4.0, 100.0])
    flags = _apply_outlier_method(rates, "IQR", 0.95)
    assert list(flags) == [False, False, False, False, True]


def test_detect_sample_outliers_below_lod():
    # sample 5 is all-<LOD -> high lod_rate -> outlier
    rows = [["Valid", "Valid"]] * 5 + [["< LOD", "< LOD"]]
    status = np.array(rows, dtype=object)
    res = _detect_sample_outliers(status, metric="below_lod", method="IQR")
    assert res["outliers"][-1]
    assert not res["outliers"][:-1].any()

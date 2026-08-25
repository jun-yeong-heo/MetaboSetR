"""Best-effort smoke tests for the QC diagnostics layer (CLAUDE.md §3)."""

import numpy as np
from conftest import SYNTHETIC

import metabosetr as m


def _prep():
    D = SYNTHETIC
    conc = m.read_webidq(str(D / "sample_conc.xlsx"), str(D / "sample_status.xlsx"))
    qc = m.read_webidq_qc(str(D / "qc_conc.xlsx"), str(D / "qc_status.xlsx"))
    return m.preprocess(conc, qc)


def test_filter_summary():
    prep = _prep()
    s = m.qc_filter_summary(prep)
    assert list(s["item"]) == ["samples", "metabolites"]
    # totals consistent, pass+fail == total per row
    assert (s["pass"] + s["fail"] == s["total"]).all()
    assert s.loc[s["item"] == "samples", "total"].iloc[0] == 6


def test_pca_plot_runs():
    prep = _prep()
    fig = m.qc_pca_plot(prep, mode="qc_only")
    assert fig is not None
    assert len(fig.axes) == 1


def test_permanova_runs():
    prep = _prep()
    res = m.qc_permanova(prep)  # two kits present in synthetic QC
    assert "test statistic" in res.index
    assert not np.isnan(res["p-value"])


def test_cv_and_lod_plots_run():
    prep = _prep()
    assert m.qc_cv_plot(prep) is not None
    assert m.qc_lod_rate_plot(prep, by="sample") is not None

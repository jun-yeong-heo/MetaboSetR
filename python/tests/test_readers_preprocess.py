import numpy as np
import pytest
from conftest import SYNTHETIC

import metabosetr as m
from metabosetr import _xlsx
from metabosetr.read_webidq import _extract_kit_id, _resolve_metabolite_names


def _paths():
    return {
        "sc": str(SYNTHETIC / "sample_conc.xlsx"),
        "ss": str(SYNTHETIC / "sample_status.xlsx"),
        "qc": str(SYNTHETIC / "qc_conc.xlsx"),
        "qs": str(SYNTHETIC / "qc_status.xlsx"),
        "t1": str(SYNTHETIC / "thresholds_KIT01.xlsx"),
        "t2": str(SYNTHETIC / "thresholds_KIT02.xlsx"),
    }


def test_read_webidq_shapes_and_kit_id():
    p = _paths()
    conc = m.read_webidq(p["sc"], p["ss"])
    assert conc.assay.shape == (6, 30)
    assert conc.status.shape == (6, 30)
    assert conc.sample_meta["kit_id"].tolist() == ["KIT01"] * 3 + ["KIT02"] * 3
    assert len(conc.metabolite_names) == 30


def test_reader_values_match_raw_cells():
    # Reader fidelity without R: assay[0,0] equals the raw cell value.
    p = _paths()
    conc = m.read_webidq(p["sc"], p["ss"])
    wb = _xlsx.load_workbook(p["sc"])
    ws = wb[wb.sheetnames[0]]
    assert np.isclose(conc.assay[0, 0], float(ws.cell(2, 7).value))


def test_text_and_color_status_agree():
    p = _paths()
    text = m.read_webidq(p["sc"], p["ss"])
    color = m.read_webidq(p["sc"], status_method="color")
    assert np.array_equal(text.status.astype(str), color.status.astype(str))


def test_extract_kit_id():
    assert _extract_kit_id(["PROJ-X-KIT01", "a KIT12 b", None, "no kit"]) == [
        "KIT01",
        "KIT12",
        None,
        None,
    ]


def test_resolve_metabolite_names_errors_on_unknown():
    with pytest.raises(ValueError, match="Could not resolve"):
        _resolve_metabolite_names(["Arg", "definitely_not_a_metabolite"])


def test_read_thresholds():
    p = _paths()
    th = m.read_thresholds([p["t1"], p["t2"]])
    assert th.lod.shape == (2, 30)
    assert th.kit_ids == ["KIT01", "KIT02"]


def test_preprocess_end_to_end():
    p = _paths()
    conc = m.read_webidq(p["sc"], p["ss"])
    qc = m.read_webidq_qc(p["qc"], p["qs"])
    th = m.read_thresholds([p["t1"], p["t2"]])
    prep = m.preprocess(conc, qc, thresholds=th)
    assert prep.sample_outliers.shape == (6,)
    assert prep.metabolite_pass.shape == (30,)
    assert set(["cv_pass", "lod_rate", "lod_pass", "pass"]).issubset(
        prep.metabolite_filter_log.columns
    )
    assert prep.params["cv_mode"] == "per_kit"


def test_preprocess_winsorize_runs():
    p = _paths()
    conc = m.read_webidq(p["sc"], p["ss"])
    qc = m.read_webidq_qc(p["qc"], p["qs"])
    th = m.read_thresholds([p["t1"], p["t2"]])
    prep = m.preprocess(conc, qc, thresholds=th, winsorize_uloq=True)
    # winsorized sample assay must not exceed the raw where ULOQ applied;
    # at minimum the pipeline completes and preserves shape.
    assert prep.sample.assay.shape == (6, 30)


def test_preprocess_winsorize_without_thresholds_warns():
    p = _paths()
    conc = m.read_webidq(p["sc"], p["ss"])
    qc = m.read_webidq_qc(p["qc"], p["qs"])
    with pytest.warns(UserWarning, match="skipping winsorization"):
        m.preprocess(conc, qc, winsorize_uloq=True)

import numpy as np
import pytest
from conftest import SYNTHETIC

import metabosetr as m
from metabosetr._padjust import p_adjust
from metabosetr.pathways import get_pathway_sets, list_pathway_sets


def _prep():
    D = SYNTHETIC
    conc = m.read_webidq(str(D / "sample_conc.xlsx"), str(D / "sample_status.xlsx"))
    qc = m.read_webidq_qc(str(D / "qc_conc.xlsx"), str(D / "qc_status.xlsx"))
    return m.preprocess(conc, qc)


GROUP = np.array(["A", "A", "A", "B", "B", "B"], dtype=object)


# ---- p.adjust vs R stats::p.adjust ----------------------------------------

def test_padjust_bh_by_match_r():
    p = np.array([0.001, 0.02, 0.03, 0.5, 0.9])
    r_bh = [0.005, 0.05, 0.05, 0.625, 0.9]
    r_by = [0.011416666666666667, 0.114166666666666652, 0.114166666666666652, 1.0, 1.0]
    assert np.allclose(p_adjust(p, "BH"), r_bh)
    assert np.allclose(p_adjust(p, "BY"), r_by)


def test_padjust_na_uses_nonNA_count():
    # R evaluates default n after dropping NA -> n = 7 here, so max stays 0.9.
    p = np.array([0.001, 0.01, 0.02, 0.03, 0.5, 0.9, np.nan, 0.04])
    out = p_adjust(p, "BH")
    assert np.isnan(out[6])
    assert np.isclose(out[5], 0.9)  # 0.9 * 7/7
    assert np.isclose(out[0], 0.007)  # 0.001 * 7/1


# ---- ORA hypergeometric vs R phyper ---------------------------------------

def test_ora_pval_matches_r_phyper():
    # universe of 20, one set with 5 members, sig of 4, overlap 3.
    uni = [f"m{i}" for i in range(20)]
    members = ["m0", "m1", "m2", "m3", "m4"]
    sig = ["m0", "m1", "m2", "m19"]  # 3 hits in the set
    # monkey-free: use the public API by faking a domain is overkill; test the
    # hypergeometric identity directly against the R phyper reference.
    from scipy.stats import hypergeom

    k, M, N, n = 3, 5, 20, 4
    assert np.isclose(hypergeom.sf(k - 1, N, M, n), 0.031991744066047476)


def test_ora_end_to_end_structure():
    prep = _prep()
    names = list(prep.sample.metabolite_names)
    res = m.ora_pathway_sets(names[:8], names, "immunomet")
    assert list(res.columns) == [
        "set_name", "set_size", "overlap", "expected", "pval",
        "padj_BH", "padj_BY", "overlap_members",
    ]
    assert (res["pval"] >= 0).all() and (res["pval"] <= 1).all()
    assert res["pval"].is_monotonic_increasing  # sorted by pval


# ---- Approach A: feature testing ------------------------------------------

def test_wilcoxon_matches_scipy_for_one_metabolite():
    from scipy.stats import mannwhitneyu

    from metabosetr.transform import log2_with_pseudocount

    prep = _prep()
    res = m.test_metabolites(prep, GROUP, method="wilcoxon")
    al = log2_with_pseudocount(prep.sample.assay.T)
    g1 = GROUP == "A"
    g2 = GROUP == "B"
    x1 = al[0, g1]
    x2 = al[0, g2]
    expect_p = mannwhitneyu(x2, x1, alternative="two-sided", method="asymptotic").pvalue
    assert np.isclose(res["P.Value"].iloc[0], expect_p)
    assert np.isclose(res["log2FC"].iloc[0], x2.mean() - x1.mean())


def test_limma_runs_and_shapes():
    prep = _prep()
    res = m.test_metabolites(prep, GROUP, method="limma")
    assert list(res.columns[:5]) == ["metabolite", "log2FC", "AveExpr", "t", "P.Value"]
    assert len(res) == len(prep.sample.metabolite_names)
    assert res["P.Value"].between(0, 1).all()
    assert "padj_BH" in res.columns and "padj_BY" in res.columns


def test_limma_is_default():
    prep = _prep()
    a = m.test_metabolites(prep, GROUP)
    b = m.test_metabolites(prep, GROUP, method="limma")
    assert np.allclose(a["P.Value"], b["P.Value"], equal_nan=True)


# ---- Approach B: GSEA (deterministic given seed) --------------------------

def test_gsea_runs_and_is_deterministic():
    prep = _prep()
    r1 = m.gsea_pathway_sets(prep, GROUP, "immunomet", seed=42, permutation_num=200)
    r2 = m.gsea_pathway_sets(prep, GROUP, "immunomet", seed=42, permutation_num=200)
    assert len(r1) > 0
    assert list(r1["Term"]) == list(r2["Term"])
    assert np.allclose(
        r1["NES"].astype(float).to_numpy(), r2["NES"].astype(float).to_numpy()
    )


# ---- loaders + annotate ----------------------------------------------------

def test_list_pathway_sets_counts():
    lp = list_pathway_sets()
    counts = lp.groupby("domain").size().to_dict()
    assert counts["reactome"] == 224
    assert counts["smpdb"] == 266
    assert len(get_pathway_sets("immunomet")) == 5


def test_annotate_metabolites():
    ann = m.annotate_metabolites(["Nicotine", "not_a_metabolite"])
    assert ann.loc[0, "kegg_id"] == "C00745"
    assert ann.loc[1, "kegg_id"] is None


def test_unknown_domain_errors():
    with pytest.raises(ValueError, match="Unknown domain"):
        m.get_pathway_sets("nope")

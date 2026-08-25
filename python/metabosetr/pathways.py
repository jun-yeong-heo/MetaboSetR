"""Pathway-set loaders (port of ``pathway-sets.R``).

DESIGN.md §6.4 / §8.5. All sets are metabolite-level (members are
metabolite ``short_name``s). Built at import time from the committed master
TSVs -- the same source ``make_sysdata.R`` uses -- so ``pathway_sets`` /
``pathway_sets_meta`` match the R sysdata objects.
"""

from __future__ import annotations

import functools

import pandas as pd

from .reference import _DATA_RAW, EXTDATA_PATHWAY_SETS, _read_tsv

# domain -> master TSV path. Member column is "shortname" for every domain.
_MASTERS = {
    "immunomet": _DATA_RAW / "pathway_sets" / "immunomet" / "immunomet_metabolites_master.tsv",
    "reactome": _DATA_RAW / "pathway_sets" / "ramp" / "reactome_master.tsv",
    "wikipathways": _DATA_RAW / "pathway_sets" / "ramp" / "wikipathways_master.tsv",
    "smpdb": _DATA_RAW / "pathway_sets" / "ramp" / "smpdb_master.tsv",
    "lion": _DATA_RAW / "pathway_sets" / "lion" / "lion_master.tsv",
    "source": _DATA_RAW / "pathway_sets" / "ramp" / "source_master.tsv",
    "health": _DATA_RAW / "pathway_sets" / "ramp" / "health_master.tsv",
}

_DOMAINS = tuple(_MASTERS.keys())


@functools.lru_cache(maxsize=None)
def _master(domain: str) -> pd.DataFrame:
    return _read_tsv(_MASTERS[domain])


def _check_domain(domain: str) -> None:
    if domain not in _MASTERS:
        raise ValueError(
            f"Unknown domain '{domain}'. Available: {', '.join(_DOMAINS)}"
        )


@functools.lru_cache(maxsize=None)
def get_pathway_sets(domain: str) -> dict:
    """Named dict ``set_name -> [members]`` for a domain (fgsea-style input).

    Members are unique with first-appearance order preserved, matching R
    ``load_set_dict``.
    """
    _check_domain(domain)
    df = _master(domain)
    out = {}
    for sn in pd.unique(df["set_name"]):
        members = pd.unique(df.loc[df["set_name"] == sn, "shortname"])
        out[sn] = list(members)
    return out


def list_pathway_sets() -> pd.DataFrame:
    """All ``(domain, set_name, n_members)`` rows across every domain."""
    rows = []
    for dom in _DOMAINS:
        for sn, members in get_pathway_sets(dom).items():
            rows.append({"domain": dom, "set_name": sn, "n_members": len(members)})
    return pd.DataFrame(rows)


def get_pathway_sets_meta(domain: str) -> pd.DataFrame:
    """Long-format ``domain, set_name, member, direction, brief_note``.

    ``direction`` is taken from the master when present, else None; no master
    carries ``brief_note``, so it is None (matching R ``load_set_long``).
    """
    _check_domain(domain)
    df = _master(domain)
    direction = df["direction"] if "direction" in df.columns else None
    return pd.DataFrame(
        {
            "domain": domain,
            "set_name": df["set_name"].to_numpy(),
            "member": df["shortname"].to_numpy(),
            "direction": direction.to_numpy() if direction is not None else None,
            "brief_note": None,
        }
    )


def pathway_set_gmt_path(domain: str) -> str:
    """Path to the installed ``<domain>_metabolites.gmt``."""
    _check_domain(domain)
    return str(EXTDATA_PATHWAY_SETS / f"{domain}_metabolites.gmt")

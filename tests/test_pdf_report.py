import pandas as pd
from build_pdf_report import TOP5, build_html


def _df():
    return pd.DataFrame([
        {"gene": "AR", "source": "UniProt", "source_id": "P10275-3",
         "change_class": "substituted segment", "nmd_flag": "NMD-candidate",
         "domains_hit": "NR LBD", "loc_flags": "", "pca_evidence": "CRPC",
         "gtex_baseline": "multi-isoform normal", "corroborating_sources": 3},
        {"gene": "AR", "source": "Literature", "source_id": "AR-V7",
         "change_class": "LBD-truncated", "nmd_flag": "", "domains_hit": "",
         "loc_flags": "", "pca_evidence": "resistance",
         "gtex_baseline": "", "corroborating_sources": 3},
        {"gene": "X", "source": "TCGA", "source_id": "X_ES_1",
         "change_class": "truncation", "nmd_flag": "no", "domains_hit": "",
         "loc_flags": "", "pca_evidence": "", "gtex_baseline": "",
         "corroborating_sources": 1},
    ])


def test_top5_well_formed():
    assert len(TOP5) == 5
    assert [c["rank"] for c in TOP5] == [1, 2, 3, 4, 5]
    assert TOP5[0]["gene"] == "AR"
    for c in TOP5:
        for k in ["gene", "variant", "event", "impact", "pca", "why", "sources"]:
            assert c[k]


def test_build_html_has_core_sections():
    figs = {k: f"/tmp/{k}.png" for k in ("sources", "func", "corr", "arv7")}
    html = build_html(_df(), figs)
    assert "Top 5 splicing cases" in html
    assert "AR-V7" in html
    assert "Cross-source corroboration" in html
    assert "Limitations" in html
    # corroborated gene AR appears in the table, singleton X is excluded
    assert ">AR<" in html

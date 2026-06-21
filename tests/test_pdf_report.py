import pandas as pd
from build_pdf_report import (
    TOP5, NOVEL5, SEQ_REPS, build_html, collect_sequences, render_seq_html,
)


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


def test_novel5_well_formed():
    assert len(NOVEL5) == 5
    assert [c["rank"] for c in NOVEL5] == [1, 2, 3, 4, 5]
    assert {c["gene"] for c in NOVEL5} == {"MAP3K7", "EXOC7", "PRUNE2",
                                           "ENAH", "STEAP3"}
    # novelty cases must distinguish established context from the novel angle
    for c in NOVEL5:
        for k in ("event", "known", "novel", "signal", "sources"):
            assert c[k]


def test_render_seq_highlights_changed_region():
    # prefix 2 identical, suffix 1 identical -> residues 3..4 changed
    html = render_seq_html("ABCDE", prefix=2, suffix=1)
    assert "<mark>CD</mark>" in html
    assert html.startswith("AB")
    assert html.endswith("E")


def test_render_seq_windows_long_sequence():
    seq = "M" * 1000 + "QWERTY" + "K" * 1000  # changed block in the middle
    html = render_seq_html(seq, prefix=1000, suffix=1000)
    assert "<mark>QWERTY</mark>" in html
    assert "identical]" in html  # long flanks elided


def test_collect_sequences_has_before_after_for_featured():
    seqs = collect_sequences()
    assert len(seqs) == len(SEQ_REPS)
    ar = next(s for s in seqs if s["gene"] == "AR")
    assert ar["available"] and ar["before"] and ar["after"]
    assert len(ar["before"]) == 920 and len(ar["after"]) == 644
    # CCND1 is literature-only: no sequence, handled gracefully
    ccnd1 = next(s for s in seqs if s["gene"] == "CCND1")
    assert ccnd1["available"] is False


def test_build_html_has_core_sections():
    figs = {k: f"/tmp/{k}.png"
            for k in ("sources", "func", "corr", "arv7", "novel")}
    html = build_html(_df(), figs)
    assert "Top 5 splicing cases" in html
    assert "AR-V7" in html
    assert "Cross-source corroboration" in html
    assert "Novel but credible candidates" in html
    assert "MAP3K7" in html and "EXOC7" in html
    assert "Limitations" in html
    assert "Before / after protein sequences" in html
    assert "seqblock" in html and "<mark>" in html
    # corroborated gene AR appears in the table, singleton X is excluded
    assert ">AR<" in html

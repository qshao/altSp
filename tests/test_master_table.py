from integrate.master_table import (
    COLUMNS, normalize_uniprot, add_corroboration,
)


def test_columns_contains_core_schema():
    for c in ["variant_id", "gene", "uniprot", "source", "change_class",
              "domains_hit", "nmd_flag", "corroborating_sources", "provenance"]:
        assert c in COLUMNS


def test_normalize_uniprot_row_shape():
    rec = {"gene": "AR", "accession": "P10275", "canonical_len": 920,
           "canonical_seq": "M" * 920,
           "is_prostate_disease": True,
           "diseases": [{"id": "Prostate cancer, hereditary, X-linked 3",
                         "desc": ""}]}
    iso = {"isoform_id": "P10275-3", "isoform_name": "3",
           "synonyms": ["AR-V7"], "after_len": 644, "after_seq": "M" * 644,
           "change_class": "substituted segment", "varseq": []}
    row = normalize_uniprot(rec, iso)
    assert row["gene"] == "AR"
    assert row["source"] == "UniProt"
    assert row["uniprot"] == "P10275"
    assert row["after_len"] == 644
    assert "AR-V7" in row["provenance"] or "AR-V7" in row["source_id"]
    assert set(COLUMNS) <= set(row)


def test_add_corroboration_counts_distinct_sources():
    rows = [
        {"gene": "ITGA6", "source": "UniProt"},
        {"gene": "ITGA6", "source": "TCGA"},
        {"gene": "AR", "source": "UniProt"},
    ]
    out = add_corroboration(rows)
    by = {(r["gene"], r["source"]): r["corroborating_sources"] for r in out}
    assert by[("ITGA6", "UniProt")] == 2
    assert by[("AR", "UniProt")] == 1

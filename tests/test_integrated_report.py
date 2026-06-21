import pandas as pd
from build_integrated_report import corroborated


def test_corroborated_filters_and_sorts():
    df = pd.DataFrame([
        {"gene": "AR", "corroborating_sources": 3},
        {"gene": "X", "corroborating_sources": 1},
        {"gene": "ITGA6", "corroborating_sources": 2},
    ])
    out = corroborated(df)
    assert list(out["gene"]) == ["AR", "ITGA6"]

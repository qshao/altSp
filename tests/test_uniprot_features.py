import json
from annotate.uniprot_features import parse_features, Feature


def test_parse_features_extracts_typed_intervals():
    rec = json.load(open("tests/fixtures/uniprot_AR_record.json"))
    feats = parse_features(rec)
    assert feats, "expected at least one feature"
    assert all(isinstance(f, Feature) for f in feats)
    assert all(f.start <= f.end for f in feats)
    types = {f.type for f in feats}
    # AR record carries Domain and/or Region features
    assert types & {"Domain", "Region", "Zinc finger", "DNA binding"}

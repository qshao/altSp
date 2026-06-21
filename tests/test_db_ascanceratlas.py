import json
from collect.db_ascanceratlas import parse_items


def test_parse_items_keeps_only_prostate_rows():
    payload = json.load(open("tests/fixtures/ascancer_sample.json"))
    rows = parse_items(payload["items"])
    assert rows, "expected >=1 prostate row"
    assert all("gene" in r and "event_type" in r for r in rows)
    assert all("prostate" in r["cancer"].lower() for r in rows)
    # the fixture's non-prostate row must be filtered out
    assert len(rows) < len(payload["items"])

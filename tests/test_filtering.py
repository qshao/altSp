from pathlib import Path
from altsplice_protein.filtering import (
    parse_isoforms, is_significant, iter_significant_events,
    collect_unique_isoforms,
)

FIX = Path(__file__).parent / "fixtures" / "mini.csv"


def test_parse_isoforms_splits_and_drops_na():
    assert parse_isoforms("A2M-001,A2M-002") == ["A2M-001", "A2M-002"]
    assert parse_isoforms("NA") == []
    assert parse_isoforms("") == []
    assert parse_isoforms(" TST-003 ") == ["TST-003"]


def test_is_significant_thresholds():
    assert is_significant(0.001, 0.3) is True
    assert is_significant(0.001, -0.3) is True
    assert is_significant(0.5, 0.3) is False     # FDR too high
    assert is_significant(0.001, 0.05) is False  # dPSI too small
    assert is_significant(None, 0.3) is False


def test_iter_significant_events_filters_and_parses():
    events = list(iter_significant_events(FIX))
    ids = {e.splice_event for e in events}
    assert ids == {"TST_ES_1", "TST_AP_1", "NC_ES_1"}  # NS_ES_1 excluded
    by_id = {e.splice_event: e for e in events}
    assert by_id["TST_ES_1"].splice_in == ["TST-001"]
    assert by_id["TST_ES_1"].splice_out == ["TST-002"]
    assert by_id["TST_AP_1"].splice_out == []          # NA -> empty


def test_collect_unique_isoforms():
    events = list(iter_significant_events(FIX))
    assert collect_unique_isoforms(events) == {
        "TST-001", "TST-002", "TST-003", "NCG-001", "NCG-002",
    }

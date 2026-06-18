from altsplice_protein.comparison import representative, compare, compare_event


def test_representative_picks_longest():
    assert representative(["AA", "AAAA", None]) == "AAAA"
    assert representative([None, None]) is None
    assert representative([]) is None


def test_compare_identical():
    assert compare("MAAA", "MAAA") == "identical"


def test_compare_clean_length_change():
    # after shorter than before, before starts with after -> clean exon removal
    assert compare("MAAAAAAAAAKDDDD", "MAAAAAAAAAK") == "length_change(-4aa)"
    assert compare("MAAAAAAAAAK", "MAAAAAAAAAKDDDD") == "length_change(+4aa)"


def test_compare_frameshift_like():
    assert compare("MAAAKDDDD", "MAAAKEEEE") == "frameshift_or_seqchange"


def test_compare_missing_sides():
    assert compare(None, "MAAA") == "noncoding_side"
    assert compare("MAAA", None) == "noncoding_side"
    assert compare(None, None) == "unresolved"


def test_compare_event_single_isoform():
    assert compare_event([], ["MAAA"], has_both_sides=False) == "single_isoform"


def test_compare_event_uses_representatives():
    res = compare_event(["MAAAAAAAAAKDDDD"], ["MAAA", "MAAAAAAAAAK"],
                        has_both_sides=True)
    assert res == "length_change(-4aa)"

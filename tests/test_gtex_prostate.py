from collect.gtex_prostate import summarize_transcripts


def test_summarize_flags_multi_isoform_normal():
    payload = {"data": [
        {"median": 12.0, "transcriptId": "t1"},
        {"median": 4.0, "transcriptId": "t2"},
        {"median": 0.1, "transcriptId": "t3"},
    ], "gencodeId": "ENSG00000169083.16"}
    s = summarize_transcripts(payload)
    assert s["n_transcripts"] == 3
    assert s["max_median"] == 12.0
    assert s["second_median"] == 4.0
    assert s["multi_isoform_normal"] is True


def test_summarize_single_isoform():
    payload = {"data": [{"median": 9.0, "transcriptId": "t1"},
                        {"median": 0.2, "transcriptId": "t2"}],
               "gencodeId": "X"}
    s = summarize_transcripts(payload)
    assert s["multi_isoform_normal"] is False

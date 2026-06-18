import altsplice_protein.__main__ as cli


def test_main_parses_args_and_calls_run(monkeypatch, capsys):
    captured = {}

    def fake_run(csv_path, data_dir, results_dir, fdr_max, dpsi_min, use_rest):
        captured.update(dict(
            csv_path=csv_path, data_dir=data_dir, results_dir=results_dir,
            fdr_max=fdr_max, dpsi_min=dpsi_min, use_rest=use_rest,
        ))
        return {"events": 3, "protein": 4}

    monkeypatch.setattr(cli, "run", fake_run)
    stats = cli.main([
        "spliceseq_info_PRAD.csv", "--fdr-max", "0.01",
        "--dpsi-min", "0.2", "--no-rest",
    ])
    assert captured["csv_path"] == "spliceseq_info_PRAD.csv"
    assert captured["fdr_max"] == 0.01
    assert captured["dpsi_min"] == 0.2
    assert captured["use_rest"] is False
    assert stats["events"] == 3
    out = capsys.readouterr().out
    assert "events: 3" in out

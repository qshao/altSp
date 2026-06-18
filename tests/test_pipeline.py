from pathlib import Path
from altsplice_protein.pipeline import run_with_files

FIX = Path(__file__).parent / "fixtures"


def test_run_with_files_end_to_end(tmp_path):
    stats = run_with_files(
        csv_path=FIX / "mini.csv",
        gtf_path=FIX / "mini.gtf",
        pep_path=FIX / "mini.pep.fa",
        results_dir=tmp_path,
        use_rest=False,
    )
    assert stats["events"] == 3
    assert stats["unique_isoforms"] == 5
    assert stats["protein"] == 4     # TST-001/002/003, NCG-001
    assert stats["noncoding"] == 1   # NCG-002

    fasta = (tmp_path / "proteins.fasta").read_text()
    assert ">TST-001|ENST001|ENSP001|TST" in fasta
    assert "MAAAAAAAAAK" in fasta
    assert "NCG-002" not in fasta    # non-coding excluded from FASTA

    tsv = (tmp_path / "events_proteins.tsv").read_text().splitlines()
    header = tsv[0].split("\t")
    rows = {r.split("\t")[0]: dict(zip(header, r.split("\t"))) for r in tsv[1:]}

    assert rows["TST_ES_1"]["comparison"] == "length_change(-4aa)"
    assert rows["TST_ES_1"]["has_complete_pair"] == "yes"
    assert rows["TST_ES_1"]["before_len"] == "15"
    assert rows["TST_ES_1"]["after_len"] == "11"

    assert rows["TST_AP_1"]["comparison"] == "single_isoform"
    assert rows["TST_AP_1"]["has_complete_pair"] == "no"

    assert rows["NC_ES_1"]["comparison"] == "noncoding_side"

    unresolved = (tmp_path / "unresolved.txt").read_text()
    assert "NCG-002" in unresolved   # reported as noncoding

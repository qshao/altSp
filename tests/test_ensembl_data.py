from pathlib import Path
from altsplice_protein.ensembl_data import (
    parse_gtf_transcripts, parse_pep_fasta, build_resolver_map,
)

FIX = Path(__file__).parent / "fixtures"
GTF = FIX / "mini.gtf"
PEP = FIX / "mini.pep.fa"


def test_parse_gtf_transcripts_only_transcript_rows():
    tx = parse_gtf_transcripts(GTF)
    assert set(tx) == {"TST-001", "TST-002", "TST-003", "NCG-001", "NCG-002"}
    assert tx["TST-001"] == ("ENST001", "TST", "protein_coding")
    assert tx["NCG-002"][2] == "retained_intron"


def test_parse_pep_fasta_keyed_by_transcript():
    pep = parse_pep_fasta(PEP)
    assert pep["ENST001"] == ("ENSP001", "MAAAAAAAAAK")
    assert pep["ENST002"][1] == "MAAAAAAAAAKDDDD"
    assert "ENST011" not in pep  # non-coding, no peptide


def test_build_resolver_map_status():
    m = build_resolver_map(GTF, PEP)
    assert m["TST-001"].status == "protein"
    assert m["TST-001"].protein_seq == "MAAAAAAAAAK"
    assert m["TST-001"].protein_id == "ENSP001"
    assert m["NCG-002"].status == "noncoding"
    assert m["NCG-002"].protein_seq is None
    assert m["TST-001"].source == "local"

from collect.literature_variants import SEED, parse_europepmc


def test_seed_has_core_variants_and_fields():
    genes = {s["gene"] for s in SEED}
    assert {"AR", "KLF6"} <= genes  # AR-V7, KLF6-SV1 present
    for s in SEED:
        for k in ("gene", "uniprot", "variant", "event_type", "effect", "query"):
            assert k in s and s[k]


def test_parse_europepmc_extracts_pmids():
    payload = {"resultList": {"result": [
        {"id": "111", "pmid": "111"}, {"id": "222", "pmid": "222"}]}}
    assert parse_europepmc(payload) == ["111", "222"]

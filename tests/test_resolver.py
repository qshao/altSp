import altsplice_protein.resolver as resolver
from altsplice_protein.models import IsoformProtein


def _map():
    return {
        "TST-001": IsoformProtein(
            name="TST-001", transcript_id="ENST001", protein_id="ENSP001",
            gene_symbol="TST", biotype="protein_coding",
            protein_seq="MAAA", source="local", status="protein",
        )
    }


def test_resolve_all_local_hit():
    out = resolve_all_only_local()
    assert out["TST-001"].source == "local"
    assert out["TST-001"].status == "protein"


def resolve_all_only_local():
    return resolver.resolve_all(["TST-001"], _map(), use_rest=False)


def test_resolve_all_miss_without_rest_is_unresolved():
    out = resolver.resolve_all(["GHOST-001"], _map(), use_rest=False)
    assert out["GHOST-001"].status == "unresolved"
    assert out["GHOST-001"].transcript_id is None


def test_resolve_all_miss_uses_rest(monkeypatch):
    def fake_rest(name):
        return IsoformProtein(
            name=name, transcript_id="ENST999", protein_id="ENSP999",
            protein_seq="MKKK", source="rest", status="protein",
        )
    monkeypatch.setattr(resolver, "resolve_rest", fake_rest)
    out = resolver.resolve_all(["GHOST-001"], _map(), use_rest=True)
    assert out["GHOST-001"].source == "rest"
    assert out["GHOST-001"].protein_seq == "MKKK"


def test_resolve_rest_parses_transcript_and_protein(monkeypatch):
    calls = {}

    def fake_get(path):
        if path.startswith("/xrefs/symbol/"):
            calls["xref"] = path
            return [
                {"type": "gene", "id": "ENSG9"},
                {"type": "transcript", "id": "ENST999"},
            ]
        if path.startswith("/sequence/id/"):
            calls["seq"] = path
            return {"id": "ENSP999", "seq": "MKKK"}
        raise AssertionError(path)

    monkeypatch.setattr(resolver, "_rest_get", fake_get)
    iso = resolver.resolve_rest("A2M-001")
    assert iso.transcript_id == "ENST999"
    assert iso.protein_seq == "MKKK"
    assert iso.status == "protein"
    assert iso.source == "rest"
    assert "ENST999" in calls["seq"]

from annotate.uniprot_features import Feature
from annotate.functional_impact import (
    changed_interval, features_in, nmd_flag, annotate_variant,
)


def test_changed_interval_basic():
    # before/after share prefix "ABC" and suffix "Z"; middle differs
    assert changed_interval("ABCDEFZ", "ABCXZ") == (4, 6)


def test_changed_interval_truncation():
    # after is a strict prefix of before -> changed region is the lost tail
    assert changed_interval("ABCDEFG", "ABCD") == (5, 7)


def test_changed_interval_identical_is_none():
    assert changed_interval("ABCDEF", "ABCDEF") is None


def test_features_in_overlap():
    feats = [Feature("Domain", 1, 5, "Kinase domain"),
             Feature("Domain", 50, 90, "SH2 domain")]
    hit = features_in((4, 6), feats, {"Domain"})
    assert hit == ["Kinase domain"]


def test_nmd_flag_truncation_not_in_last_exon():
    # large C-terminal loss, no shared suffix -> candidate
    assert nmd_flag(900, 600, "substituted segment", 0) == "NMD-candidate"


def test_nmd_flag_small_tail_change_is_no():
    assert nmd_flag(900, 895, "C-terminal truncation", 0) == "no"


def test_annotate_variant_assembles_keys():
    feats = [Feature("Domain", 4, 8, "Ligand-binding"),
             Feature("Transmembrane", 4, 6, "Helical"),
             Feature("Region", 4, 6, "Disordered")]
    out = annotate_variant("ABCDEFGHZ", "ABCXZ", "substituted segment", feats)
    assert out["changed_interval"] == (4, 8)
    assert "Ligand-binding" in out["domains_hit"]
    assert any("Transmembrane" in f for f in out["loc_flags"])
    assert out["disorder_overlap"] == "yes"
    assert out["nmd_flag"] in {"NMD-candidate", "no", "n.a."}

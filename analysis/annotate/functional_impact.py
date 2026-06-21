"""Predict how a splice change affects protein function, from UniProt features."""
from __future__ import annotations
from annotate.uniprot_features import DOMAIN_TYPES, SITE_TYPES, LOC_TYPES

NMD_TAIL_AA = 50  # ~50 nt rule, approximated at protein level


def _common_prefix(a: str, b: str) -> int:
    n = min(len(a), len(b)); i = 0
    while i < n and a[i] == b[i]:
        i += 1
    return i


def _common_suffix(a: str, b: str, used: int) -> int:
    n = min(len(a), len(b)) - used; i = 0
    while i < n and a[-1 - i] == b[-1 - i]:
        i += 1
    return i


def changed_interval(before: str, after: str):
    if not before or not after or before == after:
        return None
    pre = _common_prefix(before, after)
    suf = _common_suffix(before, after, pre)
    start = pre + 1
    end = len(before) - suf
    if end < start:
        end = start
    return (start, end)


def features_in(interval, feats, types) -> list[str]:
    if interval is None:
        return []
    s, e = interval
    out = []
    for f in feats:
        if f.type in types and not (f.end < s or f.start > e):
            out.append(f.description or f.type)
    return out


def _loc_flags(interval, feats) -> list[str]:
    if interval is None:
        return []
    s, e = interval
    out = []
    for f in feats:
        if f.type in LOC_TYPES and not (f.end < s or f.start > e):
            out.append(f"{f.type}:{f.description}" if f.description else f.type)
    return out


def _disorder(interval, feats) -> str:
    if interval is None:
        return "n.a."
    s, e = interval
    for f in feats:
        if f.type == "Region" and "disorder" in (f.description or "").lower():
            if not (f.end < s or f.start > e):
                return "yes"
    return "no"


def nmd_flag(before_len, after_len, change_class, common_suffix) -> str:
    if not before_len or not after_len:
        return "n.a."
    lost_cterm = before_len - after_len
    # truncation/substitution that removes the C-terminus and is not confined to
    # the last ~50 aa is an NMD candidate (protein-level heuristic).
    if common_suffix == 0 and lost_cterm > NMD_TAIL_AA:
        return "NMD-candidate"
    return "no"


def annotate_variant(before, after, change_class, feats) -> dict:
    iv = changed_interval(before or "", after or "")
    pre = _common_prefix(before or "", after or "") if before and after else 0
    suf = _common_suffix(before or "", after or "", pre) if before and after else 0
    return {
        "changed_interval": iv,
        "domains_hit": features_in(iv, feats, DOMAIN_TYPES),
        "regions_lost": features_in(iv, feats, SITE_TYPES),
        "loc_flags": _loc_flags(iv, feats),
        "disorder_overlap": _disorder(iv, feats),
        "nmd_flag": nmd_flag(len(before) if before else 0,
                             len(after) if after else 0, change_class, suf),
    }

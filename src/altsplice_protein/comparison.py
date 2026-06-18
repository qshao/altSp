from __future__ import annotations
from typing import Iterable


def representative(seqs: Iterable[str | None]) -> str | None:
    present = [s for s in seqs if s]
    return max(present, key=len) if present else None


def compare(before: str | None, after: str | None) -> str:
    if before is None and after is None:
        return "unresolved"
    if before is None or after is None:
        return "noncoding_side"
    if before == after:
        return "identical"
    short, long = sorted([before, after], key=len)
    if long.startswith(short):
        delta = len(after) - len(before)
        return f"length_change({delta:+d}aa)"
    return "frameshift_or_seqchange"


def compare_event(
    before_seqs: list, after_seqs: list, has_both_sides: bool
) -> str:
    if not has_both_sides:
        return "single_isoform"
    return compare(representative(before_seqs), representative(after_seqs))

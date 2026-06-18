from __future__ import annotations
import csv
from typing import Iterable, Iterator
from .models import Event


def parse_isoforms(field: str) -> list[str]:
    out: list[str] = []
    for tok in (field or "").split(","):
        tok = tok.strip()
        if tok and tok != "NA":
            out.append(tok)
    return out


def _to_float(x) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def is_significant(fdr, dpsi, fdr_max: float = 0.05, dpsi_min: float = 0.1) -> bool:
    return (
        fdr is not None and dpsi is not None
        and fdr < fdr_max and abs(dpsi) >= dpsi_min
    )


def iter_significant_events(
    csv_path, fdr_max: float = 0.05, dpsi_min: float = 0.1
) -> Iterator[Event]:
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            fdr = _to_float(row["FDR_Difference"])
            dpsi = _to_float(row["PSI_Difference"])
            if not is_significant(fdr, dpsi, fdr_max, dpsi_min):
                continue
            yield Event(
                splice_event=row["Splice_Event"],
                gene_symbol=row["Gene_Symbol"],
                splice_type=row["Splice_Type"],
                psi_difference=dpsi,
                fdr_difference=fdr,
                splice_in=parse_isoforms(row["SpliceIn_IsoName"]),
                splice_out=parse_isoforms(row["SpliceOut_IsoName"]),
            )


def collect_unique_isoforms(events: Iterable[Event]) -> set[str]:
    names: set[str] = set()
    for e in events:
        names.update(e.splice_in)
        names.update(e.splice_out)
    return names

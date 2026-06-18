from __future__ import annotations
from dataclasses import dataclass


@dataclass
class Event:
    splice_event: str
    gene_symbol: str
    splice_type: str
    psi_difference: float
    fdr_difference: float
    splice_in: list[str]   # inclusion isoform names ("after")
    splice_out: list[str]  # exclusion isoform names ("before")


@dataclass
class IsoformProtein:
    name: str
    transcript_id: str | None = None
    protein_id: str | None = None
    gene_symbol: str | None = None
    biotype: str | None = None
    protein_seq: str | None = None
    source: str = "local"       # 'local' | 'rest'
    status: str = "unresolved"  # 'protein' | 'noncoding' | 'unresolved'

"""Extract before/after protein sequences for prioritized splicing events and
compute a structural-impact summary used by the Word report.

"before" = representative splice-OUT protein (longest resolved isoform)
"after"  = representative splice-IN  protein
(matching the pipeline's convention in pipeline.py).

Impact metrics are sequence-derived (no external DB): how much of the protein is
preserved before the splice diverges, and whether the C-terminus is
lost/altered/extended. These feed the "predicted functional impact" reasoning.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import json
import pandas as pd

FASTA = "results/proteins.fasta"
EVENTS = "results/events_proteins.tsv"
PRIORITIZED = "analysis/prioritized_events.tsv"
OUT = "analysis/sequence_impact.json"


def load_fasta(path: str) -> dict[str, tuple[str, str]]:
    """protein_id -> (isoform_name, sequence)."""
    by_pid: dict[str, tuple[str, str]] = {}
    name = pid = None
    buf: list[str] = []

    def flush():
        if pid and buf:
            by_pid[pid] = (name, "".join(buf))

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                flush()
                parts = line[1:].strip().split("|")
                name = parts[0]
                pid = parts[2] if len(parts) > 2 else parts[0]
                buf = []
            else:
                buf.append(line.strip())
        flush()
    return by_pid


def representative(pid_field: str, fasta: dict) -> tuple[str | None, str | None]:
    """Longest resolved protein among a ';'-joined protein-id field."""
    best = None
    for pid in (pid_field or "").split(";"):
        pid = pid.strip()
        if pid and pid != "." and pid in fasta:
            name, seq = fasta[pid]
            if best is None or len(seq) > len(best[2]):
                best = (name, pid, seq)
    if best is None:
        return None, None
    return best[1], best[2]


def common_prefix(a: str, b: str) -> int:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return n


def common_suffix(a: str, b: str, cap: int) -> int:
    n = 0
    for x, y in zip(reversed(a), reversed(b)):
        if x != y or n >= cap:
            break
        n += 1
    return n


@dataclass
class Impact:
    splice_event: str
    gene: str
    splice_type: str
    endpoint: str
    HR: float
    direction: str
    bh_p: float
    pca_gene: bool
    comparison: str
    before_pid: str | None
    after_pid: str | None
    before_len: int
    after_len: int
    identical_prefix: int
    identical_suffix: int
    pct_preserved: float
    net_aa_change: int
    impact_class: str
    impact_detail: str
    before_seq: str
    after_seq: str


def classify(before: str, after: str, comparison: str) -> tuple[str, str, int, int, float]:
    pre = common_prefix(before, after)
    suf = common_suffix(before, after, min(len(before), len(after)) - pre)
    net = len(after) - len(before)
    pct = round(100 * pre / len(before), 1) if before else 0.0
    if comparison.startswith("length_change"):
        if suf > 0 and pre + suf >= min(len(before), len(after)):
            cls = "in-frame indel"
            det = (f"In-frame change: identical for first {pre} aa and last {suf} aa; "
                   f"net {net:+d} aa internally/terminally. Domain architecture "
                   f"largely retained unless the altered segment overlaps a functional region.")
        else:
            cls = "in-frame terminal change"
            det = (f"Identical for first {pre} aa; C-terminus altered ({net:+d} aa). "
                   f"Likely modifies a C-terminal domain/motif.")
    elif comparison == "frameshift_or_seqchange":
        lost = len(before) - pre
        if net <= -0.3 * len(before):
            cls = "truncation"
            det = (f"Sequence preserved for {pre} aa ({pct}% of the original), then "
                   f"diverges; product is {abs(net)} aa shorter. The C-terminal "
                   f"~{lost} aa of the original protein are lost or replaced — domains "
                   f"in that region are predicted non-functional.")
        else:
            cls = "frameshift / sequence swap"
            det = (f"Identical for {pre} aa ({pct}%), then the reading frame/sequence "
                   f"changes for the remaining ~{lost} aa. Any motif or domain past "
                   f"residue {pre} is altered.")
    elif comparison == "noncoding_side":
        cls = "isoform loss (one side non-coding)"
        det = ("One side of the event is non-coding / unresolved, i.e. the splice "
               "switches between a coding and a non-protein-coding product — "
               "effectively gain or loss of the entire protein.")
    else:
        cls = comparison
        det = "See sequences."
    return cls, det, pre, suf, pct


def main() -> None:
    fasta = load_fasta(FASTA)
    ev = pd.read_csv(EVENTS, sep="\t").set_index("Splice_Event")
    pri = pd.read_csv(PRIORITIZED, sep="\t")

    out: list[dict] = []
    for _, r in pri.iterrows():
        se = r["Splice_Event"]
        if se not in ev.index:
            continue
        e = ev.loc[se]
        bpid, bseq = representative(str(e["SpliceOut_protein_ids"]), fasta)
        apid, aseq = representative(str(e["SpliceIn_protein_ids"]), fasta)
        if not bseq or not aseq:
            continue  # need both sides to show before/after
        cls, det, pre, suf, pct = classify(bseq, aseq, str(e["comparison"]))
        out.append(asdict(Impact(
            splice_event=se, gene=r["Gene"], splice_type=r["Type"],
            endpoint=r["endpoint"], HR=float(r["HR"]), direction=r["direction"],
            bh_p=float(r["bh_p"]), pca_gene=bool(r["pca_gene"]),
            comparison=str(e["comparison"]), before_pid=bpid, after_pid=apid,
            before_len=len(bseq), after_len=len(aseq),
            identical_prefix=pre, identical_suffix=suf, pct_preserved=pct,
            net_aa_change=len(aseq) - len(bseq), impact_class=cls,
            impact_detail=det, before_seq=bseq, after_seq=aseq,
        )))

    with open(OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"events with both before+after protein resolved: {len(out)}")
    print("impact_class distribution:")
    print(pd.Series([o["impact_class"] for o in out]).value_counts().to_string())
    print("\nknown PCa genes with full sequences:",
          [o["gene"] for o in out if o["pca_gene"]])


if __name__ == "__main__":
    main()

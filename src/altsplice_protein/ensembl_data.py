from __future__ import annotations
import gzip
import os
import re
import urllib.request
from .models import IsoformProtein

GTF_URL = (
    "https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/"
    "Homo_sapiens.GRCh37.75.gtf.gz"
)
PEP_URL = (
    "https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/"
    "Homo_sapiens.GRCh37.75.pep.all.fa.gz"
)

_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def _open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def download(url: str, dest: str) -> str:
    if os.path.exists(dest):
        return dest
    tmp = dest + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, dest)
    return dest


def parse_gtf_transcripts(gtf_path):
    """transcript_name -> (transcript_id, gene_symbol, biotype)."""
    result: dict[str, tuple[str, str | None, str | None]] = {}
    with _open_maybe_gzip(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = dict(_ATTR_RE.findall(cols[8]))
            name = attrs.get("transcript_name")
            tid = attrs.get("transcript_id")
            if not name or not tid:
                continue
            result[name] = (
                tid, attrs.get("gene_name"), attrs.get("transcript_biotype")
            )
    return result


def parse_pep_fasta(pep_path):
    """transcript_id -> (protein_id, sequence)."""
    result: dict[str, tuple[str, str]] = {}
    pid: str | None = None
    tid: str | None = None
    seq: list[str] = []

    def flush():
        if tid and pid:
            result[tid] = (pid, "".join(seq))

    with _open_maybe_gzip(pep_path) as f:
        for line in f:
            if line.startswith(">"):
                flush()
                seq = []
                header = line[1:].split()
                pid = header[0]
                tid = None
                for tok in header[1:]:
                    if tok.startswith("transcript:"):
                        tid = tok.split(":", 1)[1]
            else:
                seq.append(line.strip())
        flush()
    return result


def build_resolver_map(gtf_path, pep_path):
    """transcript_name -> IsoformProtein (local source)."""
    tx = parse_gtf_transcripts(gtf_path)
    pep = parse_pep_fasta(pep_path)
    out: dict[str, IsoformProtein] = {}
    for name, (tid, gene, biotype) in tx.items():
        protein_id, protein_seq = (None, None)
        if tid in pep:
            protein_id, protein_seq = pep[tid]
        out[name] = IsoformProtein(
            name=name,
            transcript_id=tid,
            protein_id=protein_id,
            gene_symbol=gene,
            biotype=biotype,
            protein_seq=protein_seq,
            source="local",
            status="protein" if protein_seq else "noncoding",
        )
    return out

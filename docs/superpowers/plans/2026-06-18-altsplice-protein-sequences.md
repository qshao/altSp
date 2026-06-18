# Alternative Splicing → Protein Sequences (PRAD) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** From `spliceseq_info_PRAD.csv`, produce the "before" (spliced-out) and "after" (spliced-in) protein sequences for each significant alternative splicing event.

**Architecture:** A small Python package reads the SpliceSeq CSV, filters to significant events, extracts SpliceSeq isoform names, resolves each name to its Ensembl GRCh37 release-75 annotated protein (local GTF + pep FASTA, with REST fallback), then writes a protein FASTA, a per-event before/after TSV, and an unresolved-names report.

**Tech Stack:** Python 3 standard library (`csv`, `gzip`, `re`, `urllib`, `argparse`, `dataclasses`), pytest. Biopython is available but not required. Data: Ensembl GRCh37 release-75 GTF + `pep.all` FASTA.

## Global Constraints

- Python 3.10+ (uses `X | None` type syntax). Confirm with `python3 --version`.
- Significance filter (defaults): `FDR_Difference < 0.05` AND `abs(PSI_Difference) >= 0.1`.
- "after" = inclusion form = `SpliceIn_IsoName`; "before" = exclusion form = `SpliceOut_IsoName`.
- Annotation build is **Ensembl GRCh37 release-75** (matches SpliceSeq `GENE-001` naming). Do not substitute a GRCh38 build.
- GTF URL: `https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz`
- PEP URL: `https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz`
- REST fallback base: `https://grch37.rest.ensembl.org`
- All package modules live in `src/altsplice_protein/`; tests in `tests/`.
- Comparison categories: `identical`, `length_change(±Naa)`, `frameshift_or_seqchange`, `noncoding_side`, `single_isoform`, `unresolved`.

---

### Task 1: Project setup, data models, and event filtering

**Files:**
- Create: `.gitignore`
- Create: `pytest.ini`
- Create: `src/altsplice_protein/__init__.py` (empty)
- Create: `src/altsplice_protein/models.py`
- Create: `src/altsplice_protein/filtering.py`
- Create: `tests/__init__.py` (empty)
- Create: `tests/conftest.py`
- Create: `tests/fixtures/mini.csv`
- Test: `tests/test_filtering.py`

**Interfaces:**
- Produces: `Event` dataclass with fields `splice_event:str, gene_symbol:str, splice_type:str, psi_difference:float, fdr_difference:float, splice_in:list[str], splice_out:list[str]`.
- Produces: `IsoformProtein` dataclass with fields `name:str, transcript_id:str|None, protein_id:str|None, gene_symbol:str|None, biotype:str|None, protein_seq:str|None, source:str="local", status:str="unresolved"`.
- Produces: `parse_isoforms(field:str)->list[str]`, `is_significant(fdr, dpsi, fdr_max=0.05, dpsi_min=0.1)->bool`, `iter_significant_events(csv_path, fdr_max=0.05, dpsi_min=0.1)->Iterator[Event]`, `collect_unique_isoforms(events:Iterable[Event])->set[str]`.

- [ ] **Step 1: Initialize git and create ignore/config files**

```bash
cd /home/qshao/altSp
git init -q
```

Create `.gitignore`:

```gitignore
data/
results/
__pycache__/
*.pyc
.pytest_cache/
```

Create `pytest.ini`:

```ini
[pytest]
pythonpath = src
testpaths = tests
```

- [ ] **Step 2: Create the data models**

Create `src/altsplice_protein/__init__.py` as an empty file, then create `src/altsplice_protein/models.py`:

```python
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
```

- [ ] **Step 3: Create the test fixture CSV**

Create `tests/fixtures/mini.csv` (only the columns the code reads; covers significant ES with both sides, a non-significant row, a single-isoform AP row, and a noncoding-side ES row):

```csv
Splice_Event,Gene_Symbol,Gene_ID,Splice_Type,SpliceIn_IsoName,SpliceOut_IsoName,PSI_Difference,FDR_Difference
TST_ES_1,TST,ENSG1,ES,"TST-001","TST-002",0.3,0.001
NS_ES_1,TST,ENSG1,ES,"TST-001","TST-002",0.05,0.5
TST_AP_1,TST,ENSG1,AP,"TST-003","NA",0.2,0.001
NC_ES_1,NCG,ENSG2,ES,"NCG-001","NCG-002",-0.4,0.002
```

- [ ] **Step 4: Create `tests/__init__.py` (empty) and write the failing test**

Create `tests/test_filtering.py`:

```python
from pathlib import Path
from altsplice_protein.filtering import (
    parse_isoforms, is_significant, iter_significant_events,
    collect_unique_isoforms,
)

FIX = Path(__file__).parent / "fixtures" / "mini.csv"


def test_parse_isoforms_splits_and_drops_na():
    assert parse_isoforms("A2M-001,A2M-002") == ["A2M-001", "A2M-002"]
    assert parse_isoforms("NA") == []
    assert parse_isoforms("") == []
    assert parse_isoforms(" TST-003 ") == ["TST-003"]


def test_is_significant_thresholds():
    assert is_significant(0.001, 0.3) is True
    assert is_significant(0.001, -0.3) is True
    assert is_significant(0.5, 0.3) is False     # FDR too high
    assert is_significant(0.001, 0.05) is False  # dPSI too small
    assert is_significant(None, 0.3) is False


def test_iter_significant_events_filters_and_parses():
    events = list(iter_significant_events(FIX))
    ids = {e.splice_event for e in events}
    assert ids == {"TST_ES_1", "TST_AP_1", "NC_ES_1"}  # NS_ES_1 excluded
    by_id = {e.splice_event: e for e in events}
    assert by_id["TST_ES_1"].splice_in == ["TST-001"]
    assert by_id["TST_ES_1"].splice_out == ["TST-002"]
    assert by_id["TST_AP_1"].splice_out == []          # NA -> empty


def test_collect_unique_isoforms():
    events = list(iter_significant_events(FIX))
    assert collect_unique_isoforms(events) == {
        "TST-001", "TST-002", "TST-003", "NCG-001", "NCG-002",
    }
```

- [ ] **Step 5: Run the test to verify it fails**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_filtering.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'altsplice_protein.filtering'`.

- [ ] **Step 6: Implement `filtering.py`**

Create `src/altsplice_protein/filtering.py`:

```python
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
```

Create `tests/conftest.py` (empty placeholder for shared fixtures added later):

```python
# Shared pytest fixtures (added as needed by later tasks).
```

- [ ] **Step 7: Run the test to verify it passes**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_filtering.py -v`
Expected: PASS (4 passed).

- [ ] **Step 8: Commit**

```bash
cd /home/qshao/altSp
git add .gitignore pytest.ini src/altsplice_protein tests
git commit -m "feat: significance filtering and data models"
```

---

### Task 2: Ensembl GTF/pep parsing and resolver map

**Files:**
- Create: `src/altsplice_protein/ensembl_data.py`
- Create: `tests/fixtures/mini.gtf`
- Create: `tests/fixtures/mini.pep.fa`
- Test: `tests/test_ensembl_data.py`

**Interfaces:**
- Consumes: `IsoformProtein` from `models.py`.
- Produces module constants `GTF_URL:str`, `PEP_URL:str`.
- Produces: `download(url, dest)->str`, `parse_gtf_transcripts(gtf_path)->dict[str, tuple[str, str|None, str|None]]` mapping `transcript_name -> (transcript_id, gene_symbol, biotype)`, `parse_pep_fasta(pep_path)->dict[str, tuple[str, str]]` mapping `transcript_id -> (protein_id, sequence)`, `build_resolver_map(gtf_path, pep_path)->dict[str, IsoformProtein]` keyed by `transcript_name`.

- [ ] **Step 1: Create GTF and pep FASTA fixtures**

Create `tests/fixtures/mini.gtf` (tab-separated; includes a non-transcript line and a non-coding transcript to exercise filtering). IMPORTANT: columns are separated by literal TAB characters.

```text
#!genome-build GRCh37
1	ensembl	gene	1	100	.	+	.	gene_id "ENSG1"; gene_name "TST"; gene_biotype "protein_coding";
1	ensembl	transcript	1	100	.	+	.	gene_id "ENSG1"; transcript_id "ENST001"; gene_name "TST"; transcript_name "TST-001"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1	ensembl	transcript	1	90	.	+	.	gene_id "ENSG1"; transcript_id "ENST002"; gene_name "TST"; transcript_name "TST-002"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
1	ensembl	transcript	1	80	.	+	.	gene_id "ENSG1"; transcript_id "ENST003"; gene_name "TST"; transcript_name "TST-003"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
2	ensembl	transcript	1	100	.	+	.	gene_id "ENSG2"; transcript_id "ENST010"; gene_name "NCG"; transcript_name "NCG-001"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
2	ensembl	transcript	1	100	.	+	.	gene_id "ENSG2"; transcript_id "ENST011"; gene_name "NCG"; transcript_name "NCG-002"; gene_biotype "protein_coding"; transcript_biotype "retained_intron";
```

Verify the fixture really uses tabs:

```bash
cd /home/qshao/altSp && grep -cP "\t" tests/fixtures/mini.gtf
```
Expected: `6` (every non-comment line and the comment line counts only if it has a tab — expect at least 5). If it prints `0`, the editor inserted spaces — redo with real tabs.

Create `tests/fixtures/mini.pep.fa` (note: NCG-002 / ENST011 has no entry → non-coding):

```text
>ENSP001 pep:known chromosome:GRCh37:1:1:100:1 gene:ENSG1 transcript:ENST001 gene_biotype:protein_coding transcript_biotype:protein_coding
MAAAAAAAAAK
>ENSP002 pep:known chromosome:GRCh37:1:1:90:1 gene:ENSG1 transcript:ENST002 gene_biotype:protein_coding transcript_biotype:protein_coding
MAAAAAAAAAKDDDD
>ENSP003 pep:known chromosome:GRCh37:1:1:80:1 gene:ENSG1 transcript:ENST003 gene_biotype:protein_coding transcript_biotype:protein_coding
MAAAAAAAAAK
>ENSP010 pep:known chromosome:GRCh37:2:1:100:1 gene:ENSG2 transcript:ENST010 gene_biotype:protein_coding transcript_biotype:protein_coding
MNNNNNNNNNN
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_ensembl_data.py`:

```python
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
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_ensembl_data.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'altsplice_protein.ensembl_data'`.

- [ ] **Step 4: Implement `ensembl_data.py`**

Create `src/altsplice_protein/ensembl_data.py`:

```python
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
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_ensembl_data.py -v`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
cd /home/qshao/altSp
git add src/altsplice_protein/ensembl_data.py tests/test_ensembl_data.py tests/fixtures/mini.gtf tests/fixtures/mini.pep.fa
git commit -m "feat: Ensembl GRCh37 GTF/pep parsing and resolver map"
```

---

### Task 3: Name resolution with REST fallback

**Files:**
- Create: `src/altsplice_protein/resolver.py`
- Test: `tests/test_resolver.py`

**Interfaces:**
- Consumes: `IsoformProtein` from `models.py`; resolver map (dict `name->IsoformProtein`) from `ensembl_data.build_resolver_map`.
- Produces: `resolve_rest(name)->IsoformProtein`, `resolve_all(names:Iterable[str], resolver_map:dict, use_rest:bool=True, sleep:float=0.0)->dict[str, IsoformProtein]`. `resolve_all` returns one `IsoformProtein` per input name; local hits use the map verbatim, misses go to REST (if `use_rest`) else become `status="unresolved"`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_resolver.py` (REST is exercised via monkeypatch so the test is offline and deterministic):

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_resolver.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'altsplice_protein.resolver'`.

- [ ] **Step 3: Implement `resolver.py`**

Create `src/altsplice_protein/resolver.py`:

```python
from __future__ import annotations
import json
import time
import urllib.error
import urllib.request
from typing import Iterable
from .models import IsoformProtein

REST_BASE = "https://grch37.rest.ensembl.org"


def _rest_get(path: str):
    req = urllib.request.Request(
        REST_BASE + path, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)


def resolve_rest(name: str) -> IsoformProtein:
    iso = IsoformProtein(name=name, source="rest", status="unresolved")
    try:
        xrefs = _rest_get(
            f"/xrefs/symbol/homo_sapiens/{name}?content-type=application/json"
        )
    except urllib.error.HTTPError:
        return iso
    tid = next(
        (x["id"] for x in xrefs if x.get("type") == "transcript"), None
    )
    if not tid:
        return iso
    iso.transcript_id = tid
    try:
        seq = _rest_get(
            f"/sequence/id/{tid}?type=protein;content-type=application/json"
        )
        iso.protein_seq = seq.get("seq")
        iso.protein_id = seq.get("id")
        iso.status = "protein" if iso.protein_seq else "noncoding"
    except urllib.error.HTTPError:
        iso.status = "noncoding"  # transcript exists but has no protein
    return iso


def resolve_all(
    names: Iterable[str],
    resolver_map: dict,
    use_rest: bool = True,
    sleep: float = 0.0,
) -> dict[str, IsoformProtein]:
    resolved: dict[str, IsoformProtein] = {}
    for n in sorted(names):
        if n in resolver_map:
            resolved[n] = resolver_map[n]
        elif use_rest:
            resolved[n] = resolve_rest(n)
            if sleep:
                time.sleep(sleep)
        else:
            resolved[n] = IsoformProtein(name=n, status="unresolved")
    return resolved
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_resolver.py -v`
Expected: PASS (4 passed).

- [ ] **Step 5: Commit**

```bash
cd /home/qshao/altSp
git add src/altsplice_protein/resolver.py tests/test_resolver.py
git commit -m "feat: name resolution with Ensembl REST fallback"
```

---

### Task 4: Before/after protein comparison

**Files:**
- Create: `src/altsplice_protein/comparison.py`
- Test: `tests/test_comparison.py`

**Interfaces:**
- Produces: `representative(seqs:Iterable[str|None])->str|None` (longest non-empty sequence), `compare(before:str|None, after:str|None)->str`, `compare_event(before_seqs:list, after_seqs:list, has_both_sides:bool)->str`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_comparison.py`:

```python
from altsplice_protein.comparison import representative, compare, compare_event


def test_representative_picks_longest():
    assert representative(["AA", "AAAA", None]) == "AAAA"
    assert representative([None, None]) is None
    assert representative([]) is None


def test_compare_identical():
    assert compare("MAAA", "MAAA") == "identical"


def test_compare_clean_length_change():
    # after shorter than before, before starts with after -> clean exon removal
    assert compare("MAAAAAAAAAKDDDD", "MAAAAAAAAAK") == "length_change(-4aa)"
    assert compare("MAAAAAAAAAK", "MAAAAAAAAAKDDDD") == "length_change(+4aa)"


def test_compare_frameshift_like():
    assert compare("MAAAKDDDD", "MAAAKEEEE") == "frameshift_or_seqchange"


def test_compare_missing_sides():
    assert compare(None, "MAAA") == "noncoding_side"
    assert compare("MAAA", None) == "noncoding_side"
    assert compare(None, None) == "unresolved"


def test_compare_event_single_isoform():
    assert compare_event([], ["MAAA"], has_both_sides=False) == "single_isoform"


def test_compare_event_uses_representatives():
    res = compare_event(["MAAAAAAAAAKDDDD"], ["MAAA", "MAAAAAAAAAK"],
                        has_both_sides=True)
    assert res == "length_change(-4aa)"
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_comparison.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'altsplice_protein.comparison'`.

- [ ] **Step 3: Implement `comparison.py`**

Create `src/altsplice_protein/comparison.py`:

```python
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
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_comparison.py -v`
Expected: PASS (7 passed).

- [ ] **Step 5: Commit**

```bash
cd /home/qshao/altSp
git add src/altsplice_protein/comparison.py tests/test_comparison.py
git commit -m "feat: before/after protein comparison"
```

---

### Task 5: Pipeline orchestration and output writers

**Files:**
- Create: `src/altsplice_protein/pipeline.py`
- Test: `tests/test_pipeline.py`

**Interfaces:**
- Consumes: `iter_significant_events`, `collect_unique_isoforms` (filtering); `build_resolver_map`, `download`, `GTF_URL`, `PEP_URL` (ensembl_data); `resolve_all` (resolver); `compare_event`, `representative` (comparison).
- Produces: `write_fasta(path, resolved)`, `write_events_tsv(path, events, resolved)`, `write_unresolved(path, resolved)`, `run_with_files(csv_path, gtf_path, pep_path, results_dir, fdr_max=0.05, dpsi_min=0.1, use_rest=True)->dict`, `run(csv_path, data_dir="data", results_dir="results", fdr_max=0.05, dpsi_min=0.1, use_rest=True)->dict`. The returned dict has keys `events, unique_isoforms, protein, noncoding, unresolved`.

- [ ] **Step 1: Write the failing integration test**

Create `tests/test_pipeline.py` (runs the full pipeline on fixtures, `use_rest=False`):

```python
from pathlib import Path
from altsplice_protein.pipeline import run_with_files

FIX = Path(__file__).parent / "fixtures"


def test_run_with_files_end_to_end(tmp_path):
    stats = run_with_files(
        csv_path=FIX / "mini.csv",
        gtf_path=FIX / "mini.gtf",
        pep_path=FIX / "mini.pep.fa",
        results_dir=tmp_path,
        use_rest=False,
    )
    assert stats["events"] == 3
    assert stats["unique_isoforms"] == 5
    assert stats["protein"] == 4     # TST-001/002/003, NCG-001
    assert stats["noncoding"] == 1   # NCG-002

    fasta = (tmp_path / "proteins.fasta").read_text()
    assert ">TST-001|ENST001|ENSP001|TST" in fasta
    assert "MAAAAAAAAAK" in fasta
    assert "NCG-002" not in fasta    # non-coding excluded from FASTA

    tsv = (tmp_path / "events_proteins.tsv").read_text().splitlines()
    header = tsv[0].split("\t")
    rows = {r.split("\t")[0]: dict(zip(header, r.split("\t"))) for r in tsv[1:]}

    assert rows["TST_ES_1"]["comparison"] == "length_change(-4aa)"
    assert rows["TST_ES_1"]["has_complete_pair"] == "yes"
    assert rows["TST_ES_1"]["before_len"] == "15"
    assert rows["TST_ES_1"]["after_len"] == "11"

    assert rows["TST_AP_1"]["comparison"] == "single_isoform"
    assert rows["TST_AP_1"]["has_complete_pair"] == "no"

    assert rows["NC_ES_1"]["comparison"] == "noncoding_side"

    unresolved = (tmp_path / "unresolved.txt").read_text()
    assert "NCG-002" in unresolved   # reported as noncoding
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_pipeline.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'altsplice_protein.pipeline'`.

- [ ] **Step 3: Implement `pipeline.py`**

Create `src/altsplice_protein/pipeline.py`:

```python
from __future__ import annotations
import os
from .filtering import iter_significant_events, collect_unique_isoforms
from .ensembl_data import download, build_resolver_map, GTF_URL, PEP_URL
from .resolver import resolve_all
from .comparison import compare_event, representative


def write_fasta(path, resolved) -> None:
    with open(path, "w") as f:
        for name in sorted(resolved):
            iso = resolved[name]
            if not iso.protein_seq:
                continue
            f.write(
                f">{iso.name}|{iso.transcript_id}|{iso.protein_id}"
                f"|{iso.gene_symbol}\n"
            )
            s = iso.protein_seq
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")


def _ids(names, resolved) -> str:
    if not names:
        return "."
    return ";".join(
        (resolved[n].protein_id or ".") if n in resolved else "."
        for n in names
    )


def write_events_tsv(path, events, resolved) -> None:
    cols = [
        "Splice_Event", "Gene_Symbol", "Splice_Type", "PSI_Difference",
        "FDR_Difference", "SpliceIn_isoforms", "SpliceIn_protein_ids",
        "SpliceOut_isoforms", "SpliceOut_protein_ids", "has_complete_pair",
        "before_len", "after_len", "comparison",
    ]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for e in events:
            after_seqs = [
                resolved[n].protein_seq for n in e.splice_in if n in resolved
            ]
            before_seqs = [
                resolved[n].protein_seq for n in e.splice_out if n in resolved
            ]
            has_both = bool(e.splice_in) and bool(e.splice_out)
            comparison = compare_event(before_seqs, after_seqs, has_both)
            brep = representative(before_seqs)
            arep = representative(after_seqs)
            f.write("\t".join([
                e.splice_event, e.gene_symbol, e.splice_type,
                f"{e.psi_difference:g}", f"{e.fdr_difference:g}",
                ",".join(e.splice_in) or ".", _ids(e.splice_in, resolved),
                ",".join(e.splice_out) or ".", _ids(e.splice_out, resolved),
                "yes" if has_both else "no",
                str(len(brep)) if brep else ".",
                str(len(arep)) if arep else ".",
                comparison,
            ]) + "\n")


def write_unresolved(path, resolved) -> None:
    with open(path, "w") as f:
        f.write("isoform_name\ttranscript_id\tstatus\n")
        for name in sorted(resolved):
            iso = resolved[name]
            if iso.status != "protein":
                f.write(f"{name}\t{iso.transcript_id or '.'}\t{iso.status}\n")


def _stats(events, names, resolved) -> dict:
    return {
        "events": len(events),
        "unique_isoforms": len(names),
        "protein": sum(1 for i in resolved.values() if i.status == "protein"),
        "noncoding": sum(1 for i in resolved.values() if i.status == "noncoding"),
        "unresolved": sum(1 for i in resolved.values() if i.status == "unresolved"),
    }


def run_with_files(
    csv_path, gtf_path, pep_path, results_dir,
    fdr_max: float = 0.05, dpsi_min: float = 0.1, use_rest: bool = True,
) -> dict:
    os.makedirs(results_dir, exist_ok=True)
    resolver_map = build_resolver_map(gtf_path, pep_path)
    events = list(iter_significant_events(csv_path, fdr_max, dpsi_min))
    names = collect_unique_isoforms(events)
    resolved = resolve_all(names, resolver_map, use_rest=use_rest)
    write_fasta(os.path.join(results_dir, "proteins.fasta"), resolved)
    write_events_tsv(
        os.path.join(results_dir, "events_proteins.tsv"), events, resolved
    )
    write_unresolved(
        os.path.join(results_dir, "unresolved.txt"), resolved
    )
    return _stats(events, names, resolved)


def run(
    csv_path, data_dir: str = "data", results_dir: str = "results",
    fdr_max: float = 0.05, dpsi_min: float = 0.1, use_rest: bool = True,
) -> dict:
    os.makedirs(data_dir, exist_ok=True)
    gtf = download(GTF_URL, os.path.join(data_dir, "GRCh37.75.gtf.gz"))
    pep = download(PEP_URL, os.path.join(data_dir, "GRCh37.75.pep.all.fa.gz"))
    return run_with_files(
        csv_path, gtf, pep, results_dir, fdr_max, dpsi_min, use_rest
    )
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_pipeline.py -v`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
cd /home/qshao/altSp
git add src/altsplice_protein/pipeline.py tests/test_pipeline.py
git commit -m "feat: pipeline orchestration and output writers"
```

---

### Task 6: CLI entry point and real-data run

**Files:**
- Create: `src/altsplice_protein/__main__.py`
- Create: `README.md`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `run` from `pipeline.py`.
- Produces: `main(argv:list[str]|None=None)->dict` (parses args, calls `run`, prints stats) and a `python -m altsplice_protein` entry point.

- [ ] **Step 1: Write the failing CLI test**

Create `tests/test_cli.py` (patches `run` so no download happens):

```python
import altsplice_protein.__main__ as cli


def test_main_parses_args_and_calls_run(monkeypatch, capsys):
    captured = {}

    def fake_run(csv_path, data_dir, results_dir, fdr_max, dpsi_min, use_rest):
        captured.update(dict(
            csv_path=csv_path, data_dir=data_dir, results_dir=results_dir,
            fdr_max=fdr_max, dpsi_min=dpsi_min, use_rest=use_rest,
        ))
        return {"events": 3, "protein": 4}

    monkeypatch.setattr(cli, "run", fake_run)
    stats = cli.main([
        "spliceseq_info_PRAD.csv", "--fdr-max", "0.01",
        "--dpsi-min", "0.2", "--no-rest",
    ])
    assert captured["csv_path"] == "spliceseq_info_PRAD.csv"
    assert captured["fdr_max"] == 0.01
    assert captured["dpsi_min"] == 0.2
    assert captured["use_rest"] is False
    assert stats["events"] == 3
    out = capsys.readouterr().out
    assert "events: 3" in out
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_cli.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'altsplice_protein.__main__'`.

- [ ] **Step 3: Implement `__main__.py`**

Create `src/altsplice_protein/__main__.py`:

```python
from __future__ import annotations
import argparse
from .pipeline import run


def main(argv: list[str] | None = None) -> dict:
    p = argparse.ArgumentParser(
        prog="altsplice_protein",
        description="Resolve before/after protein sequences for significant "
                    "SpliceSeq events.",
    )
    p.add_argument("csv", help="Path to spliceseq_info_*.csv")
    p.add_argument("--data-dir", default="data",
                   help="Where Ensembl files are downloaded/cached")
    p.add_argument("--results-dir", default="results",
                   help="Where output files are written")
    p.add_argument("--fdr-max", type=float, default=0.05)
    p.add_argument("--dpsi-min", type=float, default=0.1)
    p.add_argument("--no-rest", action="store_true",
                   help="Disable Ensembl REST fallback for unresolved names")
    a = p.parse_args(argv)
    stats = run(
        a.csv, a.data_dir, a.results_dir, a.fdr_max, a.dpsi_min, not a.no_rest
    )
    for k, v in stats.items():
        print(f"{k}: {v}")
    return stats


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd /home/qshao/altSp && python3 -m pytest tests/test_cli.py -v`
Expected: PASS (1 passed).

- [ ] **Step 5: Run the full unit-test suite**

Run: `cd /home/qshao/altSp && python3 -m pytest -v`
Expected: PASS (all tests across all files).

- [ ] **Step 6: Write the README**

Create `README.md`:

```markdown
# Alternative Splicing → Protein Sequences (PRAD)

Resolves the "before" (spliced-out) and "after" (spliced-in) protein sequences
for significant alternative splicing events in `spliceseq_info_PRAD.csv`
(TCGA SpliceSeq). Isoform names are resolved against Ensembl GRCh37 release-75.

## Run

    python -m altsplice_protein spliceseq_info_PRAD.csv

First run downloads ~50 MB of Ensembl GRCh37 r75 files into `data/`.
Outputs are written to `results/`:

- `proteins.fasta` — one protein per resolved isoform.
- `events_proteins.tsv` — per-event before/after proteins and comparison.
- `unresolved.txt` — isoform names with no protein (non-coding or unmapped).

Options: `--fdr-max`, `--dpsi-min`, `--results-dir`, `--data-dir`, `--no-rest`.

## Develop

    PYTHONPATH=src python -m pytest
```

- [ ] **Step 7: Execute on the real dataset**

Run: `cd /home/qshao/altSp && PYTHONPATH=src python3 -m altsplice_protein spliceseq_info_PRAD.csv`
Expected: downloads the two Ensembl files (first run only), then prints stats. Expect roughly `events: 2678`, `unique_isoforms: 5845`, and a high `protein` count. Note the actual numbers.

- [ ] **Step 8: Validate the output**

Run these checks and confirm:

```bash
cd /home/qshao/altSp
head -1 results/events_proteins.tsv
grep -P "\tA2M_ES_20222\t|^A2M_ES_20222\t" results/events_proteins.tsv | head
grep -c ">" results/proteins.fasta
wc -l results/unresolved.txt
# spot-check a known mapping resolved to a protein:
grep "^>A4GALT-003|" results/proteins.fasta
```

Expected: the TSV header matches the 13 columns; `A2M_ES_20222` row shows `has_complete_pair=yes` with before/after protein IDs; `proteins.fasta` has thousands of entries; `A4GALT-003` appears with an `ENST`/`ENSP` id. Report the resolution coverage (protein vs noncoding vs unresolved) from the stats.

- [ ] **Step 9: Commit**

```bash
cd /home/qshao/altSp
git add src/altsplice_protein/__main__.py tests/test_cli.py README.md
git commit -m "feat: CLI entry point, README, and real-data run"
```

---

## Self-Review Notes

- **Spec coverage:** filter (Task 1); local GTF/pep resolution (Task 2); REST fallback (Task 3); before/after comparison + non-coding/single-isoform/unresolved flags (Task 4); three output files with specified headers (Task 5); validation/spot-checks (Task 6, Steps 7–8). Comparison enum refined vs. spec: `frameshift_likely` is implemented as the precise, testable `frameshift_or_seqchange` (sequences differ and the shorter is not a clean prefix of the longer); clean in-frame exon gain/loss is `length_change(±Naa)`.
- **Provenance:** annotation pinned to Ensembl GRCh37 **release-75**; URLs verified to return HTTP 200.
- **Network isolation in tests:** all unit/integration tests run offline (`use_rest=False` or monkeypatched REST); only Task 6 Step 7 touches the network.

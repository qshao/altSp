# Collect & Integrate More PCa Splicing — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add three new prostate-cancer (PCa) alternative-splicing sources (literature/EuropePMC, ASCancerAtlas, GTEx normal baseline) and a UniProt-feature-driven functional-impact analysis, then merge everything (incl. existing TCGA + UniProt) into one master variant table and one integrated report.

**Architecture:** Modular collectors write intermediate TSVs; a shared UniProt-feature annotator computes functional impact (domain mapping, NMD heuristic, disorder/localization) by overlapping each variant's changed residue interval with UniProt feature tracks; an integrator normalizes all sources into one schema with a cross-source corroboration count; a report builder emits Word + Markdown. Pure logic is unit-tested offline; network calls are disk-cached with an `--offline` replay mode.

**Tech Stack:** Python 3, requests, pandas, python-docx, pytest. UniProtKB / EuropePMC / GTEx v2 / ASCancerAtlas REST+download.

## Global Constraints

- New importable code lives under `analysis/` as packages: `analysis/collect/`, `analysis/annotate/`, `analysis/integrate/` (each with `__init__.py`).
- `pytest.ini` `pythonpath` must include both `src` and `analysis`.
- All network calls cache to `results_collected/cache/<source>/`; never re-hit an API when a cache file exists; `--offline` reads only cache/fixtures.
- A missing/blocked source is a warning + provenance note, never a crash (CancerSplicingQTL is HTTP 403 → permanently skipped with a note).
- Outputs go to `results_collected/` (data) and `docs/` (reports). `results_collected/` is git-ignored except curated final TSVs explicitly added.
- All disease/resistance/NMD/domain calls are predictions/hypotheses — label them as such in any prose.
- Commit after each task. Commit message trailer: `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.

---

### Task 1: Package scaffolding + UniProt feature-track parser/cache

**Files:**
- Modify: `pytest.ini`
- Create: `analysis/__init__.py`, `analysis/annotate/__init__.py`, `analysis/collect/__init__.py`, `analysis/integrate/__init__.py`
- Create: `analysis/annotate/uniprot_features.py`
- Create: `tests/fixtures/uniprot_AR_record.json` (one real record extracted from `uniprot_prad_splicing/proteins.json`)
- Test: `tests/test_uniprot_features.py`

**Interfaces:**
- Produces:
  - `Feature = namedtuple("Feature", "type start end description")`
  - `parse_features(record: dict) -> list[Feature]` — pulls Domain, Region, Binding site, Active site, Signal, Transmembrane, Topological domain, Motif, plus disordered Regions, from a UniProt JSON record's `features`.
  - `fetch_features(acc: str, cache_dir: str, offline: bool=False) -> list[Feature]` — GET `https://rest.uniprot.org/uniprotkb/{acc}.json`, cache raw JSON to `{cache_dir}/{acc}.json`, return `parse_features`.

- [ ] **Step 1: Add `analysis` to pytest pythonpath**

Edit `pytest.ini`:
```ini
[pytest]
pythonpath = src analysis
testpaths = tests
```

- [ ] **Step 2: Create empty package files**

```bash
touch analysis/__init__.py analysis/annotate/__init__.py analysis/collect/__init__.py analysis/integrate/__init__.py
```

- [ ] **Step 3: Extract the AR fixture record**

```bash
python -c "import json; r=[x for x in json.load(open('uniprot_prad_splicing/proteins.json')) if x['primaryAccession']=='P10275'][0]; json.dump(r, open('tests/fixtures/uniprot_AR_record.json','w'))"
```

- [ ] **Step 4: Write the failing test**

```python
# tests/test_uniprot_features.py
import json
from annotate.uniprot_features import parse_features, Feature


def test_parse_features_extracts_typed_intervals():
    rec = json.load(open("tests/fixtures/uniprot_AR_record.json"))
    feats = parse_features(rec)
    assert feats, "expected at least one feature"
    assert all(isinstance(f, Feature) for f in feats)
    # AR has a curated NR (nuclear receptor) DNA-binding / ligand-binding region
    assert all(f.start <= f.end for f in feats)
    types = {f.type for f in feats}
    # AR record carries Domain and/or Region features
    assert types & {"Domain", "Region", "Zinc finger", "DNA-binding"}
```

- [ ] **Step 5: Run test, verify it fails**

Run: `pytest tests/test_uniprot_features.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'annotate.uniprot_features'`

- [ ] **Step 6: Implement `uniprot_features.py`**

```python
"""Fetch and parse UniProt feature tracks used for functional-impact analysis."""
from __future__ import annotations
import json
import os
from collections import namedtuple

import requests

Feature = namedtuple("Feature", "type start end description")

# UniProt feature types we use, grouped by downstream purpose.
DOMAIN_TYPES = {"Domain", "Region", "Repeat", "Zinc finger", "DNA-binding",
                "Coiled coil", "Motif", "Compositional bias"}
SITE_TYPES = {"Binding site", "Active site", "Site"}
LOC_TYPES = {"Signal", "Transit peptide", "Transmembrane",
             "Topological domain", "Intramembrane"}
KEEP = DOMAIN_TYPES | SITE_TYPES | LOC_TYPES


def parse_features(record: dict) -> list[Feature]:
    out: list[Feature] = []
    for f in record.get("features", []):
        ftype = f.get("type", "")
        if ftype not in KEEP:
            continue
        loc = f.get("location", {})
        s = (loc.get("start") or {}).get("value")
        e = (loc.get("end") or {}).get("value")
        if s is None or e is None:
            continue
        out.append(Feature(ftype, int(s), int(e), f.get("description", "")))
    return out


def fetch_features(acc: str, cache_dir: str, offline: bool = False) -> list[Feature]:
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{acc}.json")
    if os.path.exists(path):
        return parse_features(json.load(open(path)))
    if offline:
        return []
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    rec = r.json()
    json.dump(rec, open(path, "w"))
    return parse_features(rec)
```

- [ ] **Step 7: Run test, verify it passes**

Run: `pytest tests/test_uniprot_features.py -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add pytest.ini analysis/__init__.py analysis/annotate analysis/collect analysis/integrate tests/fixtures/uniprot_AR_record.json tests/test_uniprot_features.py
git commit -m "feat(annotate): UniProt feature-track parser + disk cache

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Functional-impact analysis (changed interval, domains, NMD, localization, disorder)

**Files:**
- Create: `analysis/annotate/functional_impact.py`
- Test: `tests/test_functional_impact.py`

**Interfaces:**
- Consumes: `Feature` from `annotate.uniprot_features`.
- Produces:
  - `changed_interval(before: str, after: str) -> tuple[int,int]|None` — 1-based inclusive residue range that differs (using common prefix/suffix); `None` if identical or a side is empty.
  - `features_in(interval, feats, types) -> list[str]` — descriptions (or type) of features overlapping the interval, filtered to `types`.
  - `nmd_flag(before_len:int, after_len:int, change_class:str, common_suffix:int) -> str` — `"NMD-candidate"|"no"|"n.a."`.
  - `annotate_variant(before:str|None, after:str|None, change_class:str, feats:list[Feature]) -> dict` with keys `changed_interval, domains_hit, regions_lost, loc_flags, disorder_overlap, nmd_flag`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_functional_impact.py
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
```

- [ ] **Step 2: Run tests, verify they fail**

Run: `pytest tests/test_functional_impact.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `functional_impact.py`**

```python
"""Predict how a splice change affects protein function, from UniProt features."""
from __future__ import annotations
from annotate.uniprot_features import Feature, DOMAIN_TYPES, SITE_TYPES, LOC_TYPES

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
    if "truncation" in change_class and common_suffix == 0 and lost_cterm > NMD_TAIL_AA:
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
```

- [ ] **Step 4: Run tests, verify they pass**

Run: `pytest tests/test_functional_impact.py -v`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add analysis/annotate/functional_impact.py tests/test_functional_impact.py
git commit -m "feat(annotate): functional-impact (domains/NMD/localization/disorder)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: Master-table schema + normalize existing TCGA & UniProt sources

**Files:**
- Create: `analysis/integrate/master_table.py`
- Test: `tests/test_master_table.py`

**Interfaces:**
- Consumes: `annotate_variant`, `fetch_features`.
- Produces:
  - `COLUMNS: list[str]` — the master schema column order (matches the spec table).
  - `normalize_tcga(row: dict) -> dict` — one row of `analysis/prioritized_events.tsv` (+ matching `sequence_impact.json` entry) → master row.
  - `normalize_uniprot(rec: dict, iso: dict) -> dict` — one isoform of `uniprot_sequence_impact.json` → master row.
  - `add_corroboration(rows: list[dict]) -> list[dict]` — fill `corroborating_sources` = count of distinct `source` per `gene`.
  - `build_master(offline: bool=False) -> list[dict]` — load both existing sources, annotate with UniProt features, return master rows; CLI `main()` writes `results_collected/master_variants.tsv`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_master_table.py
from integrate.master_table import (
    COLUMNS, normalize_uniprot, add_corroboration,
)


def test_columns_contains_core_schema():
    for c in ["variant_id", "gene", "uniprot", "source", "change_class",
              "domains_hit", "nmd_flag", "corroborating_sources", "provenance"]:
        assert c in COLUMNS


def test_normalize_uniprot_row_shape():
    rec = {"gene": "AR", "accession": "P10275", "canonical_len": 920,
           "canonical_seq": "M" * 920,
           "is_prostate_disease": True,
           "diseases": [{"id": "Prostate cancer, hereditary, X-linked 3",
                         "desc": ""}]}
    iso = {"isoform_id": "P10275-3", "isoform_name": "3",
           "synonyms": ["AR-V7"], "after_len": 644, "after_seq": "M" * 644,
           "change_class": "substituted segment", "varseq": []}
    row = normalize_uniprot(rec, iso)
    assert row["gene"] == "AR"
    assert row["source"] == "UniProt"
    assert row["uniprot"] == "P10275"
    assert row["after_len"] == 644
    assert "AR-V7" in row["provenance"] or "AR-V7" in row["source_id"]
    assert set(COLUMNS) <= set(row)


def test_add_corroboration_counts_distinct_sources():
    rows = [
        {"gene": "ITGA6", "source": "UniProt"},
        {"gene": "ITGA6", "source": "TCGA"},
        {"gene": "AR", "source": "UniProt"},
    ]
    out = add_corroboration(rows)
    by = {(r["gene"], r["source"]): r["corroborating_sources"] for r in out}
    assert by[("ITGA6", "UniProt")] == 2
    assert by[("AR", "UniProt")] == 1
```

- [ ] **Step 2: Run tests, verify they fail**

Run: `pytest tests/test_master_table.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `master_table.py`**

```python
"""Normalize every splicing source into one master variant table."""
from __future__ import annotations
import json
import os

import pandas as pd

from annotate.uniprot_features import fetch_features
from annotate.functional_impact import annotate_variant

OUTDIR = "results_collected"
CACHE = os.path.join(OUTDIR, "cache", "uniprot_features")
TCGA_EVENTS = "analysis/prioritized_events.tsv"
TCGA_IMPACT = "analysis/sequence_impact.json"
UNI_IMPACT = "analysis/uniprot_sequence_impact.json"

COLUMNS = [
    "variant_id", "gene", "uniprot", "source", "source_id", "event_type",
    "status", "before_len", "after_len", "change_class", "changed_interval",
    "domains_hit", "regions_lost", "loc_flags", "disorder_overlap", "nmd_flag",
    "pca_evidence", "gtex_baseline", "corroborating_sources", "provenance",
]


def _blank() -> dict:
    return {c: "" for c in COLUMNS}


def _fmt(v) -> str:
    if isinstance(v, (list, tuple)):
        return "; ".join(str(x) for x in v) if v else ""
    return "" if v is None else str(v)


def normalize_uniprot(rec: dict, iso: dict, feats=None) -> dict:
    row = _blank()
    syn = "/".join(iso.get("synonyms", []))
    ann = annotate_variant(rec.get("canonical_seq"), iso.get("after_seq"),
                           iso.get("change_class", ""), feats or [])
    dis = [d["id"] for d in rec.get("diseases", [])
           if "prostate" in (d["id"] + d.get("desc", "")).lower()]
    row.update(
        variant_id=f"{rec['gene']}|UniProt|{iso['isoform_id']}",
        gene=rec["gene"], uniprot=rec.get("accession", ""), source="UniProt",
        source_id=f"{iso['isoform_id']} ({syn})" if syn else iso["isoform_id"],
        event_type="isoform", status=iso.get("status", "curated isoform"),
        before_len=rec.get("canonical_len", ""), after_len=iso.get("after_len", ""),
        change_class=iso.get("change_class", ""),
        changed_interval=_fmt(ann["changed_interval"]),
        domains_hit=_fmt(ann["domains_hit"]), regions_lost=_fmt(ann["regions_lost"]),
        loc_flags=_fmt(ann["loc_flags"]), disorder_overlap=ann["disorder_overlap"],
        nmd_flag=ann["nmd_flag"],
        pca_evidence="; ".join(dis) if dis else ("PCa-associated" if rec.get("is_prostate_disease") else ""),
        provenance=f"UniProt {rec.get('accession','')} {syn}".strip(),
    )
    return row


def normalize_tcga(row: dict, impact: dict | None, feats=None) -> dict:
    out = _blank()
    before = (impact or {}).get("before_seq")
    after = (impact or {}).get("after_seq")
    cls = (impact or {}).get("impact_class", row.get("protein_change", ""))
    ann = annotate_variant(before, after, cls, feats or [])
    out.update(
        variant_id=f"{row['Gene']}|TCGA|{row['Splice_Event']}",
        gene=row["Gene"], uniprot=(impact or {}).get("uniprot", ""), source="TCGA",
        source_id=row["Splice_Event"], event_type=row.get("Type", ""),
        status=f"{row.get('endpoint','')} HR {row.get('HR','')} ({row.get('direction','')}), BHp {row.get('bh_p','')}",
        before_len=(impact or {}).get("before_len", ""),
        after_len=(impact or {}).get("after_len", ""),
        change_class=cls, changed_interval=_fmt(ann["changed_interval"]),
        domains_hit=_fmt(ann["domains_hit"]), regions_lost=_fmt(ann["regions_lost"]),
        loc_flags=_fmt(ann["loc_flags"]), disorder_overlap=ann["disorder_overlap"],
        nmd_flag=ann["nmd_flag"],
        pca_evidence=f"cohort {row.get('endpoint','')} HR {row.get('HR','')}",
        provenance=f"TCGA SpliceSeq PRAD {row['Splice_Event']}",
    )
    return out


def add_corroboration(rows: list[dict]) -> list[dict]:
    from collections import defaultdict
    srcs = defaultdict(set)
    for r in rows:
        srcs[r["gene"]].add(r["source"])
    for r in rows:
        r["corroborating_sources"] = len(srcs[r["gene"]])
    return rows


def build_master(offline: bool = False) -> list[dict]:
    rows: list[dict] = []
    feat_cache: dict[str, list] = {}

    def feats_for(acc):
        if not acc:
            return []
        if acc not in feat_cache:
            feat_cache[acc] = fetch_features(acc, CACHE, offline=offline)
        return feat_cache[acc]

    # UniProt source
    uni = json.load(open(UNI_IMPACT))
    for rec in uni:
        feats = feats_for(rec.get("accession", ""))
        for iso in rec.get("isoforms", []):
            rows.append(normalize_uniprot(rec, iso, feats))

    # TCGA source (join events to sequence-impact by splice_event)
    impacts = {o["splice_event"]: o for o in json.load(open(TCGA_IMPACT))}
    ev = pd.read_csv(TCGA_EVENTS, sep="\t")
    for _, r in ev.iterrows():
        imp = impacts.get(r["Splice_Event"])
        acc = (imp or {}).get("uniprot", "")
        rows.append(normalize_tcga(r.to_dict(), imp, feats_for(acc)))

    return add_corroboration(rows)


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    rows = build_master(offline=args.offline)
    df = pd.DataFrame(rows, columns=COLUMNS)
    out = os.path.join(OUTDIR, "master_variants.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(f"rows: {len(df)}  genes: {df['gene'].nunique()}  -> {out}")
    print(df["source"].value_counts().to_string())


if __name__ == "__main__":
    main()
```

Note on TCGA `uniprot` mapping: `sequence_impact.json` may not carry a `uniprot` field. If absent, `feats_for("")` returns `[]` and annotation columns stay empty for TCGA rows — acceptable; UniProt-source rows carry the domain analysis. (A later improvement could map TCGA gene→accession.)

- [ ] **Step 4: Run unit tests, verify they pass**

Run: `pytest tests/test_master_table.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Build the table for real (network, cached)**

Run: `python analysis/integrate/master_table.py`
Expected: prints row/gene counts and a source breakdown; writes `results_collected/master_variants.tsv`. Sanity-check AR-V7 has non-empty `domains_hit`:
```bash
python -c "import pandas as pd; d=pd.read_csv('results_collected/master_variants.tsv',sep='\t'); print(d[d.source_id.astype(str).str.contains('AR-V7')][['gene','domains_hit','nmd_flag']].to_string())"
```

- [ ] **Step 6: Commit**

```bash
git add analysis/integrate/master_table.py tests/test_master_table.py
git commit -m "feat(integrate): master variant table + TCGA/UniProt normalization

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Literature-curated variants collector (EuropePMC provenance)

**Files:**
- Create: `analysis/collect/literature_variants.py`
- Test: `tests/test_literature_variants.py`

**Interfaces:**
- Produces:
  - `SEED: list[dict]` — curated PCa splice variants, each `{gene, uniprot, variant, event_type, effect, pca_note, query}`.
  - `parse_europepmc(payload: dict) -> list[str]` — extract PMIDs from a EuropePMC search JSON.
  - `pmids_for(query, cache_dir, offline=False) -> list[str]` — cached EuropePMC lookup (top 3 PMIDs).
  - `collect(cache_dir, offline=False) -> list[dict]` — SEED rows enriched with PMIDs; CLI `main()` writes `results_collected/raw/literature_variants.tsv`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_literature_variants.py
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
```

- [ ] **Step 2: Run tests, verify they fail**

Run: `pytest tests/test_literature_variants.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `literature_variants.py`**

```python
"""Curated PCa splice variants + EuropePMC provenance lookup."""
from __future__ import annotations
import csv
import json
import os

import requests

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Curated, well-characterized prostate-cancer splice variants (literature).
SEED = [
    dict(gene="AR", uniprot="P10275", variant="AR-V7",
         event_type="cryptic exon (CE3)", effect="LBD-truncated, constitutively active",
         pca_note="enzalutamide/abiraterone resistance biomarker",
         query="AR-V7 prostate cancer splice variant"),
    dict(gene="AR", uniprot="P10275", variant="AR-V567es",
         event_type="exon skipping (exons 5-7)", effect="LBD-truncated, constitutively active",
         pca_note="castration-resistant prostate cancer",
         query="AR-V567es prostate cancer"),
    dict(gene="KLF6", uniprot="Q99612", variant="KLF6-SV1",
         event_type="alternative splice", effect="oncogenic, antagonizes wild-type KLF6",
         pca_note="aggressive, treatment-refractory disease",
         query="KLF6-SV1 prostate cancer"),
    dict(gene="CCND1", uniprot="P24385", variant="cyclin D1b",
         event_type="intron 4 retention", effect="loss of Thr286, nuclear retention",
         pca_note="enhanced AR coactivation / proliferation",
         query="cyclin D1b prostate cancer splice"),
    dict(gene="FGFR2", uniprot="P21802", variant="FGFR2-IIIb/IIIc switch",
         event_type="mutually exclusive exons", effect="ligand-specificity switch",
         pca_note="EMT / progression",
         query="FGFR2 IIIb IIIc prostate cancer splicing"),
    dict(gene="BCL2L1", uniprot="Q07817", variant="Bcl-xS",
         event_type="alternative 5' splice site", effect="pro-apoptotic short isoform",
         pca_note="apoptosis balance / therapy response",
         query="Bcl-xS Bcl-xL prostate cancer splicing"),
    dict(gene="ERG", uniprot="P11308", variant="TMPRSS2-ERG isoforms",
         event_type="fusion-derived alternative splicing", effect="oncogenic ETS activation",
         pca_note="hallmark prostate-cancer fusion",
         query="TMPRSS2 ERG splice isoform prostate cancer"),
    dict(gene="PKM", uniprot="P14618", variant="PKM2",
         event_type="mutually exclusive exons (9/10)", effect="glycolytic isoform switch",
         pca_note="Warburg metabolism in tumour",
         query="PKM2 prostate cancer splicing"),
]


def parse_europepmc(payload: dict) -> list[str]:
    res = (payload.get("resultList") or {}).get("result", [])
    return [r.get("pmid") for r in res if r.get("pmid")]


def pmids_for(query: str, cache_dir: str, offline: bool = False) -> list[str]:
    os.makedirs(cache_dir, exist_ok=True)
    key = "".join(c if c.isalnum() else "_" for c in query)[:60]
    path = os.path.join(cache_dir, f"{key}.json")
    if os.path.exists(path):
        return parse_europepmc(json.load(open(path)))[:3]
    if offline:
        return []
    r = requests.get(EPMC, params={"query": query, "format": "json",
                                    "pageSize": 3, "resultType": "lite"}, timeout=60)
    r.raise_for_status()
    payload = r.json()
    json.dump(payload, open(path, "w"))
    return parse_europepmc(payload)[:3]


def collect(cache_dir: str, offline: bool = False) -> list[dict]:
    out = []
    for s in SEED:
        pmids = pmids_for(s["query"], cache_dir, offline=offline)
        row = dict(s)
        row["pmids"] = ";".join(pmids)
        out.append(row)
    return out


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    cache = "results_collected/cache/europepmc"
    rows = collect(cache, offline=args.offline)
    os.makedirs("results_collected/raw", exist_ok=True)
    out = "results_collected/raw/literature_variants.tsv"
    cols = ["gene", "uniprot", "variant", "event_type", "effect", "pca_note",
            "pmids", "query"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"literature variants: {len(rows)} -> {out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run unit tests, verify they pass**

Run: `pytest tests/test_literature_variants.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Run the collector for real (network, cached)**

Run: `python analysis/collect/literature_variants.py`
Expected: writes `results_collected/raw/literature_variants.tsv` with 8 rows, most carrying PMIDs.

- [ ] **Step 6: Commit**

```bash
git add analysis/collect/literature_variants.py tests/test_literature_variants.py
git commit -m "feat(collect): literature-curated PCa splice variants + EuropePMC provenance

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: ASCancerAtlas PRAD collector (download + parse, graceful)

**Files:**
- Create: `analysis/collect/db_ascanceratlas.py`
- Test: `tests/test_db_ascanceratlas.py`
- Create: `tests/fixtures/ascancer_sample.tsv` (small hand-made sample matching the real column layout discovered at runtime)

**Interfaces:**
- Produces:
  - `parse_table(text: str) -> list[dict]` — parse a delimited ASCancerAtlas dump into `{gene, event_type, cancer, as_id, note}` rows, filtered to prostate (`PRAD`/`prostate`).
  - `collect(cache_dir, offline=False) -> list[dict]` — download (cached) + parse; on any network/format failure return `[]` and print a warning (never raise). CLI `main()` writes `results_collected/raw/ascanceratlas_prad.tsv`.

**Discovery note (runtime):** Before writing `parse_table`, fetch the ASCancerAtlas `Download` page to find the actual bulk-file URL and delimiter/columns (`curl -s https://ngdc.cncb.ac.cn/ascancer/Download`). Build `tests/fixtures/ascancer_sample.tsv` from 3-5 real lines so the parser test reflects the true layout. If no prostate data or no machine-readable dump is available, implement `collect` to return `[]` with a printed note and record that in the integrated report's limitations — do not block the pipeline.

- [ ] **Step 1: Discover the real download format**

Run:
```bash
curl -s --max-time 30 https://ngdc.cncb.ac.cn/ascancer/Download | grep -ioE 'href="[^"]+"' | grep -iE 'download|txt|tsv|csv|zip|gz' | head
```
Record the chosen URL + delimiter in a comment at the top of `db_ascanceratlas.py`. Create `tests/fixtures/ascancer_sample.tsv` from real sample lines (header + 3-5 rows including >=1 prostate row).

- [ ] **Step 2: Write the failing test** (adjust column names to the discovered layout)

```python
# tests/test_db_ascanceratlas.py
from collect.db_ascanceratlas import parse_table


def test_parse_table_keeps_prostate_rows():
    text = open("tests/fixtures/ascancer_sample.tsv").read()
    rows = parse_table(text)
    assert rows, "expected >=1 prostate row"
    assert all("gene" in r and "event_type" in r for r in rows)
    assert all("prostate" in (r["cancer"].lower()) or "prad" in r["cancer"].lower()
               for r in rows)
```

- [ ] **Step 3: Run test, verify it fails**

Run: `pytest tests/test_db_ascanceratlas.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 4: Implement `db_ascanceratlas.py`** (column indices per discovery)

```python
"""ASCancerAtlas prostate (PRAD) alternative-splicing events.

Download URL + format discovered from https://ngdc.cncb.ac.cn/ascancer/Download
(see Step 1). Parser is tested offline against tests/fixtures/ascancer_sample.tsv.
Graceful: any failure -> [] + warning, never raises.
"""
from __future__ import annotations
import csv
import io
import os

import requests

DOWNLOAD_URL = "<FILL FROM STEP 1>"   # exact bulk-file URL
# Column names as they appear in the dump header (adjust to discovery):
GENE_COL, EVENT_COL, CANCER_COL, ID_COL = "symbol", "as_type", "cancer", "as_id"


def parse_table(text: str) -> list[dict]:
    rows = []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for r in reader:
        cancer = (r.get(CANCER_COL) or "").strip()
        if "prad" not in cancer.lower() and "prostate" not in cancer.lower():
            continue
        rows.append({
            "gene": (r.get(GENE_COL) or "").strip(),
            "event_type": (r.get(EVENT_COL) or "").strip(),
            "cancer": cancer,
            "as_id": (r.get(ID_COL) or "").strip(),
            "note": "ASCancerAtlas",
        })
    return rows


def collect(cache_dir: str, offline: bool = False) -> list[dict]:
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, "ascancer_dump.tsv")
    try:
        if os.path.exists(path):
            text = open(path).read()
        elif offline:
            return []
        else:
            r = requests.get(DOWNLOAD_URL, timeout=120)
            r.raise_for_status()
            text = r.text
            open(path, "w").write(text)
        return parse_table(text)
    except Exception as e:  # graceful degradation per spec
        print(f"[ascanceratlas] skipped: {e}")
        return []


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    rows = collect("results_collected/cache/ascanceratlas", offline=args.offline)
    os.makedirs("results_collected/raw", exist_ok=True)
    out = "results_collected/raw/ascanceratlas_prad.tsv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "event_type", "cancer", "as_id", "note"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"ASCancerAtlas PRAD rows: {len(rows)} -> {out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Run unit test, verify it passes**

Run: `pytest tests/test_db_ascanceratlas.py -v`
Expected: PASS

- [ ] **Step 6: Run the collector for real**

Run: `python analysis/collect/db_ascanceratlas.py`
Expected: prints PRAD row count → `results_collected/raw/ascanceratlas_prad.tsv` (or a graceful skip note if no machine-readable prostate dump exists).

- [ ] **Step 7: Commit**

```bash
git add analysis/collect/db_ascanceratlas.py tests/test_db_ascanceratlas.py tests/fixtures/ascancer_sample.tsv
git commit -m "feat(collect): ASCancerAtlas PRAD splicing events (graceful)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: GTEx prostate normal baseline

**Files:**
- Create: `analysis/collect/gtex_prostate.py`
- Test: `tests/test_gtex_prostate.py`

**Interfaces:**
- Produces:
  - `summarize_transcripts(payload: dict) -> dict` — from a GTEx medianTranscriptExpression JSON, return `{gencodeId, n_transcripts, max_median, second_median, multi_isoform_normal: bool}` (multi_isoform_normal = 2nd-highest median ≥ 1 TPM, i.e. >1 isoform abundant in normal prostate).
  - `gene_baseline(gene, gencode_id, cache_dir, offline=False) -> dict` — cached GTEx call → summary row.
  - `collect(genes: dict[str,str], cache_dir, offline=False) -> list[dict]` — `genes` maps symbol→versioned GENCODE id; CLI `main()` reads the master table's genes, resolves GENCODE ids via the GTEx gene endpoint (cached), writes `results_collected/gtex_prostate_baseline.tsv`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_gtex_prostate.py
from collect.gtex_prostate import summarize_transcripts


def test_summarize_flags_multi_isoform_normal():
    payload = {"data": [
        {"median": 12.0, "transcriptId": "t1"},
        {"median": 4.0, "transcriptId": "t2"},
        {"median": 0.1, "transcriptId": "t3"},
    ], "gencodeId": "ENSG00000169083.16"}
    s = summarize_transcripts(payload)
    assert s["n_transcripts"] == 3
    assert s["max_median"] == 12.0
    assert s["second_median"] == 4.0
    assert s["multi_isoform_normal"] is True


def test_summarize_single_isoform():
    payload = {"data": [{"median": 9.0, "transcriptId": "t1"},
                        {"median": 0.2, "transcriptId": "t2"}],
               "gencodeId": "X"}
    s = summarize_transcripts(payload)
    assert s["multi_isoform_normal"] is False
```

- [ ] **Step 2: Run tests, verify they fail**

Run: `pytest tests/test_gtex_prostate.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `gtex_prostate.py`**

```python
"""GTEx normal-prostate isoform baseline (median transcript expression)."""
from __future__ import annotations
import csv
import json
import os

import requests

API = "https://gtexportal.org/api/v2"
TISSUE = "Prostate"
MULTI_ISOFORM_TPM = 1.0


def summarize_transcripts(payload: dict) -> dict:
    data = payload.get("data", [])
    meds = sorted((d.get("median", 0.0) for d in data), reverse=True)
    top = meds[0] if meds else 0.0
    second = meds[1] if len(meds) > 1 else 0.0
    return {
        "gencodeId": payload.get("gencodeId", ""),
        "n_transcripts": len(data),
        "max_median": top,
        "second_median": second,
        "multi_isoform_normal": second >= MULTI_ISOFORM_TPM,
    }


def _cached_get(url, params, path, offline):
    if os.path.exists(path):
        return json.load(open(path))
    if offline:
        return {}
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    payload = r.json()
    json.dump(payload, open(path, "w"))
    return payload


def resolve_gencode(gene: str, cache_dir: str, offline=False) -> str:
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"gene_{gene}.json")
    payload = _cached_get(f"{API}/reference/gene", {"geneId": gene}, path, offline)
    data = payload.get("data", [])
    return data[0].get("gencodeId", "") if data else ""


def gene_baseline(gene: str, gencode_id: str, cache_dir: str, offline=False) -> dict:
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"tx_{gene}.json")
    payload = _cached_get(f"{API}/expression/medianTranscriptExpression",
                          {"gencodeId": gencode_id, "tissueSiteDetailId": TISSUE},
                          path, offline)
    s = summarize_transcripts(payload or {"gencodeId": gencode_id, "data": []})
    s["gene"] = gene
    return s


def collect(genes: dict, cache_dir: str, offline=False) -> list[dict]:
    out = []
    for gene, gid in genes.items():
        try:
            if not gid:
                gid = resolve_gencode(gene, cache_dir, offline=offline)
            if not gid:
                continue
            out.append(gene_baseline(gene, gid, cache_dir, offline=offline))
        except Exception as e:
            print(f"[gtex] {gene} skipped: {e}")
    return out


def main() -> None:
    import argparse
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--master", default="results_collected/master_variants.tsv")
    args = ap.parse_args()
    genes = {}
    if os.path.exists(args.master):
        df = pd.read_csv(args.master, sep="\t")
        genes = {g: "" for g in sorted(df["gene"].dropna().unique())}
    rows = collect(genes, "results_collected/cache/gtex", offline=args.offline)
    out = "results_collected/gtex_prostate_baseline.tsv"
    cols = ["gene", "gencodeId", "n_transcripts", "max_median", "second_median",
            "multi_isoform_normal"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"GTEx baseline rows: {len(rows)} -> {out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run unit tests, verify they pass**

Run: `pytest tests/test_gtex_prostate.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Smoke-test the real API on a handful of genes**

Run:
```bash
python -c "from collect.gtex_prostate import resolve_gencode, gene_baseline; import sys; sys.path.insert(0,'analysis'); g=resolve_gencode('AR','results_collected/cache/gtex'); print('AR gencode',g); print(gene_baseline('AR',g,'results_collected/cache/gtex'))"
```
Expected: prints a versioned GENCODE id and a summary dict. If the endpoint shape differs, adjust `resolve_gencode`/`summarize_transcripts` field names to the live JSON (cache the raw response and inspect it). If GTEx proves unusable, `collect` already degrades to skips — note it in limitations.

- [ ] **Step 6: Commit**

```bash
git add analysis/collect/gtex_prostate.py tests/test_gtex_prostate.py
git commit -m "feat(collect): GTEx normal-prostate isoform baseline

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 7: Fold new sources into the master table + GTEx context

**Files:**
- Modify: `analysis/integrate/master_table.py`
- Test: `tests/test_master_table.py` (add cases)

**Interfaces:**
- Produces (added to `master_table.py`):
  - `normalize_literature(row: dict, feats=None) -> dict`
  - `normalize_ascancer(row: dict) -> dict`
  - `attach_gtex(rows: list[dict], baseline: dict[str,dict]) -> list[dict]` — fill `gtex_baseline` column from a gene→summary map.
  - `build_master` extended to read `results_collected/raw/*.tsv` (when present) and call `attach_gtex` using `results_collected/gtex_prostate_baseline.tsv` (when present).

- [ ] **Step 1: Write the failing tests**

```python
# append to tests/test_master_table.py
from integrate.master_table import normalize_literature, attach_gtex


def test_normalize_literature_row_shape():
    row = {"gene": "AR", "uniprot": "P10275", "variant": "AR-V7",
           "event_type": "cryptic exon (CE3)", "effect": "LBD-truncated",
           "pca_note": "resistance", "pmids": "25", "query": "q"}
    out = normalize_literature(row)
    assert out["source"] == "Literature"
    assert out["gene"] == "AR" and out["uniprot"] == "P10275"
    assert "AR-V7" in out["source_id"]
    assert "25" in out["provenance"]


def test_attach_gtex_fills_column():
    rows = [{"gene": "AR", "gtex_baseline": ""}]
    base = {"AR": {"max_median": 9.0, "second_median": 4.0,
                   "multi_isoform_normal": True}}
    out = attach_gtex(rows, base)
    assert "multi-isoform" in out[0]["gtex_baseline"]
```

- [ ] **Step 2: Run tests, verify the new ones fail**

Run: `pytest tests/test_master_table.py -v`
Expected: 2 new tests FAIL (ImportError on the new names)

- [ ] **Step 3: Implement the additions in `master_table.py`**

Add near the other normalizers:
```python
def normalize_literature(row: dict, feats=None) -> dict:
    out = _blank()
    out.update(
        variant_id=f"{row['gene']}|Literature|{row['variant']}",
        gene=row["gene"], uniprot=row.get("uniprot", ""), source="Literature",
        source_id=row["variant"], event_type=row.get("event_type", ""),
        status="literature-reported",
        change_class=row.get("effect", ""),
        pca_evidence=row.get("pca_note", ""),
        provenance=f"PMID:{row.get('pmids','')}".rstrip(":"),
    )
    return out


def normalize_ascancer(row: dict) -> dict:
    out = _blank()
    out.update(
        variant_id=f"{row['gene']}|ASCancerAtlas|{row.get('as_id','')}",
        gene=row["gene"], source="ASCancerAtlas",
        source_id=row.get("as_id", ""), event_type=row.get("event_type", ""),
        status=row.get("cancer", ""), pca_evidence="ASCancerAtlas PRAD",
        provenance="ASCancerAtlas",
    )
    return out


def attach_gtex(rows: list[dict], baseline: dict) -> list[dict]:
    for r in rows:
        b = baseline.get(r["gene"])
        if not b:
            continue
        tag = "multi-isoform" if b.get("multi_isoform_normal") else "single-isoform"
        r["gtex_baseline"] = (f"{tag} normal (top {b.get('max_median')}, "
                              f"2nd {b.get('second_median')} TPM)")
    return rows
```

Extend `build_master` before `return add_corroboration(rows)`:
```python
    # Literature source (optional intermediate)
    lit_path = "results_collected/raw/literature_variants.tsv"
    if os.path.exists(lit_path):
        for r in pd.read_csv(lit_path, sep="\t").to_dict("records"):
            feats = feats_for(r.get("uniprot", ""))
            rows.append(normalize_literature(r, feats))

    # ASCancerAtlas source (optional intermediate)
    asc_path = "results_collected/raw/ascanceratlas_prad.tsv"
    if os.path.exists(asc_path):
        df_asc = pd.read_csv(asc_path, sep="\t")
        for r in df_asc.to_dict("records"):
            rows.append(normalize_ascancer(r))

    rows = add_corroboration(rows)

    # GTEx context (optional)
    gtex_path = "results_collected/gtex_prostate_baseline.tsv"
    if os.path.exists(gtex_path):
        gb = {r["gene"]: r for r in pd.read_csv(gtex_path, sep="\t").to_dict("records")}
        rows = attach_gtex(rows, gb)
    return rows
```
And remove the now-duplicated final `return add_corroboration(rows)` (corroboration is computed once, after all sources are appended). Confirm `build_master` returns `rows` exactly once.

- [ ] **Step 4: Run tests, verify they pass**

Run: `pytest tests/test_master_table.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Rebuild the full master table**

Run: `python analysis/integrate/master_table.py`
Then GTEx (needs the master to exist for its gene list):
Run: `python analysis/collect/gtex_prostate.py`
Then rebuild once more so GTEx context attaches:
Run: `python analysis/integrate/master_table.py`
Expected: source breakdown now includes UniProt, TCGA, Literature (+ ASCancerAtlas if available); some rows have `gtex_baseline` filled.

- [ ] **Step 6: Commit**

```bash
git add analysis/integrate/master_table.py tests/test_master_table.py
git commit -m "feat(integrate): fold literature/ASCancerAtlas/GTEx into master table

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 8: Integrated report (Word + Markdown) + run everything end-to-end

**Files:**
- Create: `analysis/build_integrated_report.py`
- Create: `results_collected/README.md`
- Modify: `.gitignore`
- Test: `tests/test_integrated_report.py`

**Interfaces:**
- Consumes: `results_collected/master_variants.tsv`, `gtex_prostate_baseline.tsv`.
- Produces:
  - `corroborated(df) -> df` — rows where `corroborating_sources >= 2`, sorted desc.
  - `build(master_tsv, out_docx, out_md) -> dict` — writes both reports; returns `{rows, genes, multi_source_genes}`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_integrated_report.py
import pandas as pd
from build_integrated_report import corroborated


def test_corroborated_filters_and_sorts():
    df = pd.DataFrame([
        {"gene": "AR", "corroborating_sources": 3},
        {"gene": "X", "corroborating_sources": 1},
        {"gene": "ITGA6", "corroborating_sources": 2},
    ])
    out = corroborated(df)
    assert list(out["gene"]) == ["AR", "ITGA6"]
```

- [ ] **Step 2: Run test, verify it fails**

Run: `pytest tests/test_integrated_report.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Implement `build_integrated_report.py`**

```python
"""Integrated multi-source PCa alternative-splicing report (Word + Markdown).

Three parts mirror the prior reports (sequences -> function impact -> PCa/
resistance) but now span all sources, with a cross-source corroboration section.
"""
from __future__ import annotations
from datetime import date

import pandas as pd
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH

MASTER = "results_collected/master_variants.tsv"
OUT_DOCX = "docs/Integrated_PCa_splicing_report.docx"
OUT_MD = "docs/Integrated_PCa_splicing_report.md"


def corroborated(df: pd.DataFrame) -> pd.DataFrame:
    d = df[df["corroborating_sources"].astype(int) >= 2].copy()
    return d.sort_values("corroborating_sources", ascending=False)


def _src_counts(df):
    return df["source"].value_counts().to_dict()


def build(master_tsv=MASTER, out_docx=OUT_DOCX, out_md=OUT_MD) -> dict:
    df = pd.read_csv(master_tsv, sep="\t").fillna("")
    df["corroborating_sources"] = df["corroborating_sources"].astype(int)
    multi = corroborated(df)
    multi_genes = sorted(multi["gene"].unique())
    counts = _src_counts(df)

    # ---------- Markdown ----------
    L = []
    L.append(f"# Integrated Prostate-Cancer Alternative-Splicing Report\n")
    L.append(f"_Multi-source: UniProt + TCGA SpliceSeq + literature + "
             f"ASCancerAtlas + GTEx baseline — {date.today():%Y-%m-%d}_\n")
    L.append("## Sources\n")
    for s, n in counts.items():
        L.append(f"- **{s}**: {n} variant rows")
    L.append(f"\nTotal: {len(df)} variant rows across {df['gene'].nunique()} genes.\n")
    L.append("## Cross-source corroboration (>=2 sources)\n")
    L.append("Genes supported by more than one independent source — the "
             "highest-confidence leads.\n")
    L.append("| Gene | Sources | Variants | Best PCa evidence |")
    L.append("|------|--------|----------|-------------------|")
    for g in multi_genes:
        sub = df[df["gene"] == g]
        srcs = ",".join(sorted(sub["source"].unique()))
        ev = next((x for x in sub["pca_evidence"] if x), "")
        L.append(f"| {g} | {srcs} | {len(sub)} | {ev} |")
    L.append("\n## Functional-impact highlights (domain / NMD / localization)\n")
    fi = df[(df["domains_hit"] != "") | (df["nmd_flag"] == "NMD-candidate")]
    for _, r in fi.sort_values("corroborating_sources", ascending=False).head(40).iterrows():
        L.append(f"- **{r['gene']}** ({r['source']}, {r['source_id']}): "
                 f"{r['change_class']}; domains hit: {r['domains_hit'] or 'n.a.'}; "
                 f"NMD: {r['nmd_flag']}; loc: {r['loc_flags'] or 'n.a.'}; "
                 f"GTEx: {r['gtex_baseline'] or 'n.a.'}")
    L.append("\n## Limitations\n")
    for s in [
        "Disease/resistance, NMD, and domain-disruption calls are predictions, "
        "not measurements.",
        "Cross-source corroboration is gene-level, not exon/variant-level.",
        "UniProt isoforms and literature variants are curated presence/absence; "
        "GTEx gives normal context only; cohort HRs are univariate single-cohort.",
        "CancerSplicingQTL was unreachable (HTTP 403) and is omitted.",
    ]:
        L.append(f"- {s}")
    open(out_md, "w").write("\n".join(L) + "\n")

    # ---------- Word ----------
    doc = Document()
    doc.styles["Normal"].font.name = "Calibri"
    doc.styles["Normal"].font.size = Pt(10.5)
    doc.add_heading("Integrated Prostate-Cancer Alternative-Splicing Report", 0)
    sub = doc.add_paragraph()
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sub.add_run(f"UniProt + TCGA + literature + ASCancerAtlas + GTEx — "
                f"{date.today():%Y-%m-%d}").italic = True

    doc.add_heading("Sources", 1)
    for s, n in counts.items():
        doc.add_paragraph(f"{s}: {n} variant rows", style="List Bullet")
    doc.add_paragraph(f"Total {len(df)} variant rows across "
                      f"{df['gene'].nunique()} genes.")

    doc.add_heading("Part 1 — Cross-source corroboration", 1)
    doc.add_paragraph("Genes supported by >=2 independent sources (highest-"
                      "confidence leads).")
    t = doc.add_table(rows=1, cols=4); t.style = "Light Grid Accent 1"
    for c, h in zip(t.rows[0].cells, ["Gene", "Sources", "Variants", "PCa evidence"]):
        c.paragraphs[0].add_run(h).bold = True
    for g in multi_genes:
        sub_df = df[df["gene"] == g]
        ev = next((x for x in sub_df["pca_evidence"] if x), "")
        cells = t.add_row().cells
        for c, v in zip(cells, [g, ",".join(sorted(sub_df["source"].unique())),
                                str(len(sub_df)), ev]):
            c.paragraphs[0].add_run(str(v)).font.size = Pt(8)

    doc.add_heading("Part 2 — Predicted functional impact", 1)
    doc.add_paragraph("Domain disruption, NMD candidacy, and localization "
                      "changes (UniProt-feature-driven; predictions).")
    fi2 = df[(df["domains_hit"] != "") | (df["nmd_flag"] == "NMD-candidate")]
    for _, r in fi2.sort_values("corroborating_sources", ascending=False).head(40).iterrows():
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(f"{r['gene']} ").bold = True
        p.add_run(f"({r['source']}, {r['source_id']}): {r['change_class']}; "
                  f"domains: {r['domains_hit'] or 'n.a.'}; NMD: {r['nmd_flag']}; "
                  f"loc: {r['loc_flags'] or 'n.a.'}").font.size = Pt(9)

    doc.add_heading("Part 3 — Implications for PCa & therapeutic resistance", 1)
    doc.add_paragraph("AR-V7 (LBD-truncated, constitutively active) is the "
                      "headline resistance variant; PTEN/AR-axis and the "
                      "adhesion/invasion genes corroborated across sources are the "
                      "priority leads. All links are hypotheses.")

    doc.add_heading("Limitations", 1)
    for s in [
        "Predictions (NMD, domain disruption, disease links), not measurements.",
        "Corroboration is gene-level; exon-level matching is future work.",
        "GTEx is normal context; cohort HRs are univariate single-cohort.",
        "CancerSplicingQTL omitted (HTTP 403).",
    ]:
        doc.add_paragraph(s, style="List Bullet")

    doc.save(out_docx)
    return {"rows": len(df), "genes": int(df["gene"].nunique()),
            "multi_source_genes": len(multi_genes)}


def main() -> None:
    stats = build()
    print(stats)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run unit test, verify it passes**

Run: `pytest tests/test_integrated_report.py -v`
Expected: PASS

- [ ] **Step 5: Add `.gitignore` rule + README**

Append to `.gitignore`:
```
results_collected/cache/
results_collected/raw/
```
Create `results_collected/README.md`:
```markdown
# Collected multi-source PCa splicing

Generated by analysis/collect/*, analysis/annotate/*, analysis/integrate/master_table.py.
Reproduce:
1. python analysis/collect/literature_variants.py
2. python analysis/collect/db_ascanceratlas.py
3. python analysis/integrate/master_table.py
4. python analysis/collect/gtex_prostate.py
5. python analysis/integrate/master_table.py   # re-run to attach GTEx context
6. python analysis/build_integrated_report.py

Tracked outputs: master_variants.tsv, gtex_prostate_baseline.tsv.
cache/ and raw/ are git-ignored (re-fetchable). --offline replays cache.
```

- [ ] **Step 6: Full end-to-end run**

Run, in order:
```bash
python analysis/collect/literature_variants.py
python analysis/collect/db_ascanceratlas.py
python analysis/integrate/master_table.py
python analysis/collect/gtex_prostate.py
python analysis/integrate/master_table.py
python analysis/build_integrated_report.py
```
Expected: final command prints `{'rows': ..., 'genes': ..., 'multi_source_genes': ...}` and writes both report files.

- [ ] **Step 7: Full test suite green**

Run: `pytest -q`
Expected: all tests pass (existing + new).

- [ ] **Step 8: Commit outputs**

```bash
git add analysis/build_integrated_report.py tests/test_integrated_report.py .gitignore results_collected/README.md results_collected/master_variants.tsv results_collected/gtex_prostate_baseline.tsv docs/Integrated_PCa_splicing_report.docx docs/Integrated_PCa_splicing_report.md
git commit -m "feat: integrated multi-source PCa splicing report + master table

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:** literature source (Task 4) ✓; ASCancerAtlas (Task 5) ✓; GTEx baseline (Task 6) ✓; CancerSplicingQTL dropped-with-note (Global Constraints + report limitations) ✓; domain mapping / NMD / disorder-localization (Task 2) ✓; UniProt-feature backbone (Task 1) ✓; master schema + corroboration (Tasks 3, 7) ✓; unified table + one integrated report (Task 8) ✓; caching + `--offline` + graceful degradation (all collectors) ✓; tests for pure functions (every task) ✓; incremental build order (Tasks 1-2 annotator first, 3 integrator, 4-6 collectors, 7-8 fold+report) ✓.

**Type consistency:** `Feature` namedtuple, `fetch_features`, `annotate_variant`, `changed_interval`, `COLUMNS`, `build_master`, `normalize_*`, `attach_gtex`, `summarize_transcripts`, `parse_table`, `parse_europepmc`, `corroborated` — names used consistently across producing and consuming tasks.

**Placeholder note:** `db_ascanceratlas.py` intentionally has a runtime-discovered `DOWNLOAD_URL`/columns (Task 5 Step 1) — this is a genuine discovery step, not a plan placeholder; the parser is fully specified and tested against a real-sample fixture, and the collector degrades gracefully if the source is unusable.

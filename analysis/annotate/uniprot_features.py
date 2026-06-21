"""Fetch and parse UniProt feature tracks used for functional-impact analysis."""
from __future__ import annotations
import json
import os
from collections import namedtuple

import requests

Feature = namedtuple("Feature", "type start end description")

# UniProt feature types we use, grouped by downstream purpose.
# (Exact type strings as they appear in UniProtKB JSON.)
DOMAIN_TYPES = {"Domain", "Region", "Repeat", "Zinc finger", "DNA binding",
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

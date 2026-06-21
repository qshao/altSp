"""ASCancerAtlas prostate alternative-splicing events.

Source: ASCancerAtlas curated-knowledge API (discovered from the Vue SPA at
https://ngdc.cncb.ac.cn/ascancer/Download):
    GET https://ngdc.cncb.ac.cn/ascancer/api/browse/knowledge?page=N&size=200
returns paged JSON `{"totalItems": ..., "items": [...]}`. In the knowledge
table prostate events are labelled cancer_name == "Prostate Cancer".

Parser is tested offline against tests/fixtures/ascancer_sample.json.
Graceful: any network/format failure -> [] + warning, never raises.
"""
from __future__ import annotations
import csv
import json
import math
import os
import time

import requests

API = "https://ngdc.cncb.ac.cn/ascancer/api/browse/knowledge"
HEADERS = {"User-Agent": "Mozilla/5.0 (research data collection)"}
PROSTATE_NAMES = {"Prostate Cancer", "Prostate Adenocarcinoma"}
PAGE_SIZE = 200


def parse_items(items: list[dict]) -> list[dict]:
    rows = []
    for r in items:
        cancer = (r.get("cancer_name") or "").strip()
        if cancer not in PROSTATE_NAMES and "prostate" not in cancer.lower():
            continue
        rows.append({
            "gene": (r.get("gene_name") or "").strip(),
            "event_type": (r.get("as_type") or "").strip(),
            "cancer": cancer,
            "as_id": (r.get("event_id") or r.get("as_model_id") or "").strip(),
            "sv_ensembl_id": (r.get("sv_ensembl_id") or "").strip(),
            "note": "ASCancerAtlas knowledge",
        })
    return rows


def _fetch_all(cache_dir: str, offline: bool) -> list[dict]:
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, "knowledge_items.json")
    if os.path.exists(path):
        return json.load(open(path))
    if offline:
        return []
    # page 1 to learn totalItems, then page through with per-page retries.
    first = requests.get(API, params={"page": 1, "size": PAGE_SIZE},
                         headers=HEADERS, timeout=60)
    first.raise_for_status()
    payload = first.json()
    total = int(payload.get("totalItems", 0))
    items = list(payload.get("items", []))
    pages = math.ceil(total / PAGE_SIZE) if total else 1
    for page in range(2, pages + 1):
        for _ in range(4):
            try:
                r = requests.get(API, params={"page": page, "size": PAGE_SIZE},
                                 headers=HEADERS, timeout=60)
                r.raise_for_status()
                items.extend(r.json().get("items", []))
                break
            except Exception:
                time.sleep(1.5)
        time.sleep(0.2)
    json.dump(items, open(path, "w"))
    return items


def collect(cache_dir: str, offline: bool = False) -> list[dict]:
    try:
        return parse_items(_fetch_all(cache_dir, offline))
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
    cols = ["gene", "event_type", "cancer", "as_id", "sv_ensembl_id", "note"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"ASCancerAtlas prostate rows: {len(rows)} "
          f"({len({r['gene'] for r in rows})} genes) -> {out}")


if __name__ == "__main__":
    main()

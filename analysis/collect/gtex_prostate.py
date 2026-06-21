"""GTEx normal-prostate isoform baseline (median transcript expression)."""
from __future__ import annotations
import csv
import json
import os

import requests

API = "https://gtexportal.org/api/v2"
TISSUE = "Prostate"
MULTI_ISOFORM_TPM = 1.0
HEADERS = {"User-Agent": "Mozilla/5.0 (research data collection)",
           "Accept": "application/json"}


def summarize_transcripts(payload: dict) -> dict:
    data = payload.get("data", [])
    meds = sorted((d.get("median", 0.0) for d in data), reverse=True)
    top = meds[0] if meds else 0.0
    second = meds[1] if len(meds) > 1 else 0.0
    gencode = payload.get("gencodeId", "") or (data[0].get("gencodeId", "") if data else "")
    return {
        "gencodeId": gencode,
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
    r = requests.get(url, params=params, headers=HEADERS, timeout=60)
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
                          {"gencodeId": gencode_id, "tissueSiteDetailId": TISSUE,
                           "datasetId": "gtex_v8"},
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

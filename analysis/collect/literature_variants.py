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

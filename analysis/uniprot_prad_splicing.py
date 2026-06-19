"""Search UniProt for prostate-cancer-related human proteins that undergo
alternative splicing, and save the results.

Strategy (UniProtKB REST query):
  reviewed (Swiss-Prot) AND human (9606)
  AND keyword "Alternative splicing" (KW-0025)
  AND free-text "prostate cancer"

For each hit we record the curated isoform catalogue (cc_alternative_products),
the splice-variant sequence features (ft_var_seq), disease involvement, and
function, plus a FASTA of every annotated isoform.

Usage:
  python analysis/uniprot_prad_splicing.py
  python analysis/uniprot_prad_splicing.py --query '<custom UniProt query>'
"""
from __future__ import annotations
import argparse
import json
import os
import time

import requests

API = "https://rest.uniprot.org/uniprotkb/search"
OUTDIR = "uniprot_prad_splicing"

DEFAULT_QUERY = (
    '(organism_id:9606) AND (reviewed:true) AND (keyword:KW-0025) '
    'AND ("prostate cancer")'
)

FIELDS = ",".join([
    "accession", "id", "protein_name", "gene_names", "organism_name",
    "length", "keyword", "cc_alternative_products", "ft_var_seq",
    "cc_disease", "cc_function", "sequence",
])


def paged(query: str, fmt: str, fields: str | None, size: int = 500):
    """Yield response objects across all cursor pages."""
    params = {"query": query, "format": fmt, "size": size}
    if fields:
        params["fields"] = fields
    url = API
    while url:
        r = requests.get(url, params=params if url == API else None, timeout=60)
        r.raise_for_status()
        yield r
        url = r.links.get("next", {}).get("url")
        params = None
        time.sleep(0.2)


def n_isoforms(rec: dict) -> int:
    for c in rec.get("comments", []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            return len(c.get("isoforms", []))
    return 0


def isoform_names(rec: dict) -> list[str]:
    out = []
    for c in rec.get("comments", []):
        if c.get("commentType") == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms", []):
                ids = ",".join(iso.get("isoformIds", []))
                name = (iso.get("name") or {}).get("value", "")
                out.append(f"{ids}:{name}" if name else ids)
    return out


def n_splice_features(rec: dict) -> int:
    return sum(1 for f in rec.get("features", [])
               if f.get("type") in ("Alternative sequence", "VAR_SEQ"))


def disease_terms(rec: dict) -> list[str]:
    out = []
    for c in rec.get("comments", []):
        if c.get("commentType") == "DISEASE":
            d = c.get("disease", {})
            name = d.get("diseaseId") or d.get("diseaseAccession")
            if name:
                out.append(name)
    return out


def first_function(rec: dict) -> str:
    for c in rec.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                return texts[0].get("value", "").replace("\t", " ")
    return ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--outdir", default=OUTDIR)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # 1) full JSON records (paged)
    records: list[dict] = []
    for r in paged(args.query, "json", FIELDS):
        records.extend(r.json().get("results", []))
    print(f"UniProt hits: {len(records)}")

    with open(os.path.join(args.outdir, "proteins.json"), "w") as f:
        json.dump(records, f, indent=2)

    # 2) flat TSV summary
    cols = ["accession", "gene", "protein_name", "length", "n_isoforms",
            "n_splice_variants", "is_prostate_disease", "disease_terms",
            "isoforms", "function"]
    with open(os.path.join(args.outdir, "proteins.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for rec in records:
            gene = ""
            genes = rec.get("genes", [])
            if genes:
                gene = (genes[0].get("geneName") or {}).get("value", "")
            pname = ((rec.get("proteinDescription", {})
                      .get("recommendedName", {})
                      .get("fullName", {}) or {}).get("value", ""))
            dts = disease_terms(rec)
            is_pca = any("prostate" in d.lower() for d in dts)
            row = [
                rec.get("primaryAccession", ""), gene, pname,
                str((rec.get("sequence", {}) or {}).get("length", "")),
                str(n_isoforms(rec)), str(n_splice_features(rec)),
                "yes" if is_pca else "no", "; ".join(dts) or ".",
                " | ".join(isoform_names(rec)) or ".",
                first_function(rec)[:500],
            ]
            f.write("\t".join(row) + "\n")

    # 3) isoform FASTA (all annotated isoforms).
    # NOTE: includeIsoform only returns isoforms that *also* match the query;
    # a keyword/organism query excludes them (isoforms carry no such fields).
    # So we re-fetch by accession in batches, which keeps the isoforms.
    accs = [r.get("primaryAccession", "") for r in records]
    accs = [a for a in accs if a]
    stream = "https://rest.uniprot.org/uniprotkb/stream"
    iso_text = []
    for i in range(0, len(accs), 100):
        batch = accs[i:i + 100]
        q = " OR ".join(f"accession:{a}" for a in batch)
        rr = requests.get(stream, params={
            "query": q, "format": "fasta", "includeIsoform": "true"},
            timeout=120)
        rr.raise_for_status()
        iso_text.append(rr.text)
        time.sleep(0.2)
    with open(os.path.join(args.outdir, "isoforms.fasta"), "w") as f:
        f.write("".join(iso_text))

    # 4) provenance
    with open(os.path.join(args.outdir, "query.txt"), "w") as f:
        f.write(f"endpoint: {API}\nquery: {args.query}\nfields: {FIELDS}\n"
                f"records: {len(records)}\n")

    n_pca = sum(1 for rec in records
                if any("prostate" in d.lower() for d in disease_terms(rec)))
    n_iso = sum(1 for rec in records if n_isoforms(rec) > 1)
    print(f"  with curated prostate-cancer disease term: {n_pca}")
    print(f"  with >1 annotated isoform:                 {n_iso}")
    print(f"  isoform FASTA sequences: "
          f"{open(os.path.join(args.outdir,'isoforms.fasta')).read().count('>')}")
    print(f"written to {args.outdir}/ (proteins.tsv, proteins.json, "
          f"isoforms.fasta, query.txt)")


if __name__ == "__main__":
    main()

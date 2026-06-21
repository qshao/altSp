"""Normalize every splicing source into one master variant table."""
from __future__ import annotations
import json
import os
from collections import defaultdict

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
        provenance=f"ASCancerAtlas {row.get('sv_ensembl_id','')}".strip(),
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


def add_corroboration(rows: list[dict]) -> list[dict]:
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

    # Literature source (optional intermediate)
    lit_path = "results_collected/raw/literature_variants.tsv"
    if os.path.exists(lit_path):
        for r in pd.read_csv(lit_path, sep="\t").to_dict("records"):
            feats = feats_for(r.get("uniprot", ""))
            rows.append(normalize_literature(r, feats))

    # ASCancerAtlas source (optional intermediate)
    asc_path = "results_collected/raw/ascanceratlas_prad.tsv"
    if os.path.exists(asc_path):
        for r in pd.read_csv(asc_path, sep="\t").to_dict("records"):
            rows.append(normalize_ascancer(r))

    rows = add_corroboration(rows)

    # GTEx normal-prostate context (optional)
    gtex_path = "results_collected/gtex_prostate_baseline.tsv"
    if os.path.exists(gtex_path):
        gb = {r["gene"]: r for r in pd.read_csv(gtex_path, sep="\t").to_dict("records")}
        rows = attach_gtex(rows, gb)
    return rows


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

"""Cross-reference UniProt prostate-cancer alternative-splicing genes against the
cohort-derived prioritized events.

Joins:
  uniprot_prad_splicing/proteins.tsv   (curated PCa splicing proteins)
  analysis/prioritized_events.tsv      (TCGA PRAD progression-associated events)

Output: analysis/crossref_uniprot_events.tsv  + console summary.
A gene in the overlap is one UniProt independently annotates as a
prostate-cancer-relevant alternatively-spliced protein AND that shows a
protein-altering, outcome-associated splice in this TCGA cohort.
"""
from __future__ import annotations
import pandas as pd

UNI = "uniprot_prad_splicing/proteins.tsv"
PRI = "analysis/prioritized_events.tsv"
OUT = "analysis/crossref_uniprot_events.tsv"


def main() -> None:
    uni = pd.read_csv(UNI, sep="\t")
    pri = pd.read_csv(PRI, sep="\t")

    # Keep the most significant event per gene on the cohort side.
    pri_best = (pri.sort_values("bh_p")
                .drop_duplicates("Gene", keep="first")
                .rename(columns={"Gene": "gene"}))

    uni_cols = uni[["gene", "accession", "n_isoforms", "n_splice_variants",
                    "is_prostate_disease", "disease_terms", "function"]]

    merged = pri_best.merge(uni_cols, on="gene", how="inner")
    merged = merged.sort_values(["is_prostate_disease", "bh_p"],
                                ascending=[False, True])

    out = merged[[
        "gene", "accession", "Splice_Event", "Type", "endpoint", "HR",
        "direction", "bh_p", "protein_change", "n_isoforms",
        "is_prostate_disease", "disease_terms", "function",
    ]]
    out.to_csv(OUT, sep="\t", index=False)

    print(f"UniProt PCa-splicing genes:     {uni['gene'].nunique()}")
    print(f"prioritized-event genes:        {pri['gene'].nunique() if 'gene' in pri else pri['Gene'].nunique()}")
    print(f"overlap (cross-referenced):     {len(out)}")
    print(f"  with curated PCa disease term: {(out['is_prostate_disease']=='yes').sum()}")
    print()
    show = out.copy()
    show["function"] = show["function"].str.slice(0, 70)
    show["disease_terms"] = show["disease_terms"].str.slice(0, 40)
    with pd.option_context("display.width", 240, "display.max_colwidth", 72):
        print(show.to_string(index=False))


if __name__ == "__main__":
    main()

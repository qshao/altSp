"""Link protein-changing splicing events to PRAD outcome (DFI/PFI hazard ratios).

The core pipeline filters events only on FDR/dPSI. The SpliceSeq CSV, however,
also carries per-event Cox hazard ratios against two clinical endpoints that are
the natural proxies for aggressive / therapy-resistant disease:

  DFI = Disease-Free Interval   (recurrence)
  PFI = Progression-Free Interval (progression -> proxy for treatment failure)

For each endpoint SpliceSeq reports two PSI dichotomizations:
  med* = split at median PSI
  fit* = split at best-fitting PSI cutpoint
with HR (high-PSI vs low-PSI group), a raw p-value, and a BH-adjusted p-value.

HR > 1  -> high splice-in PSI associates with WORSE outcome (event is risk-up)
HR < 1  -> high splice-in PSI associates with BETTER outcome (event is risk-down)

This script joins those stats onto the protein-changing events and ranks the
ones with a significant, protein-altering, outcome association.
"""
from __future__ import annotations
import math
import pandas as pd

CSV = "spliceseq_info_PRAD.csv"
EVENTS = "results/events_proteins.tsv"
OUT = "analysis/prioritized_events.tsv"

# Genes with established roles in prostate cancer / castration & therapy resistance.
PCA_GENES = {
    "AR", "FOXA1", "TP53", "PTEN", "RB1", "SPOP", "ERG", "ETV1", "ETV4",
    "KLK3", "NCOA1", "NCOA2", "NCOA3", "CDK12", "BRCA1", "BRCA2", "ATM",
    "CHD1", "MYC", "AKT1", "PIK3CA", "ZBTB16", "NDRG1", "FOLH1", "STEAP1",
    "STEAP2", "KLK2", "HOXB13", "CTNNB1", "APC", "GNAS", "MED12", "ZFHX3",
    "NKX3-1", "CCND1", "CDKN1B", "AURKA", "BIRC5", "GHR", "SRRM4", "HNRNPA1",
    "PTBP1", "ESRP1", "ESRP2", "GLYATL1", "EZH2", "WNT5A", "DNMT3B",
}


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main() -> None:
    df = pd.read_csv(CSV, low_memory=False)
    ev = pd.read_csv(EVENTS, sep="\t")

    # Restrict to events that actually change the protein product.
    changing = ev[~ev["comparison"].isin(["identical", "single_isoform", "unresolved"])]

    keep = [
        "Splice_Event", "HRmedDFI", "bhHRmedDFI", "HRfitDFI", "bhHRfitDFI",
        "HRmedPFI", "bhHRmedPFI", "HRfitPFI", "bhHRfitPFI",
        "nObsHRfitPFI", "nObsHRfitDFI",
    ]
    surv = df[keep].copy()
    m = changing.merge(surv, on="Splice_Event", how="left")

    rows = []
    for _, r in m.iterrows():
        # Best (smallest) BH-adjusted p across the four endpoint/cut combos,
        # carrying along the matching HR.
        combos = [
            ("DFI", f(r["bhHRmedDFI"]), f(r["HRmedDFI"])),
            ("DFI", f(r["bhHRfitDFI"]), f(r["HRfitDFI"])),
            ("PFI", f(r["bhHRmedPFI"]), f(r["HRmedPFI"])),
            ("PFI", f(r["bhHRfitPFI"]), f(r["HRfitPFI"])),
        ]
        combos = [(e, p, hr) for (e, p, hr) in combos
                  if p is not None and hr is not None and hr > 0]
        if not combos:
            continue
        endpoint, bestp, besthr = min(combos, key=lambda t: t[1])
        rows.append({
            "Splice_Event": r["Splice_Event"],
            "Gene": r["Gene_Symbol"],
            "Type": r["Splice_Type"],
            "dPSI": r["PSI_Difference"],
            "protein_change": r["comparison"],
            "endpoint": endpoint,
            "HR": round(besthr, 3),
            "direction": "risk-up" if besthr > 1 else "risk-down",
            "absLogHR": round(abs(math.log(besthr)), 3),
            "bh_p": bestp,
            "pca_gene": r["Gene_Symbol"] in PCA_GENES,
        })

    out = pd.DataFrame(rows)
    sig = out[out["bh_p"] < 0.05].copy()
    sig = sig.sort_values(["bh_p", "absLogHR"], ascending=[True, False])
    sig.to_csv(OUT, sep="\t", index=False)

    print(f"protein-changing events:            {len(changing)}")
    print(f"  with any survival HR:             {len(out)}")
    print(f"  BH-significant (p<0.05) outcome:  {len(sig)}")
    print(f"    risk-up (high PSI worse):       {(sig['direction']=='risk-up').sum()}")
    print(f"    risk-down:                      {(sig['direction']=='risk-down').sum()}")
    print(f"  involving known PCa genes:        {sig['pca_gene'].sum()}")
    print()
    print("=== Top 25 protein-changing, outcome-associated events ===")
    cols = ["Gene", "Type", "endpoint", "HR", "direction", "bh_p",
            "protein_change", "pca_gene"]
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(sig[cols].head(25).to_string(index=False))
    print()
    print("=== Known prostate-cancer / resistance genes (all sig hits) ===")
    pg = sig[sig["pca_gene"]]
    print(pg[cols].to_string(index=False) if len(pg) else "  (none)")


if __name__ == "__main__":
    main()

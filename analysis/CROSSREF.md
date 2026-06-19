# Cross-reference: UniProt PCa-splicing genes ∩ prioritized TCGA events

Joins `uniprot_prad_splicing/proteins.tsv` (238 UniProt genes that are human,
reviewed, alternatively-spliced, and prostate-cancer-associated) with
`analysis/prioritized_events.tsv` (315 genes with a protein-altering,
progression-associated splice in TCGA PRAD). Reproduce:
`python analysis/crossref_uniprot_events.py` → `analysis/crossref_uniprot_events.tsv`.

## Overlap: 7 genes

All 7 surface on the **PFI** (progression) endpoint and all are
`frameshift_or_seqchange` events — i.e. UniProt independently flags them as
PCa-relevant spliced proteins, and the cohort shows a protein-changing splice
that tracks with progression.

| Gene | Event | HR (PFI) | Direction | BH p | UniProt isoforms | Note |
|------|-------|---------|-----------|------|------------------|------|
| **PRUNE2** | ES | 0.54 | risk-down | 0.046 | 5 | Established prostate tumour suppressor (regulated by the PCA3 lncRNA); protective-isoform direction is consistent. |
| **ITGA6** | ES | 2.47 | risk-up | 0.0023 | 8 | Integrin α6 (laminin receptor); drives adhesion/invasion and metastasis. |
| **CTNND1** | ES | 0.45 | risk-down | 0.0034 | 32 | p120-catenin; master regulator of cadherin stability and cell–cell adhesion. |
| **FES** | AD | 3.59 | risk-up | 0.024 | 4 | Tyrosine kinase downstream of surface receptors; pro-proliferative. |
| **STEAP3** | ME | 0.46 | risk-down | 0.039 | 4 | Metalloreductase / TSAP6; STEAP family is a prostate antigen/therapeutic-target class. |
| **SMN1** | AA | 1.69 | risk-up | 0.045 | 4 | SMN complex — assembles the spliceosomal snRNPs (splicing-machinery feedback). |
| **PRICKLE4** | AA | 4.35 | risk-up | 0.0006 | 3 | Planar-cell-polarity component; directed migration. |

## Interpretation

- **Directionality is biologically coherent.** The two adhesion/tumour-suppressor
  genes (PRUNE2, CTNND1) and STEAP3 are **risk-down** (high splice-in → better
  PFI = protective isoform), whereas the invasion/proliferation genes (ITGA6,
  FES, PRICKLE4) are **risk-up**.
- **The curated hereditary-PCa genes do *not* appear in the overlap.** AR, PTEN,
  ELAC2, RNASEL, MSR1, CHEK2 etc. are in the UniProt set but their TCGA splice
  events had small ΔPSI / no resolved protein pair and so never reached the
  prioritized list — the same annotation-coverage limit noted for AR earlier.
- **Two strongest leads:** PRICKLE4 (lowest p, large HR) and ITGA6 (well-known
  invasion driver with 8 curated isoforms) — both warrant isoform-level
  follow-up.

These remain hypotheses; the cohort HRs are univariate and single-cohort.

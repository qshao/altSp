# How PRAD alternative splicing relates to disease and therapy resistance

Generated from `spliceseq_info_PRAD.csv` (TCGA SpliceSeq, prostate adenocarcinoma)
joined to the protein-level results in `results/`. Reproduce with
`python analysis/survival_link.py`; ranked table in `analysis/prioritized_events.tsv`.

## What we measured

The core pipeline keeps events that are differentially spliced tumor-vs-normal
(FDR < 0.05, |ΔPSI| ≥ 0.1) and asks whether the splice changes the protein. By
itself that says nothing about *resistance*. The SpliceSeq table also carries,
per event, Cox hazard ratios against two clinical endpoints — the usable proxies
for aggressive / treatment-resistant disease:

- **PFI** — progression-free interval (progression ≈ treatment failure)
- **DFI** — disease-free interval (recurrence)

`HR > 1` ⇒ high splice-in PSI tracks with **worse** outcome ("risk-up").
We restricted to events that **both** change the protein **and** carry a
BH-adjusted survival p < 0.05.

## Headline numbers

- 839 protein-changing significant events; **371 (44%) associate with outcome**
  at BH p < 0.05.
- **240 risk-up vs 131 risk-down** — splicing dysregulation is biased toward
  worse-prognosis directions.
- Signal concentrates on **PFI (261) over DFI (110)** — i.e. on *progression*,
  the endpoint nearest to therapy resistance.

## Mechanistic patterns (the interesting part)

**Intron retention is almost unidirectionally bad.** Retained-intron (RI) events
split **42 risk-up : 1 risk-down**, and essentially all are
`frameshift_or_seqchange`. Retaining an intron injects a premature stop / frame
disruption → truncated protein or NMD-mediated loss. In this cohort that loss
consistently tracks with faster progression. Alt-acceptor/donor (AA/AD) events
are similarly skewed risk-up (33:2, 28:5). Exon-skipping (ES) is balanced
(133:117) — it can go either way depending on what the cassette exon encodes.

**The protein really changes.** 344/371 hits are frameshift/sequence changes,
not tidy in-frame length tweaks — these splices remodel the protein, not just
trim it.

## Prostate-cancer / resistance genes flagged

Directly in the protein-changing significant set:

| Gene | Event | HR (PFI) | Protein effect | Note |
|------|-------|---------|----------------|------|
| **KLK2** | ES | 3.08 ↑ | frameshift/seqchange | kallikrein, AR-target, PSA-family |
| **SPOP** | ME | 4.26 ↑ | +257 aa | SPOP is a top PRAD driver (E3 substrate adaptor) |
| **SPOP** | ES | 2.14 ↑ | +112 aa | second SPOP isoform shift |

**AR (androgen receptor) — the canonical resistance gene — is present but subtle
here.** Its events have tiny ΔPSI (~0.01) and Ensembl GRCh37 r75 does not give a
clean before/after protein pair, so they fall out of the protein-change filter.
But the **alternate-terminator** events do carry outcome signal:
`AR_AT_89349` (splice-in of AR-004/AR-005, a C-terminally altered isoform)
associates with **better** PFI (HR 0.02–0.16, BH p ≈ 0.003–0.007). Alternate AR
C-termini are the structural class that, in the extreme (AR-V7), drop the
ligand-binding domain and drive enzalutamide/abiraterone resistance — worth
following up against an annotation that includes AR-V transcripts.

## Caveats

- HRs are univariate, single-cohort, not multiplicity-controlled across all four
  endpoint/cut combinations jointly — treat as hypothesis generation.
- "risk-up" = association, not causation; PSI cutpoints (median/fit) are coarse.
- r75 annotation misses several therapy-resistance splice variants (AR-V7, etc.).

## Suggested next steps

1. Re-run resolution against an annotation with AR-V / known oncogenic variants.
2. For top RI hits, confirm PTC/NMD prediction on the retained-intron sequence.
3. Cross-check the top genes (KLK2, SPOP, plus SMARCD3/RBM6/OGG1 DNA-repair &
   chromatin hits) against AR-signaling and DNA-repair pathways implicated in
   castration resistance.

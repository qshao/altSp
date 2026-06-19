# UniProt isoform sequence-impact (companion to the TCGA before/after analysis)

Applies the same three-part workflow used for the TCGA SpliceSeq events
(before/after sequences → functional impact → PCa/therapy-resistance) to the
curated UniProt prostate-cancer splicing set.

## Pipeline

1. `analysis/uniprot_sequence_impact.py` — pairs each protein's **canonical**
   sequence ("before") with every annotated **isoform** ("after") from
   `uniprot_prad_splicing/isoforms.fasta`, computes the common prefix/suffix,
   classifies the change, and attaches UniProt VAR_SEQ / FUNCTION / DISEASE.
   → `analysis/uniprot_sequence_impact.json` (238 proteins, 632 isoform pairs).
2. `analysis/uniprot_report_notes.py` — curated PCa/resistance notes + the TCGA
   cross-reference HRs for the 7 overlap genes.
3. `analysis/build_uniprot_report.py` — builds
   `docs/UniProt_PCa_splicing_report.docx` (3 parts + synthesis + 238-row
   appendix). Run: `PYTHONPATH=analysis python analysis/build_uniprot_report.py`.

## Change-class distribution (632 canonical→isoform pairs)

| Class | Count |
|-------|------:|
| substituted segment | 379 |
| internal deletion (in-frame) | 203 |
| internal insertion (in-frame) | 44 |
| C-terminal truncation | 6 |

## Why this complements the TCGA arm

- **It captures AR-V7.** UniProt isoform 3 of AR (P10275-3, 644 aa) is the
  ligand-binding-domain-truncated, constitutively active variant that drives
  enzalutamide/abiraterone resistance — the *exact* event GRCh37 r75 hid from
  the cohort analysis. UniProt fills the TCGA arm's headline blind spot.
- **12 curated PCa-disease genes** are featured with full before/after sequences:
  AR, PTEN, KLF6, CHEK2, ELAC2(HPC2), RNASEL(HPC1), MSR1, MSMB, EHBP1, HNF1B,
  MXI1, EPHB2.
- **7 genes are corroborated by both datasets** (UniProt PCa-splicing annotation
  AND a progression-associated protein-changing TCGA splice): PRICKLE4, ITGA6,
  CTNND1, FES, STEAP3, SMN1, PRUNE2 — with internally consistent direction.
- KLF6 illustrates the principle directly: its oncogenic splice isoform
  (KLF6-SV1) flips a tumour suppressor into a driver of aggressive disease.

All disease/resistance links are hypothesis-level; UniProt isoform catalogues
are curated presence/absence, not tumour-specific expression.

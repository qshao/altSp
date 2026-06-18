# Alternative Splicing → Protein Sequences (PRAD)

Resolves the "before" (spliced-out) and "after" (spliced-in) protein sequences
for significant alternative splicing events in `spliceseq_info_PRAD.csv`
(TCGA SpliceSeq). Isoform names are resolved against Ensembl GRCh37 release-75.

## Run

    python -m altsplice_protein spliceseq_info_PRAD.csv

First run downloads ~50 MB of Ensembl GRCh37 r75 files into `data/`.
Outputs are written to `results/`:

- `proteins.fasta` — one protein per resolved isoform.
- `events_proteins.tsv` — per-event before/after proteins and comparison.
- `unresolved.txt` — isoform names with no protein (non-coding or unmapped).

Options: `--fdr-max`, `--dpsi-min`, `--results-dir`, `--data-dir`, `--no-rest`.

## Develop

    PYTHONPATH=src python -m pytest

# Alternative Splicing → Protein Sequences (PRAD) — Design

**Date:** 2026-06-18
**Input:** `spliceseq_info_PRAD.csv` (TCGA SpliceSeq, prostate adenocarcinoma; 67,534 events)

## Goal

For each *significant* alternative splicing event, obtain the protein sequence
"before" and "after" the splicing change — i.e. the protein of the
spliced-out (event-excluded) transcript vs. the spliced-in (event-included)
transcript — so the protein-level consequence of the event can be compared.

## Key facts about the data

- Protein-relevant columns are **SpliceSeq isoform names** at transcript level:
  - `SpliceIn_IsoName` — isoform(s) that **include** the event (the "after" / inclusion form).
  - `SpliceOut_IsoName` — isoform(s) that **exclude** the event (the "before" / exclusion form).
  - `AltExons_IsoName` — informational (exon-tagged), not used for sequence retrieval.
- A side may list multiple comma-separated isoforms (e.g. `A2M-001,A2M-002`).
- Naming uses the legacy Havana `GENE-001` style. Verified resolvable against the
  **Ensembl GRCh37 (release 75)** annotation — the build TCGA SpliceSeq was made
  from. Confirmed: `A2M-001` → `ENST00000318602` (gene `ENSG00000175899`).
- Splice types: ES, AP, AT, AA, AD, RI, ME.
  - ES, RI, AA, AD, ME events typically name isoforms on **both** sides → complete
    before/after pair.
  - **AP** (alternate promoter) and **AT** (alternate terminator) name a single
    isoform by nature → no single-transcript before/after pair.

## Significant subset

Filter: `FDR_Difference < 0.05` AND `|PSI_Difference| >= 0.1`.

- Significant events: **2,678**
- Of those with **both** SpliceIn and SpliceOut isoforms (complete pair): **914**
  (ES 634, AD 92, AA 85, RI 80, ME 23)
- Single-isoform significant events (AP 973, AT 297, plus partials): processed too,
  flagged as not having a complete pair.
- Unique isoform names to resolve across the subset: **5,845**

## Source / resolution

**Primary — local Ensembl GRCh37 release 75 files** (downloaded once, pinned for provenance):
- `Homo_sapiens.GRCh37.75.gtf.gz` — provides `transcript_name → transcript_id (ENST)
  → protein_id (ENSP)` and `transcript_biotype`.
- `Homo_sapiens.GRCh37.75.pep.all.fa.gz` — protein sequences keyed by ENSP/ENST.

Build an in-memory map `transcript_name → {ENST, ENSP, biotype, protein_seq}` and
resolve all 5,845 names locally (instant, no rate limits, reproducible).

**Fallback — Ensembl GRCh37 REST** for any name not found locally:
- `GET grch37.rest.ensembl.org/xrefs/symbol/homo_sapiens/{name}` → transcript ID.
- `GET grch37.rest.ensembl.org/sequence/id/{ENST}?type=protein` → protein sequence.

Non-coding isoforms (no protein product) are recorded explicitly rather than dropped.

## Pipeline

1. **Filter** the CSV to the significant subset; extract per-event SpliceIn / SpliceOut
   isoform name lists and the unique union of names.
2. **Build resolver map** from the local GRCh37 r75 GTF + pep FASTA.
3. **Resolve** each unique name → protein. Classify each as: resolved-protein,
   non-coding (resolved transcript, no protein), or unresolved (no transcript).
   REST-fallback the unresolved set.
4. **Assemble outputs** per event with before/after comparison.

## Outputs (written to a `results/` directory)

- **`proteins.fasta`** — deduplicated protein sequence per isoform.
  Header: `>{isoform_name}|{ENST}|{ENSP}|{gene_symbol}`.
- **`events_proteins.tsv`** — one row per significant event:
  `Splice_Event, Gene_Symbol, Splice_Type, PSI_Difference, FDR_Difference,
   SpliceIn_isoforms, SpliceIn_protein_ids, SpliceOut_isoforms, SpliceOut_protein_ids,
   has_complete_pair, before_len, after_len, comparison`
  where `comparison ∈ {identical, length_change(Δaa), frameshift_likely,
   noncoding_side, single_isoform, unresolved}`.
- **`unresolved.txt`** — isoform names with no transcript or no protein, with reason.

For complete-pair events: "before" = SpliceOut protein(s), "after" = SpliceIn
protein(s). Where a side lists multiple isoforms, all are emitted; the comparison
uses the longest protein on each side as the representative.

## Validation

- Spot-check 3 known names (`A2M-001`, `A2M-013`, `A4GALT-003`) resolve to ENST with
  protein and match REST output.
- Report resolution coverage (% of 5,845 names resolved locally, % via REST, %
  unresolved/non-coding) before declaring done.

## Out of scope (YAGNI)

- De-novo translation of custom SpliceSeq exon graphs (we use annotated transcript
  proteins, not reconstructed reading frames).
- Events failing the significance filter.
- AP/AT cross-event pairing (single-isoform events emit their one protein only).

## Notes

- Working directory is **not** a git repository, so the design doc is written but
  not committed. (Offer `git init` if version control is wanted.)

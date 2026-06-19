# Alternative Splicing in Prostate Adenocarcinoma: Methods and Conclusions

**Question.** How might alternative splicing affect prostate cancer (PRAD) and
its resistance to therapy?

**Scope.** This is a *hypothesis-generating* re-analysis of an existing TCGA
SpliceSeq table for prostate adenocarcinoma. It does not prove causation; it
ranks splicing events that are simultaneously (a) differentially spliced in
tumor, (b) protein-altering, and (c) associated with clinical progression.

---

## 1. Data and tooling

| Input | Description |
|-------|-------------|
| `spliceseq_info_PRAD.csv` | TCGA SpliceSeq PRAD table — **67,534 events** across **12,436 genes**, with tumor/normal PSI, differential-splicing statistics, and Cox survival statistics per event. |
| Ensembl GRCh37 release-75 | GTF + peptide FASTA, used to map SpliceSeq isoform names → transcripts → protein sequences. |
| `src/altsplice_protein/` | Existing pipeline that filters significant events and resolves before/after protein products. |
| `analysis/survival_link.py` | New analysis joining protein-changing events to the survival statistics. |

Event-type composition of the raw table: ES 32,974 · AP 12,003 · AT 9,195 ·
AA 5,074 · AD 4,456 · RI 3,457 · ME 375
(ES exon-skip, AP/AT alternate promoter/terminator, AA/AD alternate
acceptor/donor, RI retained intron, ME mutually-exclusive exons).

---

## 2. What the existing pipeline does

For each event the pipeline:

1. **Filters for significance** — keeps events with FDR < 0.05 and
   |ΔPSI| ≥ 0.1 (`filtering.py`).
2. **Resolves isoforms to proteins** — maps each isoform name to an Ensembl
   transcript and its peptide, with a REST fallback (`resolver.py`,
   `ensembl_data.py`).
3. **Compares before vs after** — classifies the splice-out ("before") against
   the splice-in ("after") protein as `identical`, an in-frame
   `length_change(±N aa)`, a `frameshift_or_seqchange`, `noncoding_side`, or
   `single_isoform` (`comparison.py`).

Result on this cohort (`results/events_proteins.tsv`, 2,678 significant events):

| comparison class | n |
|------------------|---|
| single_isoform (no paired side) | 1,764 |
| **frameshift_or_seqchange** | 775 |
| identical | 73 |
| noncoding_side | 35 |
| in-frame length_change | 31 |

**Gap addressed here:** this output establishes *that the protein changes* but
says nothing about clinical impact or therapy resistance.

---

## 3. What this analysis adds

The SpliceSeq table also carries, per event, univariate Cox **hazard ratios**
against two clinical endpoints — the usable proxies for aggressive /
treatment-resistant disease:

- **PFI** — progression-free interval (progression ≈ treatment failure)
- **DFI** — disease-free interval (recurrence)

For each endpoint, two PSI dichotomizations are provided (`med` = split at
median PSI, `fit` = best-fitting cutpoint), each with an HR and a
Benjamini–Hochberg-adjusted p-value.

Interpretation: **HR > 1 ⇒ high splice-in PSI tracks with worse outcome
("risk-up"); HR < 1 ⇒ "risk-down".**

`analysis/survival_link.py`:

1. Takes the protein-changing events (drops `identical`, `single_isoform`,
   `unresolved`) — **839 events**.
2. Joins the four endpoint/cut HRs from the CSV.
3. For each event keeps the most significant endpoint/cut combination.
4. Retains events with BH-adjusted p < 0.05 and ranks by p then |log HR|.
5. Flags genes from a curated PRAD / castration-resistance gene set.

Output: `analysis/prioritized_events.tsv` (ranked) and console summary.

---

## 4. Results

### 4.1 Outcome association is widespread and directionally biased

- Of 839 protein-changing events, **371 (44%) associate with outcome** at
  BH p < 0.05.
- **240 risk-up vs 131 risk-down** — dysregulated splicing skews toward
  worse-prognosis directions.
- Signal concentrates on **PFI (261) over DFI (110)** — i.e. on *progression*,
  the endpoint nearest to therapy resistance.

### 4.2 Intron retention is near-unidirectionally harmful

Splice-type × direction among the 371 hits:

| Type | risk-down | risk-up |
|------|-----------|---------|
| RI (retained intron) | 1 | **42** |
| AA (alt acceptor) | 2 | 33 |
| AD (alt donor) | 5 | 28 |
| ES (exon skip) | 117 | 133 |
| ME (mutually-excl.) | 6 | 4 |

Retained introns and alternative splice sites are strongly skewed risk-up, and
essentially all are `frameshift_or_seqchange`. Mechanistically, retaining an
intron injects a premature stop / frame disruption → truncated protein or
NMD-mediated loss of function; in this cohort that loss consistently tracks with
faster progression. Exon-skipping is balanced — it can help or harm depending on
what the cassette exon encodes.

### 4.3 The protein genuinely changes

344 of 371 hits are frameshift/sequence changes rather than small in-frame
length tweaks — these splices remodel the protein product, not merely trim it.

### 4.4 Top protein-changing, progression-associated events

| Gene | Type | Endpoint | HR | Direction | BH p | Protein effect |
|------|------|----------|----|-----------|------|----------------|
| SEC31A | ES | PFI | 0.32 | risk-down | 4.7e-05 | frameshift/seqchange |
| QTRT1 | RI | PFI | 5.02 | risk-up | 8.4e-05 | frameshift/seqchange |
| WASH4P | RI | PFI | 3.39 | risk-up | 8.9e-05 | length −414 aa |
| DOCK7 | ES | PFI | 0.29 | risk-down | 9.8e-05 | frameshift/seqchange |
| ADHFE1 | AA | PFI | 4.51 | risk-up | 1.4e-04 | frameshift/seqchange |
| SGSM3 | RI | PFI | 4.70 | risk-up | 1.5e-04 | noncoding side |
| VCL | ES | PFI | 0.29 | risk-down | 1.5e-04 | frameshift/seqchange |
| SMARCD3 | RI | DFI | 5.87 | risk-up | 1.6e-04 | frameshift/seqchange |
| FADS3 | RI | PFI | 3.95 | risk-up | 2.2e-04 | frameshift/seqchange |

(Full ranked list in `analysis/prioritized_events.tsv`.) Several hits are
biologically suggestive: **SMARCD3** (SWI/SNF chromatin remodeling), **OGG1**
(base-excision DNA repair), **VCL** (vinculin, adhesion/invasion).

### 4.5 Known prostate-cancer / resistance genes

Directly in the protein-changing significant set:

| Gene | Event | HR (PFI) | Protein effect | Note |
|------|-------|---------|----------------|------|
| **KLK2** | ES | 3.08 ↑ | frameshift/seqchange | kallikrein, AR target (PSA family) |
| **SPOP** | ME | 4.26 ↑ | +257 aa | top PRAD driver (E3 substrate adaptor) |
| **SPOP** | ES | 2.14 ↑ | +112 aa | second SPOP isoform shift |

**AR (androgen receptor)** — the canonical resistance gene — is subtle in *this*
annotation: its events have tiny ΔPSI (~0.01) and GRCh37 r75 gives no clean
before/after protein pair, so they drop out of the protein-change filter. But
the **alternate-terminator** events still carry outcome signal —
`AR_AT_89349` (splice-in of AR-004/AR-005, a C-terminally altered isoform)
associates with **better** PFI (HR 0.02–0.16, BH p ≈ 0.003–0.007). Alternate AR
C-termini are the structural class that, at the extreme (AR-V7), drop the
ligand-binding domain and drive enzalutamide/abiraterone resistance. This is the
single most important follow-up: re-run against an annotation that includes
AR-V transcripts.

---

## 5. Conclusions

1. **Alternative splicing in PRAD is broadly coupled to disease progression.**
   Nearly half of protein-altering significant events associate with PFI/DFI,
   and the bias is toward *worse* outcome.

2. **Loss-of-function splicing (intron retention, alt splice sites) is the
   dominant harmful mode.** RI events are ~42:1 risk-up and almost entirely
   frameshifting — a coherent mechanism (truncation/NMD) rather than scattered
   noise.

3. **Progression (PFI), not just recurrence (DFI), is where the signal lives** —
   consistent with splicing contributing to the aggressive, treatment-refractory
   phenotype.

4. **Established PRAD genes surface with the expected direction** (KLK2, SPOP
   risk-up), and **AR alternate-terminus splicing carries outcome signal** even
   though full AR-V variants are invisible to GRCh37 r75 — a concrete,
   resistance-relevant lead.

These are **hypotheses for experimental follow-up**, not established drivers.

---

## 6. Limitations

- HRs are **univariate, single-cohort**, and not jointly multiplicity-controlled
  across the four endpoint/cut combinations; "best-of-four" selection is
  optimistic. Treat p-values as ranking aids, not confirmatory.
- Association ≠ causation; median/fit PSI cutpoints are coarse.
- GRCh37 r75 **misses key resistance splice variants** (AR-V7 and others), so
  the most therapy-relevant events are systematically under-detected here.
- Protein comparison uses the longest representative isoform per side; complex
  multi-isoform events are simplified.

---

## 7. Recommended next steps

1. **Re-annotate against AR-V / known oncogenic splice variants** and re-run the
   AR (and other resistance-gene) events.
2. **Confirm NMD/PTC predictions** for the top retained-intron hits (QTRT1,
   SGSM3, FADS3, SMARCD3) on the retained-intron sequence.
3. **Multivariate / multi-cohort validation** of the top events (adjust for
   stage, Gleason; replicate in an independent PRAD cohort).
4. **Pathway focus:** cross-check chromatin (SMARCD3), DNA-repair (OGG1) and
   AR-signaling (KLK2, AR) hits against mechanisms implicated in
   castration-resistant prostate cancer.

---

## 8. Reproducibility

```bash
# protein resolution (writes results/)
python -m altsplice_protein spliceseq_info_PRAD.csv

# survival linkage + prioritization (writes analysis/prioritized_events.tsv)
python analysis/survival_link.py
```

Artifacts: `results/events_proteins.tsv`, `analysis/prioritized_events.tsv`,
`analysis/FINDINGS.md` (short summary), this report.

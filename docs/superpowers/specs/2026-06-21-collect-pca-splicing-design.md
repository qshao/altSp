# Design: Collect & integrate more prostate-cancer alternative splicing, with functional-impact analysis

Date: 2026-06-21
Status: approved (brainstorming)

## Goal

Expand the prostate-cancer (PCa) alternative-splicing evidence base beyond the two
sources already in the repo — TCGA SpliceSeq PRAD cohort events
(`analysis/prioritized_events.tsv`) and curated UniProt isoforms
(`analysis/uniprot_sequence_impact.json`) — by adding three new sources, then
analyze how each splice variant changes protein function and what that implies
for prostate cancer and therapeutic resistance. Deliver one normalized master
variant table plus a single integrated report.

## Scope

### New data sources (all three requested)
1. **Literature-curated variants** — a curated seed list of well-characterized PCa
   splice variants (e.g. AR-V7, AR-V567es, KLF6-SV1, cyclin D1b/CCND1, FGFR2
   IIIb/IIIc, BCL2L1/Bcl-xS, etc.), each with PMID provenance fetched from the
   **EuropePMC REST API** (reachable, ~21k hits for the topic). This is curated
   table + provenance lookup, NOT full-text NLP mining.
2. **More splicing databases** — **ASCancerAtlas** PRAD alternative-splicing
   events via its `Download` bulk files (reachable). **CancerSplicingQTL is
   dropped** to a documented best-effort note (returns HTTP 403, bot-blocked).
3. **Normal-tissue baseline** — **GTEx v2 API** median transcript expression for
   prostate tissue, to flag genes whose alternative isoforms are also abundant in
   normal prostate (less likely to be tumor-specific). Versioned GENCODE IDs must
   be resolved via the GTEx gene endpoint (the API is param-sensitive).

### Functional-impact analysis (all three requested), UniProt-feature-driven
- **Domain mapping** — overlap the changed residue interval with UniProt
  Domain / Region / Binding-site / Active-site features → name the disrupted domain.
- **NMD / truncation** — protein-level heuristic (C-terminal loss not within the
  last ~50 aa → NMD-candidate), explicitly labeled an approximation (true NMD
  needs transcript exon structure).
- **Disorder & localization** — overlap with Signal-peptide / Transmembrane /
  Topological-domain / Disordered-region features.

### Deliverable
Unified master variant table + one integrated report (the chosen option).

## Non-goals (YAGNI)
- Full-text literature NLP mining (use curated seed + provenance instead).
- CancerSplicingQTL scraping around its bot block.
- Per-variant GTEx PSI (use gene/isoform-level normal context, not junction PSI).
- Quantitative tumor-specific expression modeling.

## Architecture

Annotation backbone = **UniProt feature tracks**: every variant from every source
resolves to a UniProt accession; functional impact = overlap of the changed
residue interval with that accession's annotated features.

```
analysis/collect/
  literature_variants.py    # curated PCa variant seed list + EuropePMC provenance -> intermediate TSV
  db_ascanceratlas.py       # ASCancerAtlas PRAD events (bulk download, cached)    -> intermediate TSV
  gtex_prostate.py          # GTEx v2 median transcript expression, prostate        -> baseline TSV
analysis/annotate/
  uniprot_features.py       # fetch + cache domain/region/signal/TM/disorder tracks per accession
  functional_impact.py      # changed-interval x feature overlap -> domains_hit, NMD, localization, disorder
analysis/integrate/
  master_table.py           # normalize ALL sources (incl. existing TCGA + UniProt) -> master_variants.tsv
analysis/build_integrated_report.py   # one Word + Markdown report
```

Each unit has one purpose and a file-based interface (TSV/JSON), so it can be
built and tested independently.

## Data flow

1. Each collector writes its own intermediate TSV under `results_collected/raw/`.
2. `uniprot_features.py` collects every accession seen across all sources and
   fetches+caches feature tracks (`results_collected/cache/uniprot_features/`).
3. `functional_impact.py` computes the changed interval and annotates each variant
   that has before/after sequences; for variants described only qualitatively
   (some literature/DB rows) it maps to a domain by name where possible and marks
   the rest as `sequence_not_resolved`.
4. `master_table.py` normalizes everything into one schema and computes a
   `corroborating_sources` count per gene/variant.
5. `build_integrated_report.py` emits the report.

### Changed-interval definition
For variants with before/after sequences (UniProt isoforms, TCGA peptides):
`changed_interval = [common_prefix+1 .. len - common_suffix]`. For variants
without sequences, record the described event qualitatively.

## Master schema (one row per variant x source)

| column | meaning |
|--------|---------|
| variant_id | stable id (gene + source + source_id) |
| gene | HGNC symbol |
| uniprot | primary accession |
| source | TCGA / UniProt / Literature / ASCancerAtlas |
| source_id | event id / isoform id / PMID etc. |
| event_type | ES/AA/AD/AP/AT/RI/ME/isoform/described |
| status | ΔPSI, HR, or qualitative status |
| before_len / after_len | aa lengths when resolved |
| change_class | truncation / substituted segment / in-frame indel / unresolved |
| changed_interval | residue range, when resolved |
| domains_hit | UniProt domains overlapping the change |
| regions_lost | named regions/motifs lost |
| loc_flags | signal / TM / topological overlap |
| disorder_overlap | fraction/flag of change in disordered region |
| nmd_flag | NMD-candidate (heuristic) / no / n.a. |
| pca_evidence | disease term / cohort HR / PMID |
| gtex_baseline | normal-prostate isoform context for the gene |
| corroborating_sources | count of independent sources for this gene/variant |
| provenance | PMID / URL / file |

Outputs: `results_collected/{master_variants.tsv, provenance.tsv,
gtex_prostate_baseline.tsv}` and `docs/Integrated_PCa_splicing_report.{docx,md}`.
The report keeps the three-part structure (sequences -> function impact ->
PCa/resistance) but multi-source, adding a cross-source corroboration section.

## Error handling
- All network calls cached to disk; never re-hit an API for cached data.
- `--offline` flag replays cache/fixtures.
- Each collector is independent; integrator consumes whatever intermediates exist.
- Missing/blocked source -> warning + provenance note, never a crash
  (CancerSplicingQTL block handled this way).

## Testing
Unit tests in existing `tests/` for the pure functions, with synthetic fixtures:
- changed-interval computation
- domain/feature overlap
- NMD heuristic
- schema normalization / corroboration count
Collectors run against cached fixtures in `--offline` mode.

## Build order (incremental)
1. `uniprot_features.py` + `functional_impact.py` + tests — immediately enriches
   the existing TCGA + UniProt variants (value before any new source).
2. `master_table.py` over existing two sources + report skeleton.
3. Collectors one at a time: literature -> ASCancerAtlas -> GTEx, folding each
   into the master table and report.

## Honesty constraints (carried from prior work)
- All disease/resistance links remain hypothesis-level.
- NMD and domain-disruption calls are predictions, not measurements.
- UniProt isoform catalogues are curated presence/absence, not tumor expression.
- GTEx provides normal context only.

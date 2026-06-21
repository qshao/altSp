"""Comprehensive PDF report for the multi-source PCa alternative-splicing analysis.

Pipeline: master_variants.tsv (+ GTEx baseline) -> matplotlib figures -> styled
HTML -> WeasyPrint PDF. Includes a curated "Top 5 splicing cases" section.

All disease / NMD / domain-disruption statements are predictions or
literature-curated annotations, never measurements from this pipeline.
"""
from __future__ import annotations
import html as _html
import json
import math
import os
import re
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd
from weasyprint import HTML

MASTER = "results_collected/master_variants.tsv"
FIGDIR = "results_collected/figures"
OUT_PDF = "docs/PCa_splicing_comprehensive_report.pdf"
UNI_IMPACT = "analysis/uniprot_sequence_impact.json"
TCGA_IMPACT = "analysis/sequence_impact.json"
FASTA = "results_collected/featured_sequences.fasta"

# Representative before/after variant for each featured gene (Sections 5 & 6).
# (section, gene, source, variant_id, note). CCND1 is literature-only with no
# sequence in the pipeline's sources.
SEQ_REPS = [
    ("Top 5", "AR", "UniProt", "P10275-3", "AR-V7 (LBD lost)"),
    ("Top 5", "ITGA6", "UniProt", "P23229-2", "β-propeller isoform"),
    ("Top 5", "CTNND1", "UniProt", "O60716-2", "p120-1AB"),
    ("Top 5", "KLF6", "UniProt", "Q99612-3", "truncated isoform"),
    ("Top 5", "CCND1", "Literature", "", "cyclin D1b — not in pipeline sources"),
    ("Novel", "MAP3K7", "UniProt", "O43318-4", "NMD-candidate truncation"),
    ("Novel", "EXOC7", "TCGA", "EXOC7_ES_43570", "DFI HR 5.456 event"),
    ("Novel", "PRUNE2", "UniProt", "Q8WUY3-4", "C-terminal isoform"),
    ("Novel", "ENAH", "TCGA", "ENAH_ES_9989", "DFI HR 3.135 event"),
    ("Novel", "STEAP3", "UniProt", "Q658P3-4", "reductase isoform"),
]

# Brand palette
C_AR = "#b1283a"      # AR / resistance red
C_BLUE = "#23527c"
C_TEAL = "#2a9d8f"
C_GOLD = "#e9a13b"
C_GREY = "#6c757d"
SRC_COLORS = {"UniProt": C_BLUE, "TCGA": C_AR, "ASCancerAtlas": C_TEAL,
              "Literature": C_GOLD}

plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150,
                     "font.family": "DejaVu Sans"})


# ----------------------------------------------------------------------------
# Top-5 curated cases (biology + cross-source corroboration driven)
# ----------------------------------------------------------------------------
TOP5 = [
    {
        "rank": 1, "gene": "AR", "variant": "AR-V7 / AR-V567es",
        "uniprot": "P10275", "sources": "UniProt + Literature + ASCancerAtlas",
        "event": "Cryptic-exon inclusion / exon skipping truncating the "
                 "ligand-binding domain (LBD).",
        "impact": "Loss of the C-terminal NR ligand-binding domain (pipeline "
                  "changed-interval overlaps “NR LBD; Nuclear receptor”; "
                  "flagged NMD-candidate). The DNA-binding domain and "
                  "transactivation region are retained → a constitutively "
                  "active, ligand-independent receptor.",
        "pca": "The premier driver of castration-resistant prostate cancer "
               "(CRPC). AR-V7 is a clinically validated biomarker of resistance "
               "to enzalutamide and abiraterone; AR-V567es drives ligand-"
               "independent signalling in CRPC.",
        "why": "Highest-confidence case: corroborated across three independent "
               "sources, the pipeline reproduces the known LBD-loss mechanism, "
               "and it is the best-established splicing-driven resistance axis "
               "in the disease.",
    },
    {
        "rank": 2, "gene": "ITGA6", "variant": "Integrin α6 isoforms (P23229-2/-4/-7)",
        "uniprot": "P23229", "sources": "UniProt + TCGA + ASCancerAtlas",
        "event": "Exon skipping / substituted segments across the extracellular "
                 "β-propeller.",
        "impact": "Changed intervals repeatedly overlap FG-GAP propeller repeats "
                  "(FG-GAP 1–7) and ligand-interaction regions; several "
                  "isoforms are flagged NMD-candidate. Predicted disruption of "
                  "laminin engagement and integrin signalling.",
        "pca": "Integrin α6 governs tumour-cell adhesion, migration and "
               "invasion; isoform switching is linked to a more invasive "
               "phenotype. Corroborated by an adverse cohort progression signal.",
        "why": "Three-source corroboration plus a clear, repeated domain-"
               "disruption + NMD signal on an adhesion/invasion gene.",
    },
    {
        "rank": 3, "gene": "CTNND1", "variant": "p120-catenin isoforms (incl. p120-1A/-1AB)",
        "uniprot": "O60716", "sources": "UniProt + TCGA + ASCancerAtlas",
        "event": "Alternative N-terminal start / internal exon usage across the "
                 "armadillo (ARM) repeat region.",
        "impact": "Isoforms remodel the ARM-repeat cadherin-binding surface and "
                  "nuclear-localization signals (pipeline shows ARM-repeat and "
                  "NLS overlap). Shifts the balance between cadherin "
                  "stabilisation and Rho-GTPase / transcriptional signalling.",
        "pca": "p120-catenin isoform switching (toward the N-terminally longer "
               "isoform 1) accompanies loss of epithelial adhesion and is "
               "associated with invasive prostate carcinoma; corroborated by a "
               "cohort progression signal.",
        "why": "Three-source corroboration on a core adhesion regulator with a "
               "concrete domain-level mechanism.",
    },
    {
        "rank": 4, "gene": "KLF6", "variant": "KLF6-SV1",
        "uniprot": "Q99612", "sources": "Literature + UniProt",
        "event": "Alternative 5′ splice-site use producing a truncated, "
                 "zinc-finger–deficient isoform.",
        "impact": "Loss of the DNA-binding C2H2 zinc-finger domain converts a "
                  "tumour-suppressor transcription factor into a dominant-"
                  "negative product that antagonises wild-type KLF6.",
        "pca": "KLF6-SV1 over-expression is associated with aggressive, "
               "treatment-refractory prostate cancer and poorer outcome — a "
               "splicing event that inactivates a tumour suppressor.",
        "why": "Literature-curated with UniProt support; a clean example of "
               "splicing-driven tumour-suppressor inactivation distinct from the "
               "AR axis.",
    },
    {
        "rank": 5, "gene": "CCND1", "variant": "Cyclin D1b",
        "uniprot": "P24385", "sources": "Literature + UniProt",
        "event": "Intron-4 retention producing an alternative C-terminus.",
        "impact": "Loss of Thr286, the residue required for nuclear export and "
                  "degradation → nuclear retention of a stabilised cyclin D1 "
                  "that also gains AR-coactivator activity.",
        "pca": "Cyclin D1b enhances AR transcriptional output and proliferation; "
               "a G870A germline polymorphism biasing its production is linked to "
               "prostate-cancer risk and progression.",
        "why": "Mechanistically crisp single-residue consequence that ties the "
               "cell-cycle machinery directly back to AR signalling.",
    },
]


# ----------------------------------------------------------------------------
# Novel-but-credible cases: gene is credibly PCa-relevant, but the *splicing*
# angle is under-studied in prostate cancer. signal = cohort HR / NMD / domain.
# ----------------------------------------------------------------------------
NOVEL5 = [
    {
        "rank": 1, "gene": "MAP3K7", "variant": "TAK1 NMD-candidate isoforms",
        "uniprot": "O43318", "sources": "UniProt + ASCancerAtlas",
        "signal": "2 NMD-candidate isoforms; kinase-domain region affected",
        "event": "Truncating / frame-disrupting splicing predicted to trigger "
                 "nonsense-mediated decay.",
        "known": "MAP3K7 (TAK1) <b>deletion</b> is among the best-characterised "
                 "events in aggressive PCa (~30–40% of tumours; CHD1 co-deletion "
                 "drives AR target genes, AR-V7, and enzalutamide resistance).",
        "novel": "That body of work is about copy-number loss — not splicing. A "
                 "splicing route to MAP3K7 loss-of-function would phenocopy the "
                 "deletion, and our NMD-candidate isoforms point exactly there.",
    },
    {
        "rank": 2, "gene": "EXOC7", "variant": "Exo70 epithelial/mesenchymal switch",
        "uniprot": "Q9UPT5", "sources": "TCGA + ASCancerAtlas",
        "signal": "Strongest cohort signal in the table (DFI HR ≈ 5.5, risk-up)",
        "event": "Exocyst Exo70 isoform switching (alternative exon usage).",
        "known": "Exo70 epithelial↔mesenchymal isoform switching is a proven, "
                 "ESRP1-controlled invasion driver in breast and colon cancer "
                 "(mesenchymal isoform recruits Arp2/3 for actin remodelling).",
        "novel": "Never tested in prostate, yet it carries the single strongest "
                 "prognostic signal across the whole cohort and is corroborated "
                 "by two sources — a high-value, untested mechanism.",
    },
    {
        "rank": 3, "gene": "PRUNE2", "variant": "PRUNE2 isoforms (PCA3 locus)",
        "uniprot": "Q8WUY3", "sources": "TCGA + UniProt",
        "signal": "4 domain-overlapping isoforms; protective PFI HR ≈ 0.54",
        "event": "Alternative splicing producing multiple PRUNE2 transcripts.",
        "known": "PRUNE2 is the validated tumour suppressor at the PCA3 locus — "
                 "PCA3 being the clinically used prostate-cancer biomarker lncRNA "
                 "that down-regulates PRUNE2.",
        "novel": "Maximally prostate-specific credibility, but the individual "
                 "PRUNE2 isoforms and their suppressor activity remain "
                 "uncharacterised — a focused splicing question on a marquee "
                 "prostate gene.",
    },
    {
        "rank": 4, "gene": "ENAH", "variant": "Mena invasion isoforms (MenaINV / 11a)",
        "uniprot": "Q8N8S7", "sources": "TCGA + ASCancerAtlas",
        "signal": "Corroborated; DFI HR ≈ 3.1 (risk-up)",
        "event": "Inclusion of invasion-associated exons (exon 4 / 11a).",
        "known": "MenaINV and Mena11a are textbook invasion / EMT splicing "
                 "switches in breast and lung carcinoma.",
        "novel": "The prostate literature is essentially empty, yet the same "
                 "actin-regulatory invasion machinery is corroborated here with "
                 "an adverse cohort signal.",
    },
    {
        "rank": 5, "gene": "STEAP3", "variant": "STEAP3 metalloreductase isoforms",
        "uniprot": "Q658P3", "sources": "TCGA + UniProt",
        "signal": "Domain-overlapping isoforms; protective PFI HR ≈ 0.46",
        "event": "Alternative splicing of the STEAP3 transmembrane reductase.",
        "known": "STEAP3 belongs to the six-transmembrane epithelial-antigen-of-"
                 "prostate family (STEAP1/2 are PCa therapy targets); it is a "
                 "p53-inducible metalloreductase up-regulated across cancers.",
        "novel": "Its splicing is uncharacterised in prostate cancer despite the "
                 "family's standing as prostate antigens — an under-explored but "
                 "well-grounded lead.",
    },
]


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def fig_sources(df, path):
    counts = df["source"].value_counts()
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    bars = ax.bar(counts.index, counts.values,
                  color=[SRC_COLORS.get(s, C_GREY) for s in counts.index])
    ax.bar_label(bars, padding=2, fontsize=9)
    ax.set_ylabel("Variant rows")
    ax.set_title("Variant rows by data source", fontweight="bold")
    ax.set_ylim(0, counts.max() * 1.15)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_funcimpact(df, path):
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.1))
    # left: NMD
    nmd = df["nmd_flag"].replace("", "n/a").replace("n.a.", "n/a")
    order = ["no", "NMD-candidate", "n/a"]
    vals = [int((nmd == k).sum()) for k in order]
    cols = [C_TEAL, C_AR, C_GREY]
    b = axes[0].bar(["retained\nORF", "NMD-\ncandidate", "n/a"], vals, color=cols)
    axes[0].bar_label(b, padding=2, fontsize=9)
    axes[0].set_title("Predicted NMD outcome", fontweight="bold")
    axes[0].set_ylabel("variant rows")
    # right: change class
    cc = df["change_class"].replace("", "(none)").value_counts().head(6)[::-1]
    axes[1].barh(range(len(cc)), cc.values, color=C_BLUE)
    axes[1].set_yticks(range(len(cc)))
    axes[1].set_yticklabels([c[:26] for c in cc.index], fontsize=8)
    for i, v in enumerate(cc.values):
        axes[1].text(v, i, f" {v}", va="center", fontsize=8)
    axes[1].set_title("Protein-change class", fontweight="bold")
    axes[1].set_xlim(0, cc.max() * 1.18)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_corroboration(df, path):
    multi = df[df["corroborating_sources"].astype(int) >= 2]
    rows = []
    for g in sorted(multi["gene"].unique()):
        sub = df[df["gene"] == g]
        rows.append((g, int(sub["corroborating_sources"].iloc[0]),
                     sorted(sub["source"].unique())))
    rows.sort(key=lambda r: (-r[1], r[0]))
    genes = [r[0] for r in rows]
    nsrc = [r[1] for r in rows]
    colors = [C_AR if n == 3 else C_BLUE for n in nsrc]
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    y = range(len(genes))[::-1]
    ax.barh(list(y), nsrc, color=colors)
    ax.set_yticks(list(y))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xlabel("number of independent sources")
    ax.set_xticks([0, 1, 2, 3])
    ax.set_title("Cross-source corroborated genes (≥2 sources)",
                 fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C_AR, label="3 sources"),
                       Patch(color=C_BLUE, label="2 sources")],
              loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_arv7(path):
    """Schematic of full-length AR vs AR-V7 domain architecture."""
    fig, ax = plt.subplots(figsize=(7.6, 2.7))
    ax.set_xlim(0, 920)
    ax.set_ylim(0, 3)
    ax.axis("off")

    def domain(x0, x1, y, color, label, txtcol="white"):
        ax.add_patch(FancyBboxPatch((x0, y), x1 - x0, 0.6,
                     boxstyle="round,pad=0.02,rounding_size=6",
                     fc=color, ec="none"))
        ax.text((x0 + x1) / 2, y + 0.3, label, ha="center", va="center",
                fontsize=8, color=txtcol, fontweight="bold")

    # Full-length AR (row top), domains approx: NTD 1-555, DBD 556-623,
    # hinge 624-665, LBD 666-920
    ax.text(-8, 2.0 + 0.3, "AR\n(full-length)", ha="right", va="center",
            fontsize=8, fontweight="bold")
    domain(1, 555, 2.0, C_GREY, "NTD (AF-1 transactivation)")
    domain(555, 623, 2.0, C_BLUE, "DBD")
    domain(623, 665, 2.0, C_GOLD, "H", txtcol="black")
    domain(665, 920, 2.0, C_TEAL, "LBD (ligand binding)")

    # AR-V7 (row bottom): NTD + DBD + short cryptic exon 3, no LBD
    ax.text(-8, 0.7 + 0.3, "AR-V7", ha="right", va="center",
            fontsize=8, fontweight="bold")
    domain(1, 555, 0.7, C_GREY, "NTD (AF-1 transactivation)")
    domain(555, 623, 0.7, C_BLUE, "DBD")
    domain(623, 644, 0.7, C_AR, "CE", txtcol="white")
    # marker for lost region
    ax.plot([665, 920], [1.0, 1.0], color=C_AR, lw=1.2, ls="--")
    ax.annotate("LBD lost → constitutively active,\nenzalutamide/"
                "abiraterone-resistant",
                xy=(790, 1.0), xytext=(560, 0.05), fontsize=8, color=C_AR,
                ha="left")
    ax.set_title("AR-V7: ligand-binding domain loss by cryptic-exon splicing",
                 fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_novelty(path):
    """Cohort prognostic signal for the novel candidates (log2 HR)."""
    # (gene, HR, endpoint, direction); MAP3K7 has no cohort HR (NMD-driven).
    data = [
        ("EXOC7", 5.456, "DFI", "risk-up"),
        ("ENAH", 3.135, "DFI", "risk-up"),
        ("PRUNE2", 0.544, "PFI", "protective"),
        ("STEAP3", 0.458, "PFI", "protective"),
    ]
    fig, ax = plt.subplots(figsize=(7.6, 2.7))
    genes = [d[0] for d in data]
    vals = [math.log2(d[1]) for d in data]
    cols = [C_AR if d[3] == "risk-up" else C_TEAL for d in data]
    y = list(range(len(genes)))[::-1]
    ax.barh(y, vals, color=cols)
    for yi, (g, hr, ep, dr), v in zip(y, data, vals):
        ax.text(v + (0.05 if v >= 0 else -0.05), yi,
                f"{ep} HR {hr:g}", va="center",
                ha="left" if v >= 0 else "right", fontsize=8)
    # MAP3K7 row at top (no HR)
    ax.text(0, len(genes) - 0.0 + 0.0, "", fontsize=8)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=9)
    ax.set_xlabel("log₂ hazard ratio  (TCGA-PRAD, univariate)")
    ax.set_xlim(-1.6, 2.9)
    ax.set_title("Novel candidates: cohort prognostic signal "
                 "(MAP3K7 omitted — NMD-driven, no cohort HR)",
                 fontweight="bold", fontsize=9.5)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C_AR, label="worse outcome"),
                       Patch(color=C_TEAL, label="better outcome")],
              loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Before/after sequences for the featured cases
# ----------------------------------------------------------------------------
def collect_sequences(uni_path=UNI_IMPACT, tcga_path=TCGA_IMPACT) -> list[dict]:
    """Assemble the representative before/after sequence record per featured gene."""
    uni = {r["gene"]: r for r in json.load(open(uni_path))}
    tcga = {o["splice_event"]: o for o in json.load(open(tcga_path))}
    out = []
    for section, gene, source, vid, note in SEQ_REPS:
        rec = {"section": section, "gene": gene, "source": source,
               "src_id": vid, "note": note, "before": "", "after": "",
               "prefix": 0, "suffix": 0, "change_class": "", "available": False}
        if source == "UniProt" and gene in uni:
            r = uni[gene]
            iso = next((i for i in r["isoforms"] if i["isoform_id"] == vid), None)
            if iso:
                rec.update(before=r["canonical_seq"], after=iso["after_seq"],
                           prefix=iso["identical_prefix"],
                           suffix=iso["identical_suffix"],
                           change_class=iso["change_class"], available=True)
        elif source == "TCGA" and vid in tcga:
            o = tcga[vid]
            rec.update(before=o["before_seq"], after=o["after_seq"],
                       prefix=o["identical_prefix"], suffix=o["identical_suffix"],
                       change_class=o["impact_class"], available=True)
        out.append(rec)
    return out


def render_seq_html(seq, prefix, suffix, full_threshold=700, flank=60,
                    change_cap=400) -> str:
    """Monospace HTML with the changed region <mark>-highlighted.

    Long identical flanks (and pathologically long changed regions) are elided
    with markers; complete sequences live in the FASTA export.
    """
    n = len(seq)
    lo, hi = prefix, n - suffix
    lo = max(0, min(lo, n))
    hi = max(lo, min(hi, n))
    esc = _html.escape

    def changed_block(s):
        if len(s) > change_cap:
            k = len(s) - 360
            return (f"<mark>{esc(s[:180])}</mark>"
                    f"<span class='elide'>…[{k} aa, changed]…</span>"
                    f"<mark>{esc(s[-180:])}</mark>")
        return f"<mark>{esc(s)}</mark>"

    if n <= full_threshold:
        return esc(seq[:lo]) + changed_block(seq[lo:hi]) + esc(seq[hi:])
    ps, pe = max(0, lo - flank), min(n, hi + flank)
    parts = []
    if ps > 0:
        parts.append(f"<span class='elide'>…[{ps} aa identical]…</span>")
    parts.append(esc(seq[ps:lo]))
    parts.append(changed_block(seq[lo:hi]))
    parts.append(esc(seq[hi:pe]))
    if pe < n:
        parts.append(f"<span class='elide'>…[{n - pe} aa identical]…</span>")
    return "".join(parts)


def write_fasta(records, path=FASTA) -> int:
    lines, n = [], 0
    for r in records:
        if not r["available"]:
            continue
        tag = f"{r['gene']}|{r['source']}|{r['src_id']}"
        for kind, seq in (("before", r["before"]), ("after", r["after"])):
            lines.append(f">{tag}|{kind} len={len(seq)}")
            lines.extend(seq[i:i + 60] for i in range(0, len(seq), 60))
            n += 1
    open(path, "w").write("\n".join(lines) + "\n")
    return n


# ----------------------------------------------------------------------------
# HTML assembly
# ----------------------------------------------------------------------------
def _abs(p):
    return "file://" + os.path.abspath(p)


def build_html(df, figs, seqs=None) -> str:
    if seqs is None:
        seqs = collect_sequences()
    counts = df["source"].value_counts().to_dict()
    n_rows, n_genes = len(df), df["gene"].nunique()
    n_nmd = int((df["nmd_flag"] == "NMD-candidate").sum())
    n_dom = int((df["domains_hit"] != "").sum())
    multi = df[df["corroborating_sources"].astype(int) >= 2]
    n_multi = multi["gene"].nunique()
    n_3 = df[df["corroborating_sources"].astype(int) == 3]["gene"].nunique()

    # corroboration table
    corr_rows = ""
    seen = set()
    order = sorted(multi["gene"].unique(),
                   key=lambda g: (-int(df[df["gene"] == g]["corroborating_sources"].iloc[0]), g))
    for g in order:
        sub = df[df["gene"] == g]
        srcs = ", ".join(sorted(sub["source"].unique()))
        ev = next((x for x in sub["pca_evidence"] if x), "")
        n = int(sub["corroborating_sources"].iloc[0])
        corr_rows += (f"<tr><td class='g'>{g}</td><td class='c'>{n}</td>"
                      f"<td>{srcs}</td><td>{len(sub)}</td>"
                      f"<td class='ev'>{ev}</td></tr>")

    # top-5 cards
    cards = ""
    for c in TOP5:
        cards += f"""
        <div class="case">
          <div class="case-head">
            <span class="rank">#{c['rank']}</span>
            <span class="cgene">{c['gene']}</span>
            <span class="cvar">{c['variant']}</span>
            <span class="cacc">{c['uniprot']}</span>
          </div>
          <div class="srcline">Evidence: {c['sources']}</div>
          <table class="ctab">
            <tr><th>Splicing event</th><td>{c['event']}</td></tr>
            <tr><th>Predicted protein impact</th><td>{c['impact']}</td></tr>
            <tr><th>Prostate-cancer relevance</th><td>{c['pca']}</td></tr>
            <tr><th>Why it ranks here</th><td>{c['why']}</td></tr>
          </table>
        </div>"""

    # novel-but-credible cards
    ncards = ""
    for c in NOVEL5:
        ncards += f"""
        <div class="case ncase">
          <div class="case-head">
            <span class="rank nrank">N{c['rank']}</span>
            <span class="cgene">{c['gene']}</span>
            <span class="cvar">{c['variant']}</span>
            <span class="cacc">{c['uniprot']}</span>
          </div>
          <div class="srcline nsrc">Evidence: {c['sources']} &nbsp;·&nbsp; {c['signal']}</div>
          <table class="ctab">
            <tr><th>Splicing event</th><td>{c['event']}</td></tr>
            <tr><th>Established PCa context</th><td>{c['known']}</td></tr>
            <tr><th>Why it is novel yet credible</th><td>{c['novel']}</td></tr>
          </table>
        </div>"""

    # before/after sequence appendix
    seqblocks = ""
    for s in seqs:
        if not s["available"]:
            seqblocks += f"""
            <div class="seqcase">
              <div class="seqhead"><b>{s['gene']}</b>
                <span class="seqmeta">{s['section']} · {s['source']} · {s['note']}</span></div>
              <div class="seqna">No before/after protein sequence in the
              pipeline's UniProt or TCGA sources for this literature-curated
              variant; see its PMID provenance in the master table.</div>
            </div>"""
            continue
        bl = len(s["before"])
        al = len(s["after"])
        lo = s["prefix"] + 1
        hi_b = bl - s["suffix"]
        before_html = render_seq_html(s["before"], s["prefix"], s["suffix"])
        after_html = render_seq_html(s["after"], s["prefix"], s["suffix"])
        seqblocks += f"""
        <div class="seqcase">
          <div class="seqhead"><b>{s['gene']}</b>
            <span class="seqmeta">{s['section']} · {s['source']} {s['src_id']}
            · {s['note']} · {s['change_class']} · before {bl} aa → after {al} aa
            · changed from residue {lo}</span></div>
          <div class="seqlabel">Before (canonical/reference)</div>
          <div class="seqblock">{before_html}</div>
          <div class="seqlabel">After (alternative isoform)</div>
          <div class="seqblock">{after_html}</div>
        </div>"""

    today = f"{date.today():%Y-%m-%d}"
    src_items = "".join(
        f"<li><b>{s}</b>: {counts.get(s,0)} variant rows</li>"
        for s in ["UniProt", "TCGA", "ASCancerAtlas", "Literature"])

    return f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<style>
  @page {{ size: A4; margin: 1.7cm 1.6cm 1.9cm 1.6cm;
    @bottom-center {{ content: "Prostate-Cancer Alternative-Splicing Analysis  ·  {today}  ·  page " counter(page) " / " counter(pages);
      font-size: 8pt; color: #888; }} }}
  body {{ font-family: 'DejaVu Sans', Arial, sans-serif; color: #1d1d1f;
    font-size: 10pt; line-height: 1.45; }}
  h1 {{ font-size: 21pt; color: {C_AR}; margin: 0 0 2px 0; }}
  h2 {{ font-size: 14pt; color: {C_BLUE}; border-bottom: 2px solid {C_BLUE};
    padding-bottom: 3px; margin-top: 22px; }}
  h3 {{ font-size: 11.5pt; color: #333; margin: 14px 0 4px; }}
  .sub {{ color: #555; font-size: 10pt; margin-bottom: 2px; }}
  .meta {{ color: #888; font-size: 8.5pt; }}
  .kpis {{ display: flex; gap: 8px; margin: 14px 0 4px; }}
  .kpi {{ flex: 1; background: #f5f6f8; border-radius: 8px; padding: 9px 6px;
    text-align: center; border-top: 3px solid {C_BLUE}; }}
  .kpi .n {{ font-size: 16pt; font-weight: bold; color: {C_AR}; }}
  .kpi .l {{ font-size: 7.8pt; color: #555; text-transform: uppercase;
    letter-spacing: .3px; }}
  ul {{ margin: 4px 0 4px 0; padding-left: 18px; }}
  li {{ margin: 1px 0; }}
  img {{ width: 100%; margin: 8px 0; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 8.6pt;
    margin: 6px 0; }}
  th, td {{ border: 1px solid #d9dce1; padding: 3px 6px; text-align: left;
    vertical-align: top; }}
  th {{ background: {C_BLUE}; color: #fff; font-weight: 600; }}
  td.g {{ font-weight: bold; color: {C_AR}; }}
  td.c {{ text-align: center; font-weight: bold; }}
  td.ev {{ font-size: 8pt; color: #444; }}
  tr:nth-child(even) td {{ background: #f7f8fa; }}
  .case {{ border: 1px solid #e0e0e0; border-left: 5px solid {C_AR};
    border-radius: 6px; padding: 8px 12px; margin: 11px 0;
    background: #fcfcfd; break-inside: avoid; }}
  .case-head {{ display: flex; align-items: baseline; gap: 9px; }}
  .rank {{ background: {C_AR}; color: #fff; font-weight: bold;
    border-radius: 50%; padding: 2px 9px; font-size: 11pt; }}
  .cgene {{ font-size: 13pt; font-weight: bold; color: #111; }}
  .cvar {{ font-size: 10.5pt; color: {C_BLUE}; font-weight: 600; }}
  .cacc {{ font-size: 8.5pt; color: #999; margin-left: auto; }}
  .srcline {{ font-size: 8.4pt; color: {C_TEAL}; font-weight: 600;
    margin: 2px 0 5px; }}
  .ncase {{ border-left: 5px solid {C_BLUE}; }}
  .nrank {{ background: {C_BLUE}; }}
  .nsrc {{ color: {C_BLUE}; }}
  .ctab th {{ background: #eef1f5; color: #333; width: 27%;
    font-weight: 600; }}
  .ctab td {{ background: #fff; }}
  .note {{ background: #fff8ec; border: 1px solid {C_GOLD}; border-radius: 6px;
    padding: 8px 11px; font-size: 8.8pt; margin: 10px 0; }}
  .limit li {{ font-size: 9pt; }}
  .pagebreak {{ break-before: page; }}
  .seqcase {{ margin: 9px 0; break-inside: avoid; }}
  .seqhead {{ font-size: 10pt; color: #111; border-bottom: 1px solid #ddd;
    padding-bottom: 2px; }}
  .seqhead b {{ color: {C_AR}; font-size: 11pt; }}
  .seqmeta {{ font-size: 7.6pt; color: #666; font-weight: normal; }}
  .seqlabel {{ font-size: 7.8pt; color: {C_BLUE}; font-weight: bold;
    margin: 4px 0 1px; text-transform: uppercase; letter-spacing: .3px; }}
  .seqblock {{ font-family: 'DejaVu Sans Mono', monospace; font-size: 6.4pt;
    line-height: 1.35; white-space: pre-wrap; word-break: break-all;
    background: #fafbfc; border: 1px solid #e6e8eb; border-radius: 4px;
    padding: 5px 7px; }}
  .seqblock mark {{ background: #ffe08a; color: #6b3a00; }}
  .elide {{ color: #b03a48; font-style: italic; background: #f3e6e8;
    padding: 0 3px; border-radius: 3px; }}
  .seqna {{ font-size: 8.4pt; color: #777; font-style: italic;
    padding: 4px 0; }}
</style></head><body>

<h1>Alternative Splicing in Prostate Cancer</h1>
<div class="sub">A multi-source integration of curated isoforms, cohort splicing
events, literature variants and normal-tissue baseline &mdash; with predicted
protein-functional impact</div>
<div class="meta">Generated {today} &nbsp;·&nbsp; Sources: UniProt &middot;
TCGA SpliceSeq PRAD &middot; ASCancerAtlas &middot; literature (EuropePMC) &middot;
GTEx v8 baseline</div>

<div class="kpis">
  <div class="kpi"><div class="n">{n_rows}</div><div class="l">Variant rows</div></div>
  <div class="kpi"><div class="n">{n_genes}</div><div class="l">Genes</div></div>
  <div class="kpi"><div class="n">{n_dom}</div><div class="l">Domain-disrupting</div></div>
  <div class="kpi"><div class="n">{n_nmd}</div><div class="l">NMD-candidates</div></div>
  <div class="kpi"><div class="n">{n_multi}</div><div class="l">Multi-source genes</div></div>
</div>

<h2>1 &nbsp; Executive summary</h2>
<p>This report integrates <b>{n_rows} alternative-splicing variant records across
{n_genes} genes</b> from four independent sources into one schema, then annotates
each with predicted protein-functional consequences (domain disruption,
nonsense-mediated-decay candidacy, localization signals) using UniProt feature
tracks, and contextualises them against the GTEx normal-prostate isoform
baseline. {n_dom} records disrupt at least one annotated domain or functional
region and {n_nmd} are predicted nonsense-mediated-decay candidates.
<b>{n_multi} genes are corroborated by &ge;2 independent sources</b> (of which
{n_3} by all-but-one), forming the highest-confidence shortlist. The androgen
receptor (AR) splicing axis &mdash; headlined by AR-V7 &mdash; is the dominant
splicing-driven mechanism of therapeutic resistance and anchors the top of our
case ranking.</p>

<ul>
{src_items}
</ul>
<img src="{_abs(figs['sources'])}"/>

<h2>2 &nbsp; Methods (brief)</h2>
<p>Every source is normalised into a 20-column master schema
(<code>master_variants.tsv</code>). For each variant the canonical and alternative
protein sequences are aligned to a common prefix/suffix; the resulting
<i>changed interval</i> is intersected with UniProt feature tracks (domains,
binding sites, zinc-fingers, DNA-binding regions, signal/transmembrane segments)
to call domain disruption and localization changes. A C-terminal truncation
losing &gt;50&nbsp;aa flags a protein-level NMD-candidate (heuristic). Genes are
cross-tabulated by source to compute corroboration. GTEx v8 median transcript
expression provides a normal-prostate isoform baseline.</p>
<div class="note"><b>Interpretation note.</b> Disease, resistance, NMD and
domain-disruption statements are <b>predictions or literature-curated
annotations</b>, not measurements produced by this pipeline. Cohort hazard ratios
are univariate, single-cohort, and some Cox fits are numerically unstable; they
are used only as supporting context. CancerSplicingQTL was unreachable
(HTTP&nbsp;403) and is omitted.</div>

<h2>3 &nbsp; Predicted functional-impact landscape</h2>
<img src="{_abs(figs['func'])}"/>
<p>In-frame substituted segments and internal deletions dominate, but a
substantial tail of frameshift and truncating events feeds the
{n_nmd}-strong NMD-candidate pool &mdash; isoforms most likely to remove protein
rather than remodel it.</p>

<h2>4 &nbsp; Cross-source corroboration</h2>
<p>Genes recovered independently by more than one source are the most defensible
leads. Three genes &mdash; <b>AR, ITGA6, CTNND1</b> &mdash; appear in three
sources each.</p>
<img src="{_abs(figs['corr'])}"/>
<table>
  <tr><th>Gene</th><th>Sources</th><th>Source set</th><th>Variants</th><th>Best PCa evidence</th></tr>
  {corr_rows}
</table>

<div class="pagebreak"></div>
<h2>5 &nbsp; Top 5 splicing cases most relevant to prostate cancer</h2>
<p>Ranked by strength and independence of evidence, fidelity of the predicted
mechanism, and established prostate-cancer / therapy-resistance relevance.</p>
<img src="{_abs(figs['arv7'])}"/>
{cards}

<div class="pagebreak"></div>
<h2>6 &nbsp; Novel but credible candidates</h2>
<p>The cases in Section 5 are, by design, well-established biology. This section
asks the complementary question: which splicing events are <b>under-studied in
prostate cancer yet credible</b>? Each gene below is independently relevant to
prostate cancer (a known suppressor, antigen family, or invasion driver) while
its <i>splicing</i> consequence is largely unexplored in the disease — surfaced
here by cross-source corroboration, a cohort signal, or a predicted
NMD/domain hit. These are hypothesis-generating leads, not established
mechanisms.</p>
<img src="{_abs(figs['novel'])}"/>
{ncards}
<div class="note"><b>Credibility caveat.</b> The supporting hazard ratios are
univariate, single-cohort TCGA-PRAD signals and the NMD/domain calls are
predictions; each candidate needs exon-level validation. EXOC7 carries the
strongest signal but should be treated as the most speculative-by-magnitude.</div>

<h2>7 &nbsp; Limitations</h2>
<ul class="limit">
  <li>Functional calls (NMD, domain disruption, disease/resistance links) are
  predictions or curated annotations, not measurements from this analysis.</li>
  <li>Cross-source corroboration is gene-level, not exon/variant-level; the same
  gene may be hit by different events in different sources.</li>
  <li>UniProt isoforms and literature variants are curated presence/absence;
  GTEx provides normal context only; TCGA cohort HRs are univariate,
  single-cohort, and occasionally numerically degenerate.</li>
  <li>The NMD rule is a protein-length heuristic, not a transcript-level
  55-nt-rule evaluation.</li>
  <li>CancerSplicingQTL was unreachable (HTTP 403) and is omitted.</li>
</ul>

<h2>8 &nbsp; Reproducibility</h2>
<p>Generated from <code>results_collected/master_variants.tsv</code> and
<code>gtex_prostate_baseline.tsv</code> via
<code>analysis/build_pdf_report.py</code>. Collection and integration steps are
documented in <code>results_collected/README.md</code>. Caches replay offline
with <code>--offline</code>.</p>

<div class="pagebreak"></div>
<h2>Appendix A &nbsp; Before / after protein sequences</h2>
<p>For each featured case (Sections 5 &amp; 6), the canonical/reference protein
sequence and the alternative isoform are shown with the
<mark>changed region highlighted</mark>. Long identical flanks &mdash; and any
unusually long changed stretch &mdash; are elided with
<span class="elide">…[N aa]…</span> markers for readability; the
<b>complete</b> before/after sequences are exported to
<code>results_collected/featured_sequences.fasta</code>. UniProt cases compare
the canonical against the curated isoform; TCGA cases compare the protein
inferred for the splice-event endpoints.</p>
{seqblocks}

</body></html>"""


def build(master_tsv=MASTER, out_pdf=OUT_PDF) -> dict:
    os.makedirs(FIGDIR, exist_ok=True)
    df = pd.read_csv(master_tsv, sep="\t").fillna("")
    df["corroborating_sources"] = df["corroborating_sources"].replace("", 0).astype(int)
    figs = {
        "sources": os.path.join(FIGDIR, "f1_sources.png"),
        "func": os.path.join(FIGDIR, "f2_funcimpact.png"),
        "corr": os.path.join(FIGDIR, "f3_corroboration.png"),
        "arv7": os.path.join(FIGDIR, "f4_arv7.png"),
        "novel": os.path.join(FIGDIR, "f5_novelty.png"),
    }
    fig_sources(df, figs["sources"])
    fig_funcimpact(df, figs["func"])
    fig_corroboration(df, figs["corr"])
    fig_arv7(figs["arv7"])
    fig_novelty(figs["novel"])
    seqs = collect_sequences()
    n_fasta = write_fasta(seqs)
    html = build_html(df, figs, seqs)
    HTML(string=html, base_url=".").write_pdf(out_pdf)
    return {"rows": len(df), "genes": int(df["gene"].nunique()),
            "pdf": out_pdf, "top5": [c["gene"] for c in TOP5],
            "fasta_seqs": n_fasta}


def main() -> None:
    print(build())


if __name__ == "__main__":
    main()

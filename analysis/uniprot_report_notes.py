"""Curated PCa / therapy-resistance notes for the featured UniProt proteins, plus
the TCGA progression cross-reference for the seven overlap genes.

UniProt already supplies each protein's function and disease text; these notes
add the prostate-cancer / treatment-resistance reasoning (Part 3) that UniProt
does not editorialize, and connect the curated isoforms to the cohort signal.
"""

# gene -> therapy/PCa-resistance interpretation (hypothesis-level, literature-based)
RESISTANCE_NOTES = {
    "AR": (
        "The androgen receptor is the central driver of prostate cancer and the "
        "target of every AR-pathway inhibitor (enzalutamide, apalutamide, "
        "abiraterone). Isoform 3 here is AR-V7 (644 aa): the alternative splice "
        "removes the C-terminal ligand-binding domain and replaces it with a "
        "short cryptic-exon tail, yielding a constitutively active receptor that "
        "no longer needs androgen and cannot be blocked by ligand-competitive "
        "drugs. AR-V7 is the best-established mechanism of resistance to "
        "enzalutamide/abiraterone and a clinical biomarker. Crucially, this "
        "splice variant is INVISIBLE to the TCGA arm of this project (GRCh37 r75 "
        "lacks the AR-V cryptic exons) — UniProt captures exactly the event the "
        "cohort analysis flagged as its top limitation."
    ),
    "PTEN": (
        "PTEN is the principal lipid-phosphatase brake on PI3K/AKT signalling and "
        "one of the most frequently lost tumour suppressors in prostate cancer; "
        "PTEN loss drives castration resistance and cross-talk that bypasses AR "
        "blockade. A splice isoform that removes or alters the phosphatase domain "
        "would phenocopy PTEN loss and predict PI3K/AKT-pathway dependence "
        "(rationale for AKT-inhibitor combinations)."
    ),
    "CHEK2": (
        "CHEK2 is a DNA-damage-response kinase; germline CHEK2 truncations raise "
        "prostate-cancer risk and, like other DDR defects, may sensitise tumours "
        "to PARP inhibitors and platinum. Splice isoforms that disrupt the kinase "
        "domain are predicted loss-of-function and therefore DDR-deficiency "
        "candidates."
    ),
    "ELAC2": (
        "ELAC2/HPC2 was the first hereditary-prostate-cancer susceptibility gene "
        "mapped; it is a tRNA-processing endonuclease. Altered isoforms are "
        "candidates for the inherited-risk component of disease."
    ),
    "RNASEL": (
        "RNASEL/HPC1 mediates the interferon antiviral/apoptotic response; "
        "germline defects underlie hereditary prostate cancer. Splice changes "
        "that impair the ribonuclease are plausible risk alleles."
    ),
    "MSR1": (
        "MSR1 (macrophage scavenger receptor) variants segregate with hereditary "
        "prostate cancer, linking innate-immune/macrophage biology to disease "
        "predisposition; isoform changes may modulate receptor function."
    ),
    "MSMB": (
        "MSMB encodes PSP94/beta-microseminoprotein, an abundant prostate "
        "secretory protein and the source of the PSA-independent prostate-cancer "
        "biomarker MSMB/-57C>T locus; isoform usage affects secreted-protein "
        "output."
    ),
    "EHBP1": (
        "EHBP1 sits at a hereditary-prostate-cancer (HPC12) susceptibility locus "
        "and links endocytic trafficking to disease risk."
    ),
    "HNF1B": (
        "HNF1B is a transcription factor at a replicated prostate-cancer risk "
        "locus (HPC11); isoform balance can shift its transcriptional output."
    ),
    "MXI1": (
        "MXI1 antagonizes MYC; loss of MXI1 restraint promotes proliferation. "
        "Isoforms that weaken MYC antagonism would favour tumour growth."
    ),
    "EPHB2": (
        "EPHB2 is a receptor tyrosine kinase with tumour-suppressor behaviour in "
        "prostate epithelium; truncating isoforms are candidate loss-of-function "
        "events."
    ),
    "KLF6": (
        "KLF6 is a tumour-suppressor transcription factor; a well-known "
        "alternative-splice isoform (KLF6-SV1) is oncogenic and associated with "
        "aggressive, treatment-refractory prostate cancer — a direct example of "
        "splicing flipping a suppressor into a driver."
    ),
    # cross-reference overlap genes (also progression-associated in the cohort)
    "PRUNE2": (
        "PRUNE2 is an established prostate tumour suppressor regulated by the "
        "PCA3 lncRNA. In the TCGA cohort its protein-changing exon-skip is "
        "risk-DOWN (protective-isoform direction), consistent with suppressor "
        "biology."
    ),
    "ITGA6": (
        "Integrin alpha-6 (laminin receptor) drives adhesion, invasion and "
        "metastasis; its protein-changing exon-skip is risk-UP in the cohort, "
        "matching a pro-invasive role."
    ),
    "CTNND1": (
        "p120-catenin stabilizes E-cadherin junctions; its cohort exon-skip is "
        "risk-DOWN, consistent with adhesion maintenance."
    ),
    "FES": (
        "FES tyrosine kinase is pro-proliferative; its alt-donor event is risk-UP "
        "in the cohort."
    ),
    "STEAP3": (
        "STEAP3 metalloreductase belongs to the STEAP family of prostate "
        "antigens/therapeutic targets; its cohort event is risk-DOWN."
    ),
    "SMN1": (
        "SMN1 assembles spliceosomal snRNPs — a splicing-machinery feedback node; "
        "its cohort alt-acceptor event is risk-UP."
    ),
    "PRICKLE4": (
        "PRICKLE4 is a planar-cell-polarity/directed-migration component and the "
        "strongest cohort lead (lowest BH p, HR 4.35, risk-UP)."
    ),
}

# TCGA cross-reference (from analysis/crossref_uniprot_events.tsv): gene -> dict
TCGA_XREF = {
    "PRICKLE4": dict(event="PRICKLE4_AA_76141", type="AA", endpoint="PFI",
                     HR=4.348, direction="risk-up", bh_p=5.7e-04),
    "ITGA6": dict(event="ITGA6_ES_55968", type="ES", endpoint="PFI",
                  HR=2.474, direction="risk-up", bh_p=2.3e-03),
    "CTNND1": dict(event="CTNND1_ES_15935", type="ES", endpoint="PFI",
                   HR=0.45, direction="risk-down", bh_p=3.4e-03),
    "FES": dict(event="FES_AD_32502", type="AD", endpoint="PFI",
                HR=3.589, direction="risk-up", bh_p=2.4e-02),
    "STEAP3": dict(event="STEAP3_ME_95656", type="ME", endpoint="PFI",
                   HR=0.458, direction="risk-down", bh_p=3.9e-02),
    "SMN1": dict(event="SMN1_AA_72422", type="AA", endpoint="PFI",
                 HR=1.691, direction="risk-up", bh_p=4.5e-02),
    "PRUNE2": dict(event="PRUNE2_ES_86643", type="ES", endpoint="PFI",
                   HR=0.544, direction="risk-down", bh_p=4.6e-02),
}

# Display order for the featured section: hereditary/curated PCa-disease genes
# first (AR leads — it carries AR-V7), then the cohort-overlap genes.
FEATURED_PCA_DISEASE = [
    "AR", "PTEN", "KLF6", "CHEK2", "ELAC2", "RNASEL", "MSR1", "MSMB",
    "EHBP1", "HNF1B", "MXI1", "EPHB2",
]
FEATURED_XREF = ["PRICKLE4", "ITGA6", "CTNND1", "FES", "STEAP3", "SMN1", "PRUNE2"]


def resistance_note(gene: str) -> str:
    return RESISTANCE_NOTES.get(gene, "")

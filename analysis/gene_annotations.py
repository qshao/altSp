"""Curated, literature-level functional annotations for the featured genes.

These are established protein functions and prostate-cancer (PCa) context used to
reason about splicing consequences. They are background knowledge, not derived
from this dataset; the report frames all downstream links as hypotheses.
"""

# gene -> (known_function, prostate_cancer_relevance)
ANNOTATIONS: dict[str, tuple[str, str]] = {
    "ACLY": (
        "ATP-citrate lyase: the enzyme that converts cytosolic citrate to "
        "acetyl-CoA, the rate-limiting source of carbon for de novo fatty-acid "
        "and cholesterol synthesis. Central node of lipogenic metabolism.",
        "De novo lipogenesis is a hallmark of prostate cancer and is driven by "
        "androgen-receptor (AR) signalling; ACLY is frequently up-regulated and "
        "supports membrane and steroid-precursor supply. A protein-altering "
        "splice that changes ACLY activity could rewire the lipid metabolism on "
        "which AR-driven and castration-resistant tumours depend.",
    ),
    "MAPT": (
        "Microtubule-associated protein tau: binds and stabilises "
        "microtubules; its own pre-mRNA is heavily alternatively spliced.",
        "Tau expression modulates the binding of taxanes (docetaxel, "
        "cabazitaxel) to microtubules; high/altered tau is a reported marker of "
        "taxane resistance. A tau isoform switch is therefore directly relevant "
        "to chemotherapy resistance in advanced prostate cancer.",
    ),
    "SPOP": (
        "Speckle-type POZ protein: substrate-recognition subunit of a Cullin-3 "
        "E3 ubiquitin ligase. Targets AR, the AR coactivator SRC-3/NCOA3, and "
        "BET proteins (BRD4) for degradation.",
        "SPOP is the most frequently mutated gene in primary prostate cancer. "
        "Loss or alteration of functional SPOP stabilises AR and its "
        "coactivators, amplifying AR output — a core driver of progression and "
        "of resistance to AR-pathway inhibitors. An isoform that disrupts "
        "substrate binding would phenocopy inactivating mutation.",
    ),
    "KLK2": (
        "Kallikrein-2: a secreted, AR-regulated serine protease of the same "
        "family as PSA/KLK3; activates pro-PSA and components of the tumour "
        "micro-environment proteolytic cascade.",
        "KLK2 is a canonical AR target and prostate biomarker (basis of some "
        "PSA-family and PSMA-adjacent assays). An altered protease product "
        "reflects rewired AR-driven secretory programmes and could affect "
        "tumour-microenvironment remodelling.",
    ),
    "SMARCD3": (
        "BAF60c: a subunit of the SWI/SNF (BAF) ATP-dependent chromatin-"
        "remodelling complex that controls enhancer accessibility and "
        "lineage-specific transcription.",
        "SWI/SNF subunits are recurrently altered across cancers and shape AR "
        "and lineage-plasticity programmes implicated in neuroendocrine "
        "trans-differentiation — a route to treatment-resistant prostate "
        "cancer. A near-complete loss of this subunit would impair normal BAF "
        "function.",
    ),
    "OGG1": (
        "8-oxoguanine DNA glycosylase: initiates base-excision repair of "
        "oxidative DNA damage (8-oxoG), guarding genome integrity.",
        "DNA-repair capacity governs both mutational burden and response to "
        "DNA-damaging therapy. Altered OGG1 could raise oxidative mutagenesis "
        "(fuelling progression) and modulate sensitivity to radiation and "
        "PARP/platinum strategies used in DNA-repair-deficient prostate cancer.",
    ),
    "RBM6": (
        "RNA-binding motif protein 6: an RNA-binding/splicing-associated factor "
        "at the 3p21.3 tumour-suppressor locus.",
        "Splicing regulators that are themselves mis-spliced can propagate "
        "genome-wide splicing changes. Truncation of an RBM-family regulator "
        "may shift the splicing landscape toward oncogenic isoforms.",
    ),
    "RAP1GAP": (
        "Rap1 GTPase-activating protein: switches off Rap1 signalling, "
        "restraining integrin activation, cell adhesion and migration.",
        "RAP1GAP loss is reported in several carcinomas and de-represses Rap1, "
        "promoting invasion and migration — phenotypes tied to metastatic, "
        "treatment-refractory disease.",
    ),
    "VCL": (
        "Vinculin: a core focal-adhesion and adherens-junction protein linking "
        "integrins/cadherins to the actin cytoskeleton; restrains motility.",
        "Loss of adhesion-complex integrity promotes invasion and metastasis. "
        "An altered vinculin could weaken cell–matrix anchorage, favouring "
        "dissemination.",
    ),
    "SVIL": (
        "Supervillin: a membrane-associated actin-binding protein coupling the "
        "cytoskeleton to the plasma membrane and to myosin II contractility.",
        "Cytoskeletal remodelling underlies invasion and the cellular plasticity "
        "associated with progression.",
    ),
    "DOCK7": (
        "Dedicator of cytokinesis 7: a guanine-nucleotide exchange factor "
        "activating Rac1/Cdc42 to drive cytoskeletal dynamics and polarity.",
        "Rac/Cdc42 GEF activity feeds migration and invasion programmes.",
    ),
    "ADHFE1": (
        "Hydroxyacid-oxoacid transhydrogenase: produces the oncometabolite "
        "D-2-hydroxyglutarate (D-2HG), which can inhibit α-ketoglutarate-"
        "dependent dioxygenases and reshape the epigenome.",
        "Oncometabolite accumulation links metabolism to epigenetic "
        "dysregulation; altered ADHFE1 activity could perturb this axis.",
    ),
    "SEC31A": (
        "Outer-coat subunit of the COPII vesicle complex mediating "
        "endoplasmic-reticulum-to-Golgi protein export.",
        "Secretory-pathway capacity supports the heavy secretory load of "
        "prostate epithelium (PSA/KLK family) and surface-receptor trafficking.",
    ),
    "EXOC7": (
        "Exo70: a subunit of the exocyst complex that targets secretory "
        "vesicles to the plasma membrane and promotes invadopodia formation.",
        "Exocyst-driven membrane trafficking supports directed secretion and "
        "invasion.",
    ),
    "QTRT1": (
        "Queuine tRNA-ribosyltransferase catalytic subunit: installs the "
        "queuosine modification in tRNA anticodons, tuning translational "
        "fidelity and speed.",
        "tRNA-modification enzymes influence the translation of growth- and "
        "stress-response proteins; their dysregulation is increasingly linked "
        "to tumour proliferation.",
    ),
    "HAUS5": (
        "Subunit of the augmin/HAUS complex that nucleates spindle "
        "microtubules; required for accurate mitosis.",
        "Mitotic-spindle integrity affects proliferation and the response to "
        "anti-mitotic taxane chemotherapy.",
    ),
    "INO80E": (
        "Subunit of the INO80 ATP-dependent chromatin-remodelling complex "
        "governing nucleosome positioning, transcription and DNA repair.",
        "Chromatin-remodeller dysfunction contributes to transcriptional "
        "reprogramming and lineage plasticity in resistant disease.",
    ),
    "KIF12": (
        "Kinesin-family motor protein involved in microtubule-based transport.",
        "Motor-protein alterations can affect mitosis and intracellular "
        "trafficking.",
    ),
    "KLHDC2": (
        "Kelch-domain-containing 2: a substrate receptor of a Cullin-2 E3 "
        "ubiquitin ligase recognising C-terminal degrons.",
        "Altered E3-ligase substrate selection can stabilise or destabilise "
        "growth-regulatory proteins.",
    ),
    "FADS3": (
        "Fatty-acid desaturase 3: introduces double bonds during "
        "polyunsaturated-fatty-acid and sphingolipid biosynthesis.",
        "Lipid-desaturation flux is part of the lipogenic phenotype of prostate "
        "cancer.",
    ),
    "NRBP2": (
        "Nuclear receptor-binding protein 2: a pseudokinase adaptor implicated "
        "in cell-survival signalling.",
        "Survival-pathway adaptors can modulate apoptosis sensitivity.",
    ),
    "FAM73B": (
        "MIGA2 / mitoguardin-2: an outer-mitochondrial-membrane protein "
        "promoting mitochondrial fusion and lipid-droplet/ER contacts.",
        "Mitochondrial-dynamics and lipid handling support the metabolic "
        "demands of tumour growth.",
    ),
    "ZNF691": (
        "C2H2 zinc-finger transcription factor (function incompletely "
        "characterised).",
        "Transcription-factor isoform changes can rewire gene-expression "
        "programmes.",
    ),
    "PRICKLE4": (
        "Planar-cell-polarity pathway component (Prickle family).",
        "Planar-cell-polarity signalling influences directed migration.",
    ),
    "WASH4P": (
        "WASH-family pseudogene-derived product; WASH complex drives "
        "actin nucleation on endosomes (annotation uncertain for this locus).",
        "Endosomal actin dynamics affect receptor recycling and motility.",
    ),
}


def get(gene: str) -> tuple[str, str]:
    return ANNOTATIONS.get(
        gene,
        ("Function not specifically curated here; consequence inferred from the "
         "sequence change alone.",
         "Relevance to prostate cancer not specifically curated; included on the "
         "strength of its statistical progression association."),
    )

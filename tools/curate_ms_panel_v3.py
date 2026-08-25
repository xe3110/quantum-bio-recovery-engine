"""Curate the v3 multiple-sclerosis candidate panel.

Every record is a *directional mechanistic hypothesis* assembled from review
literature (see ``SOURCES``).  ``target_effects`` are signed, unitless values
in [-1, 1] describing the direction and relative magnitude by which a drug is
hypothesised to move a transcript's activity.  They are deliberately NOT
affinities, exposures, doses, or clinical response estimates.

Run ``python -m tools.curate_ms_panel_v3`` to regenerate
``data/drugs/ms_panel_v3.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SOURCES = [
    "PMID:39376160",  # Dalla Costa G, et al. Expert Rev Clin Pharmacol 2024
    "PMID:39177889",  # Olejnik P, et al. Pharmacol Rep 2024
    "PMID:38719481",  # Bsteh G, et al. Pharmacol Rev 2024
    "PMID:38480630",  # Carlson AK, et al. Drugs 2024
    "PMCID:PMC8218202",  # Zheng S, et al. DrugComb, Nucleic Acids Res 2021
]

# Axes a candidate is hypothesised to act on. Used for axis-coverage scoring.
AXES = ("immunomodulation", "cns_innate", "remyelination", "neuroprotection", "metabolic_repair")

# safety_burden keys, each 0 (no meaningful signal) .. 1 (major, ranking-relevant)
RISKS = ("infection", "malignancy", "cardiac", "hepatic", "teratogenicity", "autoimmunity", "ocular")


def drug(
    name, mechanism_class, moa, tier, route, cns, compartment, axes, half_life_class,
    safety_classes, risks, uncertainty, pathways, targets, indication="RRMS", note="",
):
    """Build one panel record, validating the fields that scoring depends on."""
    assert tier in {"approved", "phase_3", "phase_2", "preclinical"}, (name, tier)
    assert 0.0 <= cns <= 1.0, name
    assert compartment in {"peripheral", "cns", "both"}, (name, compartment)
    assert set(axes) <= set(AXES), (name, axes)
    assert set(risks) <= set(RISKS), (name, risks)
    assert half_life_class in {"short", "intermediate", "long", "reconstitution"}, (name, half_life_class)
    assert all(-1.0 <= v <= 1.0 for v in targets.values()), name
    return {
        "name": name,
        "mechanism_class": mechanism_class,
        "primary_moa": moa,
        "evidence_tier": tier,
        "ms_indication": indication,
        "route": route,
        "cns_penetration": cns,
        "compartment": compartment,
        "therapeutic_axes": list(axes),
        "half_life_class": half_life_class,
        "safety_classes": list(safety_classes),
        "safety_burden": {k: round(float(risks.get(k, 0.0)), 2) for k in RISKS},
        "target_uncertainty": uncertainty,
        "pathways": list(pathways),
        "relevant_pathways": list(pathways),
        "target_effects": {g: round(float(v), 3) for g, v in targets.items()},
        "curation_note": note,
    }


# Agents that engage the same biological target despite carrying different
# mechanism_class labels (modality differs, pharmacodynamics do not). Pairs
# within a family are pharmacodynamically redundant and must not be ranked as
# combinations. Anything not listed falls back to its own mechanism_class.
TARGET_FAMILY = {
    "Natalizumab": "alpha4_integrin",
    "Natalizumab-biosimilar": "alpha4_integrin",
    "Firategrast": "alpha4_integrin",
    "Teriflunomide": "DHODH",
    "Vidofludimus calcium": "DHODH",
    "Cladribine": "purine_antimetabolite",
    "Azathioprine": "purine_antimetabolite",
    "Ocrelizumab": "CD20_CD19_B_depletion",
    "Ofatumumab": "CD20_CD19_B_depletion",
    "Ublituximab": "CD20_CD19_B_depletion",
    "Rituximab": "CD20_CD19_B_depletion",
    "Inebilizumab": "CD20_CD19_B_depletion",
    "Dimethyl fumarate": "Nrf2_activation",
    "Diroximel fumarate": "Nrf2_activation",
    "Clemastine": "muscarinic_OPC",
    "Benztropine": "muscarinic_OPC",
    "Lamotrigine": "sodium_channel",
    "Phenytoin": "sodium_channel",
    "Oxcarbazepine": "sodium_channel",
    "Tofacitinib": "JAK",
    "Baricitinib": "JAK",
    "Methylprednisolone": "glucocorticoid",
    "Clobetasol": "glucocorticoid",
}

DRUGS: list[dict] = []

# ---------------------------------------------------------------------------
# Tier 1 — approved MS disease-modifying therapies
# ---------------------------------------------------------------------------
DRUGS += [
    drug("Alemtuzumab", "CD52_depletion", "Anti-CD52 lymphocyte depletion", "approved",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "reconstitution",
         ["broad_immunosuppression", "infusion_reaction"],
         {"infection": 0.8, "malignancy": 0.5, "autoimmunity": 0.9, "cardiac": 0.3},
         0.25, ["adaptive_immunity", "B_cell_immunity"],
         {"CD52": -0.85, "IL2RA": -0.45, "IFNG": -0.35, "IL17A": -0.35, "CD40": -0.25,
          "MS4A1": -0.40, "CD19": -0.35, "TBX21": -0.30, "RORC": -0.30, "FOXP3": 0.20},
         note="Profound depletion followed by reconstitution; high secondary autoimmunity."),

    drug("Cladribine", "purine_antimetabolite", "Selective lymphocyte depletion", "approved",
         "oral", 0.30, "both", ["immunomodulation"], "reconstitution",
         ["broad_immunosuppression"],
         {"infection": 0.6, "malignancy": 0.4, "teratogenicity": 0.9},
         0.30, ["adaptive_immunity", "B_cell_immunity"],
         {"IL2RA": -0.40, "MS4A1": -0.45, "CD19": -0.40, "CD40": -0.30, "IFNG": -0.25,
          "CXCL13": -0.30, "CD27": -0.35}),

    drug("Dimethyl fumarate", "Nrf2_activator", "Nrf2 pathway activation, immune shift", "approved",
         "oral", 0.55, "both", ["immunomodulation", "metabolic_repair"], "short",
         ["lymphopenia", "gi_intolerance"],
         {"infection": 0.4, "hepatic": 0.3},
         0.25, ["oxidative_stress", "adaptive_immunity", "metabolism"],
         {"NFE2L2": 0.75, "HMOX1": 0.70, "NQO1": 0.65, "GPX1": 0.40, "GSR": 0.40,
          "TNF": -0.35, "IL1B": -0.30, "NOS2": -0.40, "IL17A": -0.30, "IL10": 0.30}),

    drug("Diroximel fumarate", "Nrf2_activator", "Nrf2 activation, improved GI tolerability", "approved",
         "oral", 0.55, "both", ["immunomodulation", "metabolic_repair"], "short",
         ["lymphopenia"],
         {"infection": 0.4, "hepatic": 0.25},
         0.30, ["oxidative_stress", "adaptive_immunity", "metabolism"],
         {"NFE2L2": 0.72, "HMOX1": 0.68, "NQO1": 0.62, "GPX1": 0.38,
          "TNF": -0.32, "IL1B": -0.28, "NOS2": -0.38, "IL17A": -0.28}),

    drug("Fingolimod", "S1P_modulator", "Non-selective S1P receptor modulator", "approved",
         "oral", 0.60, "both", ["immunomodulation", "cns_innate"], "long",
         ["s1p_modulation", "bradyarrhythmia"],
         {"infection": 0.5, "cardiac": 0.6, "ocular": 0.5, "hepatic": 0.3},
         0.25, ["lymphocyte_trafficking", "adaptive_immunity", "remyelination"],
         {"S1PR1": -0.80, "S1PR5": -0.45, "SELL": -0.45, "IFNG": -0.30,
          "IL17A": -0.35, "CCL2": -0.25, "MBP": 0.15, "OLIG2": 0.10}),

    drug("Ozanimod", "S1P_modulator", "S1P1/S1P5-selective modulator", "approved",
         "oral", 0.50, "both", ["immunomodulation"], "long",
         ["s1p_modulation", "bradyarrhythmia"],
         {"infection": 0.4, "cardiac": 0.3, "ocular": 0.35, "hepatic": 0.3},
         0.25, ["lymphocyte_trafficking", "adaptive_immunity"],
         {"S1PR1": -0.78, "S1PR5": -0.50, "SELL": -0.40, "IFNG": -0.28, "IL17A": -0.30}),

    drug("Ponesimod", "S1P_modulator", "S1P1-selective modulator, rapid reversibility", "approved",
         "oral", 0.45, "peripheral", ["immunomodulation"], "short",
         ["s1p_modulation", "bradyarrhythmia"],
         {"infection": 0.4, "cardiac": 0.3, "hepatic": 0.25},
         0.25, ["lymphocyte_trafficking", "adaptive_immunity"],
         {"S1PR1": -0.82, "SELL": -0.42, "IFNG": -0.28, "IL17A": -0.28}),

    drug("Siponimod", "S1P_modulator", "S1P1/S1P5 modulator with CNS activity", "approved",
         "oral", 0.65, "both", ["immunomodulation", "cns_innate", "neuroprotection"], "intermediate",
         ["s1p_modulation", "bradyarrhythmia"],
         {"infection": 0.4, "cardiac": 0.4, "ocular": 0.35, "hepatic": 0.3, "malignancy": 0.2},
         0.25, ["lymphocyte_trafficking", "microglial_activation", "remyelination", "neuroprotection"],
         {"S1PR1": -0.75, "S1PR5": -0.60, "SELL": -0.40, "AIF1": -0.30, "CD68": -0.30,
          "IL17A": -0.30, "MBP": 0.25, "PLP1": 0.20, "OLIG2": 0.20, "NEFL": -0.25},
         indication="SPMS"),

    drug("Glatiramer acetate", "immune_deviation", "Random copolymer; Th1->Th2 deviation", "approved",
         "subcutaneous", 0.05, "peripheral", ["immunomodulation"], "short",
         ["injection_site"],
         {"infection": 0.1},
         0.35, ["adaptive_immunity"],
         {"IFNG": -0.35, "IL17A": -0.25, "TNF": -0.20, "IL10": 0.55, "TGFB1": 0.45,
          "FOXP3": 0.40, "TBX21": -0.30, "BDNF": 0.30}),

    drug("Interferon beta-1a", "type_I_interferon", "Type I interferon immunomodulation", "approved",
         "intramuscular", 0.10, "peripheral", ["immunomodulation"], "short",
         ["flu_like", "hepatotoxicity"],
         {"infection": 0.15, "hepatic": 0.5},
         0.30, ["adaptive_immunity", "blood_brain_barrier"],
         {"IL17A": -0.35, "IL23A": -0.30, "IFNG": -0.25, "IL10": 0.40, "MMP9": -0.45,
          "VCAM1": -0.25, "CLDN5": 0.30, "STAT1": 0.35, "TNF": -0.20}),

    drug("Interferon beta-1b", "type_I_interferon", "Type I interferon immunomodulation", "approved",
         "subcutaneous", 0.10, "peripheral", ["immunomodulation"], "short",
         ["flu_like", "hepatotoxicity"],
         {"infection": 0.15, "hepatic": 0.5},
         0.30, ["adaptive_immunity", "blood_brain_barrier"],
         {"IL17A": -0.33, "IL23A": -0.28, "IFNG": -0.25, "IL10": 0.38, "MMP9": -0.42,
          "VCAM1": -0.25, "CLDN5": 0.28, "STAT1": 0.33}),

    drug("Peginterferon beta-1a", "type_I_interferon", "Pegylated type I interferon", "approved",
         "subcutaneous", 0.10, "peripheral", ["immunomodulation"], "intermediate",
         ["flu_like", "hepatotoxicity"],
         {"infection": 0.15, "hepatic": 0.45},
         0.30, ["adaptive_immunity", "blood_brain_barrier"],
         {"IL17A": -0.34, "IL23A": -0.29, "IFNG": -0.25, "IL10": 0.39, "MMP9": -0.43,
          "CLDN5": 0.29, "STAT1": 0.34}),

    drug("Mitoxantrone", "topoisomerase_inhibitor", "Type II topoisomerase inhibition", "approved",
         "infusion", 0.20, "peripheral", ["immunomodulation"], "long",
         ["broad_immunosuppression", "cardiotoxicity"],
         {"infection": 0.6, "malignancy": 0.8, "cardiac": 0.9, "teratogenicity": 0.8},
         0.35, ["adaptive_immunity", "B_cell_immunity"],
         {"IL2RA": -0.40, "MS4A1": -0.40, "IFNG": -0.35, "TNF": -0.30, "CD40": -0.25,
          "CASP3": 0.30}, indication="SPMS",
         note="Cumulative-dose-limited; retained mainly as a mechanistic anchor."),

    drug("Natalizumab", "alpha4_integrin_antibody", "Anti-VLA-4 blockade of CNS entry", "approved",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["pml_risk", "infusion_reaction"],
         {"infection": 0.9, "autoimmunity": 0.2},
         0.20, ["lymphocyte_trafficking", "blood_brain_barrier"],
         {"ITGA4": -0.90, "ITGB1": -0.55, "VCAM1": -0.50, "ICAM1": -0.30, "MMP9": -0.35,
          "CCL2": -0.25, "CXCL10": -0.30, "CLDN5": 0.25},
         note="PML risk drives the dominant safety burden, JCV-status dependent."),

    drug("Ocrelizumab", "CD20_depletion", "Anti-CD20 B-cell depletion", "approved",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "long",
         ["b_cell_depletion", "hypogammaglobulinemia"],
         {"infection": 0.6, "malignancy": 0.3},
         0.20, ["B_cell_immunity", "adaptive_immunity"],
         {"MS4A1": -0.90, "CD19": -0.75, "CD27": -0.55, "CXCL13": -0.55, "CD79A": -0.60,
          "CD79B": -0.58, "POU2AF1": -0.50, "AICDA": -0.40, "CD40": -0.25},
         indication="RRMS/PPMS"),

    drug("Ofatumumab", "CD20_depletion", "Anti-CD20, subcutaneous self-administered", "approved",
         "subcutaneous", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["b_cell_depletion", "hypogammaglobulinemia"],
         {"infection": 0.5, "malignancy": 0.25},
         0.20, ["B_cell_immunity", "adaptive_immunity"],
         {"MS4A1": -0.88, "CD19": -0.72, "CD27": -0.52, "CXCL13": -0.52, "CD79A": -0.58,
          "CD79B": -0.55, "POU2AF1": -0.48, "AICDA": -0.38}),

    drug("Ublituximab", "CD20_depletion", "Glycoengineered anti-CD20", "approved",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["b_cell_depletion", "hypogammaglobulinemia"],
         {"infection": 0.5, "malignancy": 0.25},
         0.25, ["B_cell_immunity", "adaptive_immunity"],
         {"MS4A1": -0.90, "CD19": -0.74, "CD27": -0.54, "CXCL13": -0.53, "CD79A": -0.59,
          "CD79B": -0.56, "POU2AF1": -0.49}),

    drug("Rituximab", "CD20_depletion", "Anti-CD20, off-label in MS", "approved",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "long",
         ["b_cell_depletion", "hypogammaglobulinemia"],
         {"infection": 0.55, "malignancy": 0.3},
         0.25, ["B_cell_immunity", "adaptive_immunity"],
         {"MS4A1": -0.85, "CD19": -0.70, "CD27": -0.50, "CXCL13": -0.50, "CD79A": -0.55,
          "CD79B": -0.52, "POU2AF1": -0.45}, indication="off_label",
         note="Widely used off-label; included for mechanism-class breadth."),

    drug("Teriflunomide", "DHODH_inhibitor", "Pyrimidine synthesis blockade", "approved",
         "oral", 0.40, "both", ["immunomodulation", "metabolic_repair"], "long",
         ["hepatotoxicity", "teratogen"],
         {"infection": 0.3, "hepatic": 0.6, "teratogenicity": 0.95},
         0.25, ["adaptive_immunity", "metabolism"],
         {"DHODH": -0.85, "IFNG": -0.30, "IL17A": -0.30, "IL2RA": -0.30, "TNF": -0.20,
          "MTOR": -0.20, "TBX21": -0.25}),
]

# ---------------------------------------------------------------------------
# Tier 2 — late-stage clinical candidates (phase 2/3)
# ---------------------------------------------------------------------------
DRUGS += [
    drug("Tolebrutinib", "BTK_inhibitor", "CNS-penetrant covalent BTK inhibitor", "phase_3",
         "oral", 0.80, "both", ["immunomodulation", "cns_innate"], "short",
         ["b_cell_modulation", "hepatotoxicity"],
         {"infection": 0.35, "hepatic": 0.7},
         0.35, ["B_cell_immunity", "microglial_activation", "innate_immunity"],
         {"BTK": -0.85, "CD79A": -0.40, "CD79B": -0.38, "MS4A1": -0.20, "AIF1": -0.45,
          "CD68": -0.40, "TYROBP": -0.35, "IL1B": -0.35, "TNF": -0.30, "NLRP3": -0.30},
         indication="SPMS/nrSPMS",
         note="CNS penetration is the differentiating claim; hepatic signal is dose-limiting."),

    drug("Fenebrutinib", "BTK_inhibitor", "Non-covalent reversible BTK inhibitor", "phase_3",
         "oral", 0.60, "both", ["immunomodulation", "cns_innate"], "short",
         ["b_cell_modulation", "hepatotoxicity"],
         {"infection": 0.35, "hepatic": 0.5},
         0.35, ["B_cell_immunity", "microglial_activation"],
         {"BTK": -0.80, "CD79A": -0.38, "CD79B": -0.36, "AIF1": -0.35, "CD68": -0.32,
          "IL1B": -0.30, "TNF": -0.25}),

    drug("Evobrutinib", "BTK_inhibitor", "Covalent BTK inhibitor", "phase_3",
         "oral", 0.45, "both", ["immunomodulation"], "short",
         ["b_cell_modulation", "hepatotoxicity"],
         {"infection": 0.3, "hepatic": 0.55},
         0.40, ["B_cell_immunity", "innate_immunity"],
         {"BTK": -0.78, "CD79A": -0.35, "CD79B": -0.33, "AIF1": -0.25, "IL1B": -0.25},
         note="Phase 3 did not meet its primary endpoint; retained as a mechanism control."),

    drug("Remibrutinib", "BTK_inhibitor", "Covalent BTK inhibitor", "phase_3",
         "oral", 0.50, "both", ["immunomodulation"], "short",
         ["b_cell_modulation", "hepatotoxicity"],
         {"infection": 0.3, "hepatic": 0.45},
         0.45, ["B_cell_immunity", "innate_immunity"],
         {"BTK": -0.80, "CD79A": -0.36, "CD79B": -0.34, "IL1B": -0.25, "TNF": -0.22}),

    drug("Orelabrutinib", "BTK_inhibitor", "Covalent BTK inhibitor", "phase_2",
         "oral", 0.55, "both", ["immunomodulation", "cns_innate"], "short",
         ["b_cell_modulation"],
         {"infection": 0.3, "hepatic": 0.4},
         0.45, ["B_cell_immunity", "microglial_activation"],
         {"BTK": -0.78, "CD79A": -0.35, "AIF1": -0.30, "CD68": -0.28, "IL1B": -0.25}),

    drug("Frexalimab", "CD40L_antagonist", "Anti-CD40L costimulation blockade", "phase_3",
         "infusion", 0.10, "peripheral", ["immunomodulation"], "intermediate",
         ["costimulation_blockade"],
         {"infection": 0.4},
         0.40, ["adaptive_immunity", "B_cell_immunity"],
         {"CD40LG": -0.85, "CD40": -0.70, "CD27": -0.30, "AICDA": -0.45, "CXCL13": -0.35,
          "IFNG": -0.30, "IL17A": -0.25, "POU2AF1": -0.30}),

    drug("Vidofludimus calcium", "DHODH_inhibitor", "Second-generation DHODH inhibitor", "phase_3",
         "oral", 0.40, "both", ["immunomodulation", "metabolic_repair"], "intermediate",
         ["hepatotoxicity"],
         {"infection": 0.25, "hepatic": 0.35, "teratogenicity": 0.7},
         0.40, ["adaptive_immunity", "metabolism"],
         {"DHODH": -0.80, "IL17A": -0.35, "IFNG": -0.30, "RORC": -0.30, "TNF": -0.20}),

    drug("Ibudilast", "PDE_inhibitor", "PDE4/10 and MIF inhibition, glial modulation", "phase_2",
         "oral", 0.70, "cns", ["cns_innate", "neuroprotection"], "short",
         ["gi_intolerance"],
         {"infection": 0.1, "hepatic": 0.2},
         0.45, ["microglial_activation", "innate_immunity", "neuroprotection"],
         {"AIF1": -0.50, "CD68": -0.45, "TYROBP": -0.35, "TNF": -0.40, "IL1B": -0.40,
          "NLRP3": -0.35, "NEFL": -0.30, "CASP3": -0.25, "BDNF": 0.30},
         indication="progressive_MS"),

    drug("Masitinib", "tyrosine_kinase_inhibitor", "c-Kit/CSF1R mast cell and microglial TKI", "phase_3",
         "oral", 0.50, "both", ["cns_innate", "immunomodulation"], "intermediate",
         ["neutropenia", "hepatotoxicity"],
         {"infection": 0.35, "hepatic": 0.45, "cardiac": 0.2},
         0.50, ["microglial_activation", "innate_immunity"],
         {"CSF1R": -0.65, "AIF1": -0.40, "CD68": -0.40, "TYROBP": -0.35, "TNF": -0.30,
          "IL1B": -0.28, "ITGAM": -0.30}, indication="progressive_MS"),

    drug("Temelimab", "anti_HERV_W", "Anti-HERV-W envelope monoclonal antibody", "phase_2",
         "infusion", 0.15, "both", ["cns_innate", "remyelination"], "intermediate",
         ["biologic_immunogenicity"],
         {"infection": 0.1},
         0.60, ["microglial_activation", "remyelination"],
         {"TLR4": -0.55, "AIF1": -0.35, "CD68": -0.30, "TNF": -0.25, "MBP": 0.30,
          "PLP1": 0.25, "OLIG2": 0.25, "CSPG4": -0.25},
         note="Targets a putative viral driver; efficacy evidence remains preliminary."),

    drug("Opicinumab", "anti_LINGO1", "Anti-LINGO-1 remyelination antibody", "phase_2",
         "infusion", 0.20, "cns", ["remyelination"], "intermediate",
         ["biologic_immunogenicity"],
         {"infection": 0.05},
         0.60, ["remyelination", "neuroprotection"],
         {"LINGO1": -0.85, "MBP": 0.50, "PLP1": 0.45, "MAG": 0.40, "MYRF": 0.45,
          "OLIG1": 0.35, "OLIG2": 0.35, "SOX10": 0.30, "NEFL": -0.20},
         indication="remyelination",
         note="Phase 2 endpoints were not met; retained as a remyelination mechanism anchor."),

    drug("Elezanumab", "anti_RGMa", "Anti-RGMa neurorepair antibody", "phase_2",
         "infusion", 0.20, "cns", ["remyelination", "neuroprotection"], "intermediate",
         ["biologic_immunogenicity"],
         {"infection": 0.05},
         0.60, ["remyelination", "neuroprotection"],
         {"NEFL": -0.35, "MBP": 0.35, "PLP1": 0.30, "OLIG2": 0.30, "BDNF": 0.35,
          "NTRK2": 0.30, "CASP3": -0.25}, indication="progressive_MS"),

    drug("Clemastine", "H1_antihistamine", "M1 muscarinic antagonism promoting OPC differentiation", "phase_2",
         "oral", 0.75, "cns", ["remyelination"], "short",
         ["anticholinergic", "sedation"],
         {"infection": 0.05},
         0.50, ["remyelination"],
         {"MBP": 0.55, "PLP1": 0.50, "MAG": 0.40, "CNP": 0.40, "MYRF": 0.45,
          "SOX10": 0.35, "OLIG1": 0.35, "OLIG2": 0.35, "GPR17": -0.30},
         indication="remyelination"),

    drug("Bexarotene", "RXR_agonist", "Retinoid X receptor agonist promoting remyelination", "phase_2",
         "oral", 0.55, "both", ["remyelination"], "intermediate",
         ["hypothyroidism", "hyperlipidemia"],
         {"teratogenicity": 0.9, "hepatic": 0.35},
         0.55, ["remyelination", "metabolism"],
         {"MBP": 0.45, "PLP1": 0.40, "MYRF": 0.40, "SOX10": 0.35, "OLIG2": 0.30,
          "LDLR": 0.30, "CSPG4": -0.25}, indication="remyelination",
         note="Remyelination signal offset by a poor tolerability profile in trial."),

    drug("Simvastatin", "HMGCR_inhibitor", "HMG-CoA reductase inhibition", "phase_3",
         "oral", 0.35, "both", ["metabolic_repair", "neuroprotection"], "short",
         ["myopathy"],
         {"hepatic": 0.3, "teratogenicity": 0.6},
         0.45, ["metabolism", "neuroprotection", "innate_immunity"],
         {"HMGCR": -0.80, "LDLR": -0.40, "TNF": -0.25, "IL1B": -0.22, "NEFL": -0.30,
          "MTOR": -0.20, "NOS2": -0.25}, indication="SPMS"),

    drug("Metformin", "AMPK_activator", "AMPK activation, metabolic reprogramming", "phase_2",
         "oral", 0.40, "both", ["metabolic_repair", "remyelination"], "short",
         ["gi_intolerance", "lactic_acidosis"],
         {"hepatic": 0.1},
         0.50, ["metabolism", "remyelination", "adaptive_immunity"],
         {"PRKAA1": 0.75, "MTOR": -0.50, "PPARGC1A": 0.45, "TFAM": 0.35, "SIRT1": 0.40,
          "IL17A": -0.25, "FOXP3": 0.30, "MBP": 0.25, "OLIG2": 0.25}),

    drug("Alpha-lipoic acid", "antioxidant", "Mitochondrial antioxidant cofactor", "phase_2",
         "oral", 0.50, "both", ["metabolic_repair", "neuroprotection"], "short",
         ["gi_intolerance"],
         {},
         0.55, ["oxidative_stress", "metabolism", "neuroprotection"],
         {"NFE2L2": 0.45, "GPX1": 0.40, "GSR": 0.40, "CAT": 0.35, "SOD1": 0.35,
          "TXN": 0.35, "NEFL": -0.25, "MMP9": -0.25}, indication="SPMS"),

    drug("High-dose biotin", "metabolic_cofactor", "Carboxylase cofactor, myelin lipid synthesis", "phase_3",
         "oral", 0.60, "cns", ["metabolic_repair", "remyelination"], "short",
         ["assay_interference"],
         {},
         0.65, ["metabolism", "remyelination"],
         {"PPARGC1A": 0.35, "NDUFA1": 0.35, "ATP5F1A": 0.35, "SLC2A1": 0.30,
          "MBP": 0.30, "PLP1": 0.25}, indication="progressive_MS",
         note="Phase 3 (SPI2) was negative; retained as a metabolic mechanism control."),

    drug("Laquinimod", "immunomodulator", "Aryl hydrocarbon receptor-linked immunomodulation", "phase_3",
         "oral", 0.55, "both", ["immunomodulation", "cns_innate", "neuroprotection"], "intermediate",
         ["hepatotoxicity"],
         {"hepatic": 0.5, "cardiac": 0.4, "teratogenicity": 0.6},
         0.50, ["adaptive_immunity", "microglial_activation", "neuroprotection"],
         {"AIF1": -0.40, "CD68": -0.35, "TNF": -0.30, "IL17A": -0.30, "IL10": 0.30,
          "NEFL": -0.30, "BDNF": 0.35, "CASP3": -0.25},
         note="Development halted at higher doses for cardiovascular findings."),

    drug("Firategrast", "alpha4_integrin_small_molecule", "Oral reversible alpha4 integrin antagonist", "phase_2",
         "oral", 0.30, "peripheral", ["immunomodulation"], "short",
         ["integrin_blockade"],
         {"infection": 0.5},
         0.55, ["lymphocyte_trafficking", "blood_brain_barrier"],
         {"ITGA4": -0.70, "ITGB1": -0.40, "VCAM1": -0.35, "MMP9": -0.25}),

    drug("Inebilizumab", "CD19_depletion", "Anti-CD19 B and plasmablast depletion", "phase_2",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "long",
         ["b_cell_depletion", "hypogammaglobulinemia"],
         {"infection": 0.6, "malignancy": 0.25},
         0.35, ["B_cell_immunity", "adaptive_immunity"],
         {"CD19": -0.90, "MS4A1": -0.60, "CD27": -0.60, "POU2AF1": -0.55, "AICDA": -0.45,
          "CXCL13": -0.45, "TNFSF13B": -0.25}, indication="off_label",
         note="Approved in NMOSD, not MS; included for CD19-vs-CD20 mechanism contrast."),

    drug("Tocilizumab", "IL6R_antagonist", "Anti-IL-6 receptor blockade", "phase_2",
         "infusion", 0.10, "peripheral", ["immunomodulation"], "intermediate",
         ["cytokine_blockade", "hepatotoxicity"],
         {"infection": 0.6, "hepatic": 0.4},
         0.45, ["innate_immunity", "adaptive_immunity"],
         {"IL6": -0.85, "STAT3": -0.60, "RORC": -0.40, "IL17A": -0.40, "TNF": -0.20,
          "C3": -0.20}, indication="off_label"),

    drug("Natalizumab-biosimilar", "alpha4_integrin_antibody", "Biosimilar anti-VLA-4", "phase_3",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["pml_risk", "infusion_reaction"],
         {"infection": 0.9},
         0.30, ["lymphocyte_trafficking", "blood_brain_barrier"],
         {"ITGA4": -0.88, "ITGB1": -0.54, "VCAM1": -0.48, "ICAM1": -0.29, "MMP9": -0.34},
         note="Mechanistic duplicate of natalizumab; acts as a positive control pair."),
]

# ---------------------------------------------------------------------------
# Tier 3 — repurposing and preclinical mechanism-diverse candidates
# ---------------------------------------------------------------------------
DRUGS += [
    drug("Minocycline", "tetracycline_microglial", "Microglial deactivation, MMP inhibition", "phase_2",
         "oral", 0.70, "cns", ["cns_innate", "neuroprotection"], "short",
         ["photosensitivity", "gi_intolerance"],
         {"infection": 0.1, "hepatic": 0.25, "autoimmunity": 0.3},
         0.45, ["microglial_activation", "blood_brain_barrier", "innate_immunity"],
         {"AIF1": -0.55, "CD68": -0.50, "ITGAM": -0.40, "MMP9": -0.60, "TNF": -0.35,
          "IL1B": -0.35, "NOS2": -0.40, "CASP3": -0.30, "CLDN5": 0.25}),

    drug("Hydroxychloroquine", "lysosomotropic_agent", "TLR/lysosomal modulation of microglia", "phase_2",
         "oral", 0.55, "both", ["cns_innate", "immunomodulation"], "long",
         ["retinopathy", "cardiac"],
         {"ocular": 0.6, "cardiac": 0.3},
         0.55, ["microglial_activation", "innate_immunity"],
         {"TLR4": -0.55, "AIF1": -0.40, "CD68": -0.35, "IL1B": -0.30, "TNF": -0.25,
          "NLRP3": -0.30}, indication="PPMS"),

    drug("Riluzole", "glutamate_modulator", "Glutamate release inhibition, Na channel block", "phase_2",
         "oral", 0.65, "cns", ["neuroprotection"], "short",
         ["hepatotoxicity"],
         {"hepatic": 0.5},
         0.50, ["neuroprotection"],
         {"GRIN1": -0.55, "GRIA1": -0.50, "NEFL": -0.30, "CASP3": -0.30, "BAX": -0.25,
          "BCL2": 0.25}),

    drug("Memantine", "NMDA_antagonist", "Uncompetitive NMDA receptor antagonism", "phase_2",
         "oral", 0.75, "cns", ["neuroprotection"], "intermediate",
         ["cns_adverse_effects"],
         {},
         0.55, ["neuroprotection"],
         {"GRIN1": -0.70, "GRIA1": -0.35, "CASP3": -0.25, "NEFL": -0.25}),

    drug("Lamotrigine", "sodium_channel_blocker", "Voltage-gated sodium channel blockade", "phase_2",
         "oral", 0.70, "cns", ["neuroprotection"], "short",
         ["rash_sjs"],
         {"hepatic": 0.2},
         0.55, ["neuroprotection"],
         {"NEFL": -0.35, "NEFH": -0.30, "CASP3": -0.25, "GRIA1": -0.25},
         indication="SPMS",
         note="Phase 2 in SPMS did not show benefit on the primary endpoint."),

    drug("Phenytoin", "sodium_channel_blocker", "Sodium channel blockade, axonal protection", "phase_2",
         "oral", 0.70, "cns", ["neuroprotection"], "intermediate",
         ["rash_sjs", "hepatotoxicity"],
         {"hepatic": 0.4, "teratogenicity": 0.8},
         0.55, ["neuroprotection"],
         {"NEFL": -0.35, "NEFH": -0.32, "CASP3": -0.22}, indication="optic_neuritis"),

    drug("Oxcarbazepine", "sodium_channel_blocker", "Sodium channel blockade", "phase_2",
         "oral", 0.70, "cns", ["neuroprotection"], "short",
         ["hyponatremia", "rash_sjs"],
         {"hepatic": 0.2, "teratogenicity": 0.6},
         0.60, ["neuroprotection"],
         {"NEFL": -0.33, "NEFH": -0.30, "GRIA1": -0.22}, indication="SPMS"),

    drug("Amiloride", "ASIC1_blocker", "Acid-sensing ion channel 1 blockade", "phase_2",
         "oral", 0.45, "cns", ["neuroprotection"], "short",
         ["hyperkalemia"],
         {"cardiac": 0.25},
         0.60, ["neuroprotection"],
         {"NEFL": -0.30, "CASP3": -0.25, "BAX": -0.22, "BCL2": 0.20}, indication="SPMS"),

    drug("Idebenone", "mitochondrial_agent", "Short-chain quinone, electron transport support", "phase_2",
         "oral", 0.55, "cns", ["metabolic_repair", "neuroprotection"], "short",
         ["gi_intolerance"],
         {"hepatic": 0.15},
         0.60, ["metabolism", "oxidative_stress", "neuroprotection"],
         {"NDUFA1": 0.45, "ATP5F1A": 0.45, "TFAM": 0.35, "PPARGC1A": 0.35, "SOD2": 0.30,
          "NEFL": -0.25}, indication="PPMS"),

    drug("Nicotinamide riboside", "NAD_precursor", "NAD+ precursor supporting axonal metabolism", "phase_2",
         "oral", 0.50, "both", ["metabolic_repair", "neuroprotection"], "short",
         [],
         {},
         0.65, ["metabolism", "neuroprotection"],
         {"SIRT1": 0.55, "PPARGC1A": 0.45, "TFAM": 0.40, "NDUFA1": 0.35, "NEFL": -0.20}),

    drug("Resveratrol", "SIRT1_activator", "Sirtuin activation, antioxidant signalling", "preclinical",
         "oral", 0.40, "both", ["metabolic_repair", "neuroprotection"], "short",
         ["gi_intolerance"],
         {},
         0.75, ["metabolism", "oxidative_stress"],
         {"SIRT1": 0.60, "NFE2L2": 0.40, "PPARGC1A": 0.35, "TNF": -0.25, "IL1B": -0.22},
         note="Preclinical MS evidence is mixed, including reports of harm in EAE."),

    drug("Pioglitazone", "PPARG_agonist", "PPAR-gamma agonism, metabolic and glial effects", "phase_2",
         "oral", 0.40, "both", ["metabolic_repair", "cns_innate"], "intermediate",
         ["fluid_retention", "fracture_risk"],
         {"cardiac": 0.4, "hepatic": 0.25, "malignancy": 0.2},
         0.60, ["metabolism", "microglial_activation"],
         {"PPARGC1A": 0.40, "PRKAA1": 0.30, "AIF1": -0.35, "CD68": -0.30, "TNF": -0.30,
          "IL1B": -0.28, "NOS2": -0.30}),

    drug("Fasudil", "ROCK_inhibitor", "Rho-kinase inhibition, axonal outgrowth", "preclinical",
         "infusion", 0.50, "cns", ["remyelination", "neuroprotection"], "short",
         ["hypotension"],
         {"cardiac": 0.3},
         0.75, ["remyelination", "neuroprotection"],
         {"MBP": 0.35, "PLP1": 0.30, "OLIG2": 0.30, "CSPG4": -0.35, "NEFL": -0.25,
          "BDNF": 0.25}),

    drug("Quetiapine", "atypical_antipsychotic", "Repurposed OPC differentiation signal", "preclinical",
         "oral", 0.75, "cns", ["remyelination"], "short",
         ["sedation", "metabolic_syndrome"],
         {"cardiac": 0.3},
         0.80, ["remyelination"],
         {"MBP": 0.40, "PLP1": 0.35, "CNP": 0.30, "OLIG1": 0.30, "SOX10": 0.25},
         note="Rodent remyelination signal only; sedation limits translatability."),

    drug("Miconazole", "antifungal_repurposed", "Repurposed OPC differentiation enhancer", "preclinical",
         "topical", 0.30, "cns", ["remyelination"], "short",
         ["hepatotoxicity"],
         {"hepatic": 0.3},
         0.85, ["remyelination"],
         {"MBP": 0.40, "PLP1": 0.38, "MYRF": 0.35, "OLIG1": 0.30, "SOX10": 0.30},
         note="Systemic exposure and CNS delivery are unresolved barriers."),

    drug("Clobetasol", "corticosteroid_repurposed", "Repurposed glucocorticoid OPC signal", "preclinical",
         "topical", 0.35, "both", ["remyelination", "immunomodulation"], "short",
         ["glucocorticoid_class"],
         {"infection": 0.4},
         0.85, ["remyelination", "adaptive_immunity"],
         {"MBP": 0.38, "PLP1": 0.35, "MYRF": 0.32, "TNF": -0.30, "IL1B": -0.28}),

    drug("Benztropine", "muscarinic_antagonist", "M1/M3 antagonism driving OPC differentiation", "preclinical",
         "oral", 0.70, "cns", ["remyelination"], "short",
         ["anticholinergic"],
         {},
         0.75, ["remyelination"],
         {"MBP": 0.45, "PLP1": 0.40, "MAG": 0.35, "MYRF": 0.38, "SOX10": 0.30,
          "GPR17": -0.28}),

    drug("Methylprednisolone", "corticosteroid", "High-dose glucocorticoid for relapse", "approved",
         "infusion", 0.60, "both", ["immunomodulation"], "short",
         ["glucocorticoid_class", "hyperglycemia"],
         {"infection": 0.5, "cardiac": 0.2},
         0.30, ["adaptive_immunity", "blood_brain_barrier", "innate_immunity"],
         {"TNF": -0.50, "IL1B": -0.50, "IL6": -0.50, "IFNG": -0.40, "IL17A": -0.35,
          "MMP9": -0.45, "CLDN5": 0.35, "ICAM1": -0.35, "VCAM1": -0.30, "IL10": 0.30},
         indication="relapse_treatment",
         note="Relapse rescue only; not a disease-modifying maintenance therapy."),

    drug("Sirolimus", "mTOR_inhibitor", "mTORC1 inhibition, Treg expansion", "phase_2",
         "oral", 0.35, "peripheral", ["immunomodulation", "metabolic_repair"], "long",
         ["broad_immunosuppression", "metabolic"],
         {"infection": 0.6, "hepatic": 0.3, "malignancy": 0.2},
         0.55, ["metabolism", "adaptive_immunity"],
         {"MTOR": -0.85, "FOXP3": 0.45, "IL17A": -0.35, "RORC": -0.35, "IFNG": -0.25,
          "PRKAA1": 0.25}, indication="off_label"),

    drug("Mycophenolate mofetil", "IMPDH_inhibitor", "Inosine monophosphate dehydrogenase inhibition", "phase_2",
         "oral", 0.20, "peripheral", ["immunomodulation"], "short",
         ["broad_immunosuppression", "gi_intolerance"],
         {"infection": 0.6, "teratogenicity": 0.9, "malignancy": 0.3},
         0.50, ["adaptive_immunity", "B_cell_immunity"],
         {"IL2RA": -0.40, "MS4A1": -0.35, "CD19": -0.30, "IFNG": -0.30, "IL17A": -0.28},
         indication="off_label"),

    drug("Azathioprine", "purine_antimetabolite", "Purine synthesis inhibition", "approved",
         "oral", 0.25, "peripheral", ["immunomodulation"], "short",
         ["broad_immunosuppression", "myelosuppression"],
         {"infection": 0.55, "malignancy": 0.5, "hepatic": 0.4, "teratogenicity": 0.7},
         0.45, ["adaptive_immunity"],
         {"IL2RA": -0.35, "IFNG": -0.30, "IL17A": -0.25, "MS4A1": -0.25, "CD40": -0.20},
         indication="off_label"),

    drug("Methotrexate", "DHFR_inhibitor", "Dihydrofolate reductase inhibition", "approved",
         "oral", 0.20, "peripheral", ["immunomodulation"], "short",
         ["broad_immunosuppression", "hepatotoxicity"],
         {"infection": 0.5, "hepatic": 0.6, "teratogenicity": 0.95},
         0.50, ["adaptive_immunity", "metabolism"],
         {"IFNG": -0.30, "IL17A": -0.28, "TNF": -0.30, "IL2RA": -0.25, "IL6": -0.25},
         indication="off_label"),

    drug("Cyclophosphamide", "alkylating_agent", "DNA alkylation, broad lymphoablation", "approved",
         "infusion", 0.25, "peripheral", ["immunomodulation"], "short",
         ["broad_immunosuppression", "urotoxicity"],
         {"infection": 0.8, "malignancy": 0.8, "teratogenicity": 0.95, "cardiac": 0.3},
         0.45, ["adaptive_immunity", "B_cell_immunity"],
         {"IL2RA": -0.45, "MS4A1": -0.45, "CD19": -0.40, "IFNG": -0.40, "IL17A": -0.35,
          "CASP3": 0.30}, indication="off_label",
         note="High-burden rescue therapy; useful as a safety-penalty stress test."),

    drug("Tofacitinib", "JAK_inhibitor", "JAK1/3 inhibition", "preclinical",
         "oral", 0.40, "peripheral", ["immunomodulation"], "short",
         ["broad_immunosuppression", "thrombosis"],
         {"infection": 0.6, "malignancy": 0.4, "cardiac": 0.4, "hepatic": 0.3},
         0.60, ["adaptive_immunity", "innate_immunity"],
         {"STAT1": -0.55, "STAT3": -0.55, "STAT4": -0.50, "IFNG": -0.45, "IL17A": -0.40,
          "IL23A": -0.35, "IL6": -0.30, "IL2RA": -0.35}),

    drug("Baricitinib", "JAK_inhibitor", "JAK1/2 inhibition", "preclinical",
         "oral", 0.45, "both", ["immunomodulation", "cns_innate"], "short",
         ["broad_immunosuppression", "thrombosis"],
         {"infection": 0.6, "malignancy": 0.35, "cardiac": 0.4},
         0.65, ["adaptive_immunity", "microglial_activation"],
         {"STAT1": -0.55, "STAT3": -0.50, "IFNG": -0.45, "IL6": -0.35, "AIF1": -0.30,
          "CD68": -0.25, "IL23A": -0.30}),

    drug("Ustekinumab", "IL12_IL23_antagonist", "Anti-p40 IL-12/IL-23 blockade", "phase_2",
         "subcutaneous", 0.05, "peripheral", ["immunomodulation"], "long",
         ["cytokine_blockade"],
         {"infection": 0.4},
         0.45, ["adaptive_immunity"],
         {"IL12B": -0.85, "IL23A": -0.70, "IL17A": -0.45, "RORC": -0.40, "IFNG": -0.35,
          "TBX21": -0.30}, indication="off_label",
         note="Phase 2 in MS was negative; a valuable negative mechanistic control."),

    drug("Secukinumab", "IL17A_antagonist", "Anti-IL-17A blockade", "phase_2",
         "subcutaneous", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["cytokine_blockade"],
         {"infection": 0.35},
         0.50, ["adaptive_immunity"],
         {"IL17A": -0.90, "RORC": -0.35, "IL23A": -0.25, "CCR6": -0.30},
         indication="off_label"),

    drug("Eculizumab", "complement_C5_inhibitor", "Terminal complement C5 blockade", "preclinical",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["complement_blockade", "meningococcal_risk"],
         {"infection": 0.85},
         0.60, ["innate_immunity"],
         {"C5": -0.90, "C3": -0.30, "C1QA": -0.20, "CFH": 0.20}, indication="off_label"),

    drug("Belimumab", "BAFF_antagonist", "Anti-BAFF B-cell survival blockade", "preclinical",
         "infusion", 0.05, "peripheral", ["immunomodulation"], "long",
         ["b_cell_depletion"],
         {"infection": 0.4},
         0.60, ["B_cell_immunity"],
         {"TNFSF13B": -0.90, "CD19": -0.35, "MS4A1": -0.30, "CXCL13": -0.30,
          "POU2AF1": -0.30}, indication="off_label"),

    drug("Cholecalciferol", "vitamin_D", "Vitamin D receptor immunomodulation", "phase_3",
         "oral", 0.55, "both", ["immunomodulation", "metabolic_repair"], "long",
         ["hypercalcemia"],
         {},
         0.60, ["adaptive_immunity", "metabolism"],
         {"FOXP3": 0.40, "IL10": 0.35, "TGFB1": 0.30, "IL17A": -0.30, "IFNG": -0.25,
          "RORC": -0.25}, note="Adjunct evidence; not established as monotherapy DMT."),

    drug("Epigallocatechin gallate", "polyphenol_antioxidant", "Catechin antioxidant and immune modulation", "phase_2",
         "oral", 0.45, "both", ["metabolic_repair", "neuroprotection"], "short",
         ["hepatotoxicity"],
         {"hepatic": 0.45},
         0.75, ["oxidative_stress", "neuroprotection"],
         {"NFE2L2": 0.40, "HMOX1": 0.35, "CAT": 0.30, "SOD1": 0.30, "TNF": -0.25,
          "NEFL": -0.20}),

    drug("Daclizumab", "CD25_antagonist", "Anti-CD25 (IL-2R alpha) modulation", "approved",
         "subcutaneous", 0.05, "peripheral", ["immunomodulation"], "intermediate",
         ["hepatotoxicity", "autoimmune_encephalitis"],
         {"infection": 0.4, "hepatic": 0.9, "autoimmunity": 0.95},
         0.35, ["adaptive_immunity"],
         {"IL2RA": -0.90, "IFNG": -0.30, "IL17A": -0.30, "FOXP3": -0.25},
         indication="withdrawn",
         note="Withdrawn worldwide for fatal immune hepatitis/encephalitis. Retained as "
              "a safety-penalty positive control: any ranking that promotes it is wrong."),
]


def validate(drugs: list[dict]) -> dict[str, object]:
    """Cross-check the panel against the signature and report coverage stats."""
    import csv

    names = [d["name"] for d in drugs]
    duplicates = {n for n in names if names.count(n) > 1}
    if duplicates:
        raise ValueError(f"Duplicate panel entries: {sorted(duplicates)}")

    with (ROOT / "data/ms_expression_v3.csv").open(newline="") as handle:
        universe = {row["gene"] for row in csv.DictReader(handle)}

    orphans = {}
    for record in drugs:
        missing = sorted(set(record["target_effects"]) - universe)
        if missing:
            orphans[record["name"]] = missing
    if orphans:
        raise ValueError(
            "Target genes absent from data/ms_expression_v3.csv: "
            + json.dumps(orphans, indent=2)
        )

    covered = {g for d in drugs for g in d["target_effects"]}
    return {
        "n_drugs": len(drugs),
        "n_signature_genes": len(universe),
        "n_genes_targeted": len(covered),
        "signature_gene_coverage": round(len(covered) / len(universe), 3),
        "mechanism_classes": len({d["mechanism_class"] for d in drugs}),
        "evidence_tiers": {
            tier: sum(1 for d in drugs if d["evidence_tier"] == tier)
            for tier in ("approved", "phase_3", "phase_2", "preclinical")
        },
    }


def main() -> None:
    for record in DRUGS:
        record["target_family"] = TARGET_FAMILY.get(record["name"], record["mechanism_class"])
    stats = validate(DRUGS)
    payload = {
        "metadata": {
            "panel_version": "3.0.0",
            "generated_by": "tools/curate_ms_panel_v3.py",
            "purpose": (
                "Mechanism-curated candidate panel for research prioritisation. "
                "Not a prescribing, efficacy, or co-administration dataset."
            ),
            "schema": {
                "target_effects": "signed unitless directional hypotheses in [-1,1]; not affinities or doses",
                "cns_penetration": "0-1 ordinal confidence that the agent reaches CNS parenchyma at therapeutic exposure",
                "compartment": "dominant site of action: peripheral | cns | both",
                "therapeutic_axes": list(AXES),
                "safety_burden": "0-1 ordinal, ranking-relevant risk magnitude per domain: " + ", ".join(RISKS),
                "half_life_class": "short | intermediate | long | reconstitution",
                "target_family": (
                    "biological target group; agents sharing one are pharmacodynamically "
                    "redundant regardless of modality and are never ranked as a combination"
                ),
            },
            "sources": SOURCES,
            "curation_rule": (
                "Approved MS disease-modifying therapies, late-stage clinical candidates, and "
                "mechanism-diverse repurposing/preclinical agents. Records marked as negative or "
                "safety controls are retained deliberately and must not be filtered out. "
                "Every record requires independent verification before publication use."
            ),
            "controls": {
                "positive_redundancy": ["Natalizumab", "Natalizumab-biosimilar"],
                "negative_efficacy": ["Ustekinumab", "High-dose biotin", "Opicinumab", "Evobrutinib"],
                "safety_penalty": ["Daclizumab", "Cyclophosphamide", "Mitoxantrone"],
            },
            "statistics": stats,
        },
        "drugs": sorted(DRUGS, key=lambda d: d["name"]),
    }
    out = ROOT / "data/drugs/ms_panel_v3.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {out.relative_to(ROOT)}")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()

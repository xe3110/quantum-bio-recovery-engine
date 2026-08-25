"""Curate and validate the Parkinson's disease candidate panel.

The second disease in the registry, and the reason it exists: the design stack
claimed to be disease-agnostic on the strength of one populated entry. This
panel exercises the claim against a disease with different pathways, different
therapeutic axes, and -- the part that actually tests the abstraction -- a
completely different safety vocabulary. Nothing in Parkinson's care is
constrained by infection or malignancy risk; it is constrained by dyskinesia,
impulse-control disorders, orthostatic hypotension, and psychosis.

As with the MS curation tool, this **fails the build** if any target gene is
absent from the disease signature, so a typo cannot enter the panel silently.

Run: python -m tools.curate_pd_panel
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIGNATURE = ROOT / "data/pd_expression.csv"
OUT = ROOT / "data/drugs/parkinsons_panel_v1.json"

RISK_DOMAINS = (
    "dyskinesia", "impulse_control", "orthostatic_hypotension", "psychiatric",
    "cardiac", "hepatic", "somnolence", "gastrointestinal",
)
AXES = (
    "symptomatic_dopaminergic", "synuclein_proteostasis", "mitochondrial_rescue",
    "neuroinflammation_control", "trophic_support",
)


def d(**kw):
    """Build one panel record, filling the risk vector with declared zeros."""
    burden = {k: 0.0 for k in RISK_DOMAINS}
    burden.update(kw.pop("risk", {}))
    kw["safety_burden"] = burden
    kw.setdefault("target_family", kw["mechanism_class"])
    kw.setdefault("pathways", [])
    kw["relevant_pathways"] = kw["pathways"]
    return kw


DRUGS = [
 # --- symptomatic dopaminergic replacement -----------------------------
 d(name="Levodopa", mechanism_class="dopamine_precursor", primary_moa="Dopamine precursor restoring striatal dopamine",
   evidence_tier="approved", indication="all_stages", route="oral", cns_penetration=0.8, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.1,
   pathways=["dopaminergic_neurotransmission"],
   risk={"dyskinesia":0.9,"psychiatric":0.4,"orthostatic_hypotension":0.4,"gastrointestinal":0.5,"somnolence":0.3},
   target_effects={"TH":0.3,"DDC":0.5,"SLC18A2":0.35,"DRD1":0.7,"DRD2":0.7,"SLC6A3":0.2}),
 d(name="Carbidopa", mechanism_class="peripheral_AADC_inhibitor", primary_moa="Peripheral decarboxylase inhibition sparing levodopa",
   evidence_tier="approved", indication="adjunct", route="oral", cns_penetration=0.05, compartment="peripheral",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.1, pathways=["dopaminergic_neurotransmission"],
   risk={"gastrointestinal":0.2,"orthostatic_hypotension":0.2},
   target_effects={"DDC":-0.6}),
 d(name="Pramipexole", mechanism_class="D2_D3_agonist", primary_moa="Non-ergot D2/D3 receptor agonist",
   evidence_tier="approved", indication="all_stages", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="intermediate",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.15, pathways=["dopaminergic_neurotransmission"],
   risk={"impulse_control":0.8,"psychiatric":0.5,"somnolence":0.7,"orthostatic_hypotension":0.5,"dyskinesia":0.3},
   target_effects={"DRD2":0.85,"DRD1":0.2}),
 d(name="Ropinirole", mechanism_class="D2_D3_agonist", primary_moa="Non-ergot D2/D3 receptor agonist",
   evidence_tier="approved", indication="all_stages", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.15, pathways=["dopaminergic_neurotransmission"],
   risk={"impulse_control":0.8,"psychiatric":0.5,"somnolence":0.7,"orthostatic_hypotension":0.5,"dyskinesia":0.3},
   target_effects={"DRD2":0.85,"DRD1":0.2}),
 d(name="Rotigotine", mechanism_class="D2_D3_agonist", primary_moa="Transdermal non-ergot dopamine agonist",
   evidence_tier="approved", indication="all_stages", route="transdermal", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="intermediate",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.2, pathways=["dopaminergic_neurotransmission"],
   risk={"impulse_control":0.7,"psychiatric":0.4,"somnolence":0.6,"orthostatic_hypotension":0.4},
   target_effects={"DRD2":0.8,"DRD1":0.35}),
 d(name="Apomorphine", mechanism_class="D2_D3_agonist_rescue", primary_moa="Rapid-onset dopamine agonist for off episodes",
   evidence_tier="approved", indication="off_episodes", route="subcutaneous", cns_penetration=0.75, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["dopaminergic_stimulation","injection_reaction"], target_uncertainty=0.2,
   pathways=["dopaminergic_neurotransmission"],
   risk={"orthostatic_hypotension":0.9,"psychiatric":0.6,"gastrointestinal":0.8,"somnolence":0.6,"dyskinesia":0.4},
   target_effects={"DRD2":0.9,"DRD1":0.6}),
 d(name="Pergolide", mechanism_class="ergot_dopamine_agonist", primary_moa="Ergot-derived dopamine agonist, withdrawn for valvulopathy",
   evidence_tier="approved", indication="withdrawn", route="oral", cns_penetration=0.65, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="intermediate",
   safety_classes=["dopaminergic_stimulation","serotonergic_valvulopathy"], target_uncertainty=0.3,
   pathways=["dopaminergic_neurotransmission"],
   risk={"cardiac":0.95,"impulse_control":0.6,"psychiatric":0.5,"somnolence":0.5,"orthostatic_hypotension":0.5},
   target_effects={"DRD2":0.8,"DRD1":0.4}),
 # --- enzymatic modulation of dopamine turnover -------------------------
 d(name="Selegiline", mechanism_class="MAOB_inhibitor", primary_moa="Irreversible MAO-B inhibition",
   evidence_tier="approved", indication="early", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic","mitochondrial_rescue"], half_life_class="intermediate",
   safety_classes=["monoamine_potentiation"], target_uncertainty=0.2,
   pathways=["dopaminergic_neurotransmission","oxidative_stress"],
   risk={"psychiatric":0.4,"orthostatic_hypotension":0.3,"dyskinesia":0.3,"somnolence":0.2},
   target_effects={"MAOB":-0.85,"SLC6A3":0.15,"SOD2":0.2}),
 d(name="Rasagiline", mechanism_class="MAOB_inhibitor", primary_moa="Irreversible MAO-B inhibition",
   evidence_tier="approved", indication="early", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic","mitochondrial_rescue"], half_life_class="intermediate",
   safety_classes=["monoamine_potentiation"], target_uncertainty=0.2,
   pathways=["dopaminergic_neurotransmission","apoptosis"],
   risk={"psychiatric":0.3,"orthostatic_hypotension":0.3,"dyskinesia":0.3},
   target_effects={"MAOB":-0.9,"BAX":-0.25,"BCL2":0.25}),
 d(name="Safinamide", mechanism_class="MAOB_glutamate_modulator", primary_moa="Reversible MAO-B inhibition plus glutamate release modulation",
   evidence_tier="approved", indication="adjunct", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic","neuroinflammation_control"], half_life_class="intermediate",
   safety_classes=["monoamine_potentiation"], target_uncertainty=0.25,
   pathways=["dopaminergic_neurotransmission"],
   risk={"dyskinesia":0.5,"psychiatric":0.3,"hepatic":0.2},
   target_effects={"MAOB":-0.8,"CACNA1D":-0.2,"CASP3":-0.2}),
 d(name="Entacapone", mechanism_class="COMT_inhibitor", primary_moa="Peripheral COMT inhibition extending levodopa exposure",
   evidence_tier="approved", indication="adjunct", route="oral", cns_penetration=0.1, compartment="peripheral",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.15,
   pathways=["dopaminergic_neurotransmission"],
   risk={"dyskinesia":0.5,"gastrointestinal":0.4,"orthostatic_hypotension":0.2},
   target_effects={"COMT":-0.85}),
 d(name="Opicapone", mechanism_class="COMT_inhibitor", primary_moa="Long-acting peripheral COMT inhibition",
   evidence_tier="approved", indication="adjunct", route="oral", cns_penetration=0.1, compartment="peripheral",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="long",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.15,
   pathways=["dopaminergic_neurotransmission"],
   risk={"dyskinesia":0.5,"gastrointestinal":0.3},
   target_effects={"COMT":-0.9}),
 d(name="Tolcapone", mechanism_class="COMT_inhibitor", primary_moa="Central and peripheral COMT inhibition; hepatotoxicity risk",
   evidence_tier="approved", indication="restricted", route="oral", cns_penetration=0.45, compartment="both",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["dopaminergic_stimulation","hepatotoxicity"], target_uncertainty=0.2,
   pathways=["dopaminergic_neurotransmission"],
   risk={"hepatic":0.95,"dyskinesia":0.5,"gastrointestinal":0.4},
   target_effects={"COMT":-0.95}),
 d(name="Istradefylline", mechanism_class="A2A_antagonist", primary_moa="Adenosine A2A receptor antagonism",
   evidence_tier="approved", indication="adjunct", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic","neuroinflammation_control"], half_life_class="long",
   safety_classes=["dopaminergic_stimulation"], target_uncertainty=0.3,
   pathways=["dopaminergic_neurotransmission","neuroinflammation"],
   risk={"dyskinesia":0.5,"psychiatric":0.3,"gastrointestinal":0.2},
   target_effects={"DRD2":0.3,"AIF1":-0.2,"TNF":-0.2}),
 d(name="Amantadine", mechanism_class="NMDA_antagonist", primary_moa="NMDA antagonism reducing levodopa-induced dyskinesia",
   evidence_tier="approved", indication="dyskinesia", route="oral", cns_penetration=0.75, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic","neuroinflammation_control"], half_life_class="intermediate",
   safety_classes=["anticholinergic_burden"], target_uncertainty=0.3,
   pathways=["dopaminergic_neurotransmission","neuroinflammation"],
   risk={"psychiatric":0.6,"orthostatic_hypotension":0.3,"somnolence":0.2},
   target_effects={"DRD2":0.25,"SLC6A3":0.2,"AIF1":-0.25,"IL1B":-0.2}),
 d(name="Trihexyphenidyl", mechanism_class="anticholinergic", primary_moa="Muscarinic antagonism for tremor",
   evidence_tier="approved", indication="tremor", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["symptomatic_dopaminergic"], half_life_class="short",
   safety_classes=["anticholinergic_burden"], target_uncertainty=0.35,
   pathways=["dopaminergic_neurotransmission"],
   risk={"psychiatric":0.7,"somnolence":0.4,"gastrointestinal":0.3},
   target_effects={"DRD2":0.15}),
 # --- disease-modifying candidates: synuclein and lysosome ---------------
 d(name="Ambroxol", mechanism_class="GCase_chaperone", primary_moa="Pharmacological chaperone raising glucocerebrosidase activity",
   evidence_tier="phase_2", indication="disease_modifying", route="oral", cns_penetration=0.5, compartment="cns",
   therapeutic_axes=["synuclein_proteostasis"], half_life_class="short",
   safety_classes=["lysosomal_modulation"], target_uncertainty=0.4,
   pathways=["lysosomal_autophagy","alpha_synuclein_proteostasis"],
   risk={"gastrointestinal":0.3},
   target_effects={"GBA1":0.75,"SNCA":-0.4,"LAMP1":0.35,"SCARB2":0.3,"CTSD":0.25,"TFEB":0.2}),
 d(name="Venglustat", mechanism_class="GCS_inhibitor", primary_moa="Glucosylceramide synthase inhibition reducing substrate load",
   evidence_tier="phase_2", indication="disease_modifying", route="oral", cns_penetration=0.55, compartment="cns",
   therapeutic_axes=["synuclein_proteostasis"], half_life_class="long",
   safety_classes=["lysosomal_modulation"], target_uncertainty=0.45,
   pathways=["lysosomal_autophagy"],
   risk={"gastrointestinal":0.3,"psychiatric":0.2},
   target_effects={"GBA1":0.4,"SNCA":-0.35,"LAMP2":0.25}),
 d(name="Prasinezumab", mechanism_class="anti_synuclein_antibody", primary_moa="Monoclonal antibody against aggregated alpha-synuclein",
   evidence_tier="phase_2", indication="disease_modifying", route="infusion", cns_penetration=0.05, compartment="peripheral",
   therapeutic_axes=["synuclein_proteostasis"], half_life_class="long",
   safety_classes=["infusion_reaction"], target_uncertainty=0.5,
   pathways=["alpha_synuclein_proteostasis"],
   risk={"gastrointestinal":0.1},
   target_effects={"SNCA":-0.6,"SNCAIP":-0.25,"NEFL":-0.2}),
 d(name="Cinpanemab", mechanism_class="anti_synuclein_antibody", primary_moa="Monoclonal antibody against alpha-synuclein; failed phase 2",
   evidence_tier="phase_2", indication="failed", route="infusion", cns_penetration=0.05, compartment="peripheral",
   therapeutic_axes=["synuclein_proteostasis"], half_life_class="long",
   safety_classes=["infusion_reaction"], target_uncertainty=0.55,
   pathways=["alpha_synuclein_proteostasis"],
   risk={"gastrointestinal":0.1},
   target_effects={"SNCA":-0.55,"SNCAIP":-0.2}),
 d(name="Nilotinib", mechanism_class="ABL_inhibitor", primary_moa="c-Abl inhibition promoting autophagic synuclein clearance",
   evidence_tier="phase_2", indication="failed", route="oral", cns_penetration=0.2, compartment="both",
   therapeutic_axes=["synuclein_proteostasis"], half_life_class="intermediate",
   safety_classes=["kinase_inhibition","QT_prolongation"], target_uncertainty=0.5,
   pathways=["alpha_synuclein_proteostasis","lysosomal_autophagy"],
   risk={"cardiac":0.7,"hepatic":0.4,"gastrointestinal":0.4},
   target_effects={"SNCA":-0.4,"MAP1LC3B":0.35,"SQSTM1":-0.3,"TFEB":0.25,"MTOR":-0.2}),
 # --- disease-modifying candidates: kinase and mitochondria --------------
 d(name="LRRK2 inhibitor (BIIB122)", mechanism_class="LRRK2_inhibitor", primary_moa="Selective LRRK2 kinase inhibition",
   evidence_tier="phase_2", indication="disease_modifying", route="oral", cns_penetration=0.65, compartment="cns",
   therapeutic_axes=["synuclein_proteostasis","mitochondrial_rescue"], half_life_class="intermediate",
   safety_classes=["kinase_inhibition","pulmonary_phospholipidosis"], target_uncertainty=0.4,
   pathways=["kinase_signalling","lysosomal_autophagy"],
   risk={"gastrointestinal":0.2,"psychiatric":0.1},
   target_effects={"LRRK2":-0.9,"RAB10":-0.7,"RAB29":-0.5,"SNCA":-0.3,"LAMP1":0.3,"MAP1LC3B":0.3,"VPS35":0.2}),
 d(name="Exenatide", mechanism_class="GLP1_agonist", primary_moa="GLP-1 receptor agonism with neurotrophic and anti-inflammatory effects",
   evidence_tier="phase_3", indication="disease_modifying", route="subcutaneous", cns_penetration=0.3, compartment="both",
   therapeutic_axes=["trophic_support","neuroinflammation_control","mitochondrial_rescue"], half_life_class="intermediate",
   safety_classes=["incretin_effects"], target_uncertainty=0.45,
   pathways=["trophic_support","neuroinflammation","mitochondrial_function"],
   risk={"gastrointestinal":0.6,"psychiatric":0.1},
   target_effects={"BDNF":0.45,"PPARGC1A":0.4,"TNF":-0.35,"IL1B":-0.35,"AIF1":-0.3,"NFKB1":-0.3,"CASP3":-0.25,"MTOR":-0.2}),
 d(name="Ursodeoxycholic acid", mechanism_class="mitochondrial_rescue_agent", primary_moa="Bile acid improving mitochondrial function",
   evidence_tier="phase_2", indication="disease_modifying", route="oral", cns_penetration=0.3, compartment="both",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="intermediate",
   safety_classes=["hepatobiliary"], target_uncertainty=0.5,
   pathways=["mitochondrial_function","apoptosis"],
   risk={"gastrointestinal":0.4,"hepatic":0.2},
   target_effects={"NDUFS1":0.35,"ATP5F1A":0.35,"PPARGC1A":0.3,"BAX":-0.25,"BCL2":0.25}),
 d(name="Deferiprone", mechanism_class="iron_chelator", primary_moa="Brain-penetrant iron chelation; failed phase 2/3 on motor outcome",
   evidence_tier="phase_3", indication="failed", route="oral", cns_penetration=0.55, compartment="cns",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="short",
   safety_classes=["myelosuppression"], target_uncertainty=0.35,
   pathways=["iron_metabolism","oxidative_stress"],
   risk={"gastrointestinal":0.4,"hepatic":0.3,"cardiac":0.2},
   target_effects={"FTH1":-0.6,"FTL":-0.6,"TFRC":-0.4,"ACSL4":-0.3,"GPX4":0.3,"SLC40A1":0.25}),
 d(name="Isradipine", mechanism_class="CaV1_3_blocker", primary_moa="L-type calcium channel blockade; failed phase 3 (STEADY-PD III)",
   evidence_tier="phase_3", indication="failed", route="oral", cns_penetration=0.4, compartment="both",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="short",
   safety_classes=["vasodilation"], target_uncertainty=0.35,
   pathways=["calcium_homeostasis","mitochondrial_function"],
   risk={"orthostatic_hypotension":0.7,"cardiac":0.3},
   target_effects={"CACNA1D":-0.8,"ITPR1":-0.3,"CALB1":0.2}),
 d(name="Inosine", mechanism_class="urate_precursor", primary_moa="Raises serum urate as an antioxidant; failed phase 3 (SURE-PD3)",
   evidence_tier="phase_3", indication="failed", route="oral", cns_penetration=0.4, compartment="both",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="short",
   safety_classes=["urate_elevation"], target_uncertainty=0.5,
   pathways=["oxidative_stress"],
   risk={"gastrointestinal":0.3},
   target_effects={"SOD1":0.25,"CAT":0.25,"GPX4":0.2,"TXN":0.2}),
 d(name="Coenzyme Q10", mechanism_class="mitochondrial_cofactor", primary_moa="Electron transport cofactor; failed phase 3 (QE3)",
   evidence_tier="phase_3", indication="failed", route="oral", cns_penetration=0.25, compartment="both",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="intermediate",
   safety_classes=["nutraceutical"], target_uncertainty=0.55,
   pathways=["mitochondrial_function","oxidative_stress"],
   risk={"gastrointestinal":0.2},
   target_effects={"NDUFS1":0.3,"NDUFA9":0.3,"ATP5F1A":0.25,"SOD2":0.2}),
 d(name="Creatine", mechanism_class="energetic_buffer", primary_moa="Phosphocreatine energy buffering; failed phase 3 (NET-PD LS-1)",
   evidence_tier="phase_3", indication="failed", route="oral", cns_penetration=0.3, compartment="both",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="short",
   safety_classes=["nutraceutical"], target_uncertainty=0.55,
   pathways=["mitochondrial_function"],
   risk={"gastrointestinal":0.3},
   target_effects={"ATP5F1A":0.3,"TFAM":0.2}),
 d(name="Pioglitazone", mechanism_class="PPARG_agonist", primary_moa="PPAR-gamma agonism; failed phase 2 futility (FS-ZONE)",
   evidence_tier="phase_2", indication="failed", route="oral", cns_penetration=0.4, compartment="both",
   therapeutic_axes=["mitochondrial_rescue","neuroinflammation_control"], half_life_class="intermediate",
   safety_classes=["fluid_retention"], target_uncertainty=0.45,
   pathways=["mitochondrial_function","neuroinflammation"],
   risk={"cardiac":0.5,"gastrointestinal":0.2},
   target_effects={"PPARGC1A":0.45,"TFAM":0.3,"TNF":-0.3,"IL1B":-0.3,"NFKB1":-0.3,"AIF1":-0.25}),
 # --- neuroinflammation --------------------------------------------------
 d(name="Minocycline", mechanism_class="tetracycline_microglial", primary_moa="Microglial deactivation; failed phase 2 futility",
   evidence_tier="phase_2", indication="failed", route="oral", cns_penetration=0.7, compartment="cns",
   therapeutic_axes=["neuroinflammation_control"], half_life_class="intermediate",
   safety_classes=["tetracycline_class"], target_uncertainty=0.4,
   pathways=["neuroinflammation"],
   risk={"gastrointestinal":0.4,"psychiatric":0.2},
   target_effects={"AIF1":-0.6,"CD68":-0.5,"IL1B":-0.45,"TNF":-0.4,"NLRP3":-0.35,"CASP1":-0.3,"GFAP":-0.3,"CASP3":-0.3}),
 d(name="NLRP3 inhibitor (inzomelid-class)", mechanism_class="NLRP3_inhibitor", primary_moa="Inflammasome blockade upstream of IL-1beta",
   evidence_tier="phase_2", indication="disease_modifying", route="oral", cns_penetration=0.6, compartment="cns",
   therapeutic_axes=["neuroinflammation_control"], half_life_class="intermediate",
   safety_classes=["innate_immune_suppression"], target_uncertainty=0.45,
   pathways=["neuroinflammation"],
   risk={"gastrointestinal":0.2,"hepatic":0.2},
   target_effects={"NLRP3":-0.85,"CASP1":-0.7,"IL1B":-0.7,"IL6":-0.35,"AIF1":-0.35,"GFAP":-0.25}),
 d(name="Sargramostim", mechanism_class="GM_CSF", primary_moa="GM-CSF driving regulatory T-cell responses",
   evidence_tier="phase_2", indication="disease_modifying", route="subcutaneous", cns_penetration=0.1, compartment="peripheral",
   therapeutic_axes=["neuroinflammation_control"], half_life_class="short",
   safety_classes=["immune_stimulation"], target_uncertainty=0.55,
   pathways=["neuroinflammation"],
   risk={"gastrointestinal":0.3,"cardiac":0.2},
   target_effects={"TNF":-0.3,"IL1B":-0.3,"HLA-DRA":-0.25,"CCL2":-0.25}),
 # --- trophic support ----------------------------------------------------
 d(name="GDNF (intraputaminal)", mechanism_class="GDNF_delivery", primary_moa="Direct glial cell line-derived neurotrophic factor infusion",
   evidence_tier="phase_2", indication="failed", route="intraputaminal", cns_penetration=1.0, compartment="cns",
   therapeutic_axes=["trophic_support"], half_life_class="short",
   safety_classes=["neurosurgical"], target_uncertainty=0.5,
   pathways=["trophic_support","dopaminergic_neurotransmission"],
   risk={"psychiatric":0.3,"gastrointestinal":0.1},
   target_effects={"GDNF":0.9,"RET":0.7,"TH":0.5,"SLC6A3":0.4,"NR4A2":0.35,"BDNF":0.25}),
 d(name="Nicotinamide riboside", mechanism_class="NAD_precursor", primary_moa="NAD+ precursor supporting mitochondrial metabolism",
   evidence_tier="phase_2", indication="disease_modifying", route="oral", cns_penetration=0.5, compartment="both",
   therapeutic_axes=["mitochondrial_rescue"], half_life_class="short",
   safety_classes=["nutraceutical"], target_uncertainty=0.5,
   pathways=["mitochondrial_function"],
   risk={"gastrointestinal":0.2},
   target_effects={"PPARGC1A":0.4,"TFAM":0.35,"NDUFS1":0.3,"ATP5F1A":0.3,"SOD2":0.25}),
]

CONTROLS = {
    "positive_redundancy": ["Pramipexole", "Ropinirole"],
    "negative_efficacy": ["Creatine", "Coenzyme Q10", "Isradipine", "Inosine", "Cinpanemab"],
    "safety_penalty": ["Tolcapone", "Pergolide", "Apomorphine"],
}


def main() -> None:
    with SIGNATURE.open(newline="") as handle:
        signature_genes = {row["gene"] for row in csv.DictReader(handle)}

    unknown = sorted({
        gene for drug in DRUGS for gene in drug["target_effects"]
        if gene not in signature_genes
    })
    if unknown:
        raise SystemExit(
            f"Panel references genes absent from the PD signature: {unknown}\n"
            "Add them to data/pd_expression.csv or correct the symbol."
        )
    for drug in DRUGS:
        stray = set(drug["safety_burden"]) - set(RISK_DOMAINS)
        if stray:
            raise SystemExit(f"{drug['name']} declares unknown risk domains: {sorted(stray)}")
        stray_axes = set(drug["therapeutic_axes"]) - set(AXES)
        if stray_axes:
            raise SystemExit(f"{drug['name']} declares unknown axes: {sorted(stray_axes)}")
        for gene, value in drug["target_effects"].items():
            if not -1.0 <= value <= 1.0:
                raise SystemExit(f"{drug['name']} effect on {gene} outside [-1,1]: {value}")
    names = [d["name"] for d in DRUGS]
    if len(names) != len(set(names)):
        raise SystemExit("Duplicate drug name in panel")
    for group, members in CONTROLS.items():
        missing = [m for m in members if m not in names]
        if missing:
            raise SystemExit(f"Control group {group} names absent agents: {missing}")

    targeted = {g for d in DRUGS for g in d["target_effects"]}
    payload = {
        "metadata": {
            "panel_version": "1.0.0",
            "generated_by": "tools/curate_pd_panel.py",
            "disease": "parkinsons",
            "purpose": "Mechanism-curated candidate panel for research prioritisation. Not a prescribing, efficacy, or co-administration dataset.",
            "curation_rule": "Approved symptomatic therapies plus disease-modifying candidates, weighted deliberately toward agents that failed in phase 2/3. Parkinson's has an unusually well-documented record of neuroprotection failures, which makes it a strong source of negative controls.",
            "controls": CONTROLS,
            "statistics": {
                "n_drugs": len(DRUGS),
                "n_signature_genes": len(signature_genes),
                "n_genes_targeted": len(targeted),
                "signature_gene_coverage": round(len(targeted) / len(signature_genes), 3),
                "mechanism_classes": len({d["mechanism_class"] for d in DRUGS}),
                "evidence_tiers": {
                    tier: sum(1 for d in DRUGS if d["evidence_tier"] == tier)
                    for tier in ("approved", "phase_3", "phase_2", "preclinical")
                },
            },
            "caveat": "Directional target effects are curated hypotheses, not measured pharmacology. Every record requires independent verification before publication use.",
        },
        "drugs": DRUGS,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    stats = payload["metadata"]["statistics"]
    print(f"Wrote {OUT.relative_to(ROOT)}")
    print(f"  {stats['n_drugs']} agents, {stats['mechanism_classes']} mechanism classes")
    print(f"  {stats['n_genes_targeted']}/{stats['n_signature_genes']} signature genes targeted "
          f"({stats['signature_gene_coverage']})")
    print(f"  evidence tiers: {stats['evidence_tiers']}")


if __name__ == "__main__":
    main()

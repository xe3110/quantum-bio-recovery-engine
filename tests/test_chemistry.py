"""Behavioural tests for the chemistry model.

The design engine is only as trustworthy as its ability to say what a
molecule *is*.  These tests pin the three things everything else rests on:
that SMILES survive a parse-write round trip, that formula and mass agree
with the literature for real drugs, and that invalid chemistry is refused
rather than silently scored.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from core.chemistry.descriptors import compute_descriptors, tpsa
from core.chemistry.druglikeness import (
    acidic_centres, basic_pka_estimate, cns_mpo, lipinski, profile_molecule,
    structural_alerts,
)
from core.chemistry.fingerprint import assess_novelty, similarity, tanimoto, circular_fingerprint
from core.chemistry.molecule import (
    SmilesError, attach, cap_attachments, parse_smiles, write_smiles,
)

ROOT = Path(__file__).resolve().parents[1]


def _curated_structures():
    """Every reference set any registered disease declares, not just MS's.

    Hardcoding one disease's file meant a second disease could add a novelty
    reference set that nothing validated -- and an unvalidated reference set is
    worse than none, because every novelty claim measured against it inherits
    the error silently.
    """
    from core.models.disease import available_diseases, load_disease

    entries = []
    for identifier in available_diseases():
        disease = load_disease(identifier)
        if disease.structures_path is None:
            continue
        payload = json.loads(Path(disease.structures_path).read_text())
        for entry in payload["structures"]:
            if entry.get("smiles"):
                entries.append({**entry, "disease": identifier})
    return entries


WITH_STRUCTURE = _curated_structures()

ROUND_TRIP_CORPUS = [
    "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "C1=CC=CN=C1", "c1cc[nH]c1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "C[N+](C)(C)CCO", "O=S(=O)(N)c1ccccc1",
    "C1CC2CCC1CC2", "c1ccc2ccccc2c1", "N#CC(=O)N", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "c1ccc(cc1)-c1ccccc1", "CN1CCN(CC1)c1ccc(cc1)C(=O)N", "c1ccsc1", "c1ccoc1",
]


@pytest.mark.parametrize("smiles", ROUND_TRIP_CORPUS)
def test_smiles_round_trips_without_changing_the_molecule(smiles):
    """Writing then re-reading must preserve composition exactly.

    A writer that drops a ring-closure digit produces a different molecule
    that still parses, which would corrupt every downstream descriptor
    silently. Formula and mass are the check that catches it.
    """
    original = parse_smiles(smiles)
    reparsed = parse_smiles(write_smiles(original))
    assert reparsed.molecular_formula() == original.molecular_formula()
    assert reparsed.molecular_weight() == pytest.approx(original.molecular_weight())


@pytest.mark.parametrize(
    "entry", WITH_STRUCTURE,
    ids=[f"{e['disease']}-{e['name']}" for e in WITH_STRUCTURE],
)
def test_curated_structures_match_their_literature_formula_and_mass(entry):
    """Every curated structure must reproduce its published formula and mass.

    This is what keeps a mistyped SMILES out of the novelty reference set: a
    wrong structure would make every novelty claim measured against it wrong
    too, without any other symptom.
    """
    descriptors = compute_descriptors(entry["smiles"])
    assert descriptors.formula == entry["formula"], entry["name"]
    assert descriptors.molecular_weight == pytest.approx(entry["mass"], abs=0.15)


@pytest.mark.parametrize("smiles,expected", [
    ("CCO", 20.23), ("c1ccccc1", 0.0), ("c1ccncc1", 12.89), ("Nc1ccccc1", 26.02),
    ("CC(=O)Oc1ccccc1C(=O)O", 63.60), ("CC(=O)Nc1ccc(O)cc1", 49.33),
    ("c1cc[nH]c1", 15.79), ("c1ccoc1", 13.14),
])
def test_tpsa_matches_published_ertl_values(smiles, expected):
    assert tpsa(parse_smiles(smiles)) == pytest.approx(expected, abs=0.02)


def test_aromatic_atom_with_exocyclic_double_bond_is_accepted():
    """Caffeine's ring carbonyls must not read as pentavalent carbon.

    An aromatic atom is normally credited a delocalised pi bond; an atom
    already carrying an explicit double bond is not. Getting this wrong
    rejects every purinone, xanthine, and pyridone in chemical space.
    """
    caffeine = parse_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    assert caffeine.molecular_formula() == "C8H10N4O2"


@pytest.mark.parametrize("smiles", ["C(C)(C)(C)(C)C", "c1ccccc1(C)(C)", "CC(=O)(=O)C"])
def test_over_valent_structures_are_refused(smiles):
    with pytest.raises(SmilesError):
        parse_smiles(smiles)


@pytest.mark.parametrize("smiles", ["c1ccccc", "CC(", "CC)", "C1CC", "[C@@H", ""])
def test_malformed_smiles_are_refused(smiles):
    with pytest.raises(SmilesError):
        parse_smiles(smiles)


def test_attachment_refuses_a_heteroatom_junction():
    """Joining two heteroatoms would build a peroxide or a hydrazine.

    Combinatorial assembly reaches these constantly -- an ether linker onto an
    alkoxy-terminated arm is a peroxide -- so the refusal has to live in the
    join itself, not in a downstream filter.
    """
    ether = parse_smiles("*O*")
    alkoxy_arm = parse_smiles("*OC(c1ccccc1)c1ccccc1")
    with pytest.raises(SmilesError, match="junction"):
        attach(ether, alkoxy_arm)


def test_attachment_joins_carbon_to_carbon_and_preserves_both_parts():
    core = parse_smiles("*c1ccc(*)cc1")
    arm = parse_smiles("*C(=O)N")
    joined = attach(core, arm)
    assert len(joined.attachment_points) == 1
    capped = cap_attachments(joined)
    assert not capped.attachment_points
    assert capped.molecular_formula() == "C7H7NO"


def test_capping_removes_every_dangling_attachment():
    molecule = cap_attachments(parse_smiles("*c1cc(*)cc(*)c1"))
    assert not molecule.attachment_points
    assert molecule.molecular_formula() == "C6H6"


def test_lipinski_flags_a_known_rule_of_five_violator():
    """Sirolimus-scale molecules must register as violating, not squeak through."""
    descriptors = compute_descriptors("CCCCCCCCCCCCCCCCCCCCC(=O)NCCCCCCCCCCCCCCCCO")
    assert lipinski(descriptors)["violations"] >= 1


def test_cns_mpo_prefers_a_cns_drug_over_a_large_polar_molecule():
    """Memantine is a CNS drug; methotrexate is not. The score must agree."""
    memantine = profile_molecule("CC12CC3CC(C)(C1)CC(N)(C3)C2")
    methotrexate = profile_molecule(
        "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)NC(CCC(=O)O)C(=O)O"
    )
    assert memantine.cns_mpo["total"] > methotrexate.cns_mpo["total"]


def test_cns_mpo_alone_does_not_catch_an_acid_but_the_acid_check_does():
    """A documented blind spot in Wager's score, covered by a separate check.

    MPO takes the most *basic* pKa as its ionisation term, so a small
    carboxylic acid scores well while being anionic at pH 7.4 and barred from
    the brain. The delivery gate has to consult both.
    """
    statin_acid = profile_molecule("CC(O)CC(O)CC(=O)O")
    assert statin_acid.cns_mpo["total"] > 4.0
    assert statin_acid.anionic_at_physiological_ph
    assert "carboxylic_acid" in statin_acid.acidic_centres


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=O)O", "carboxylic_acid"),
    ("OS(=O)(=O)c1ccccc1", "sulfonic_acid"),
    ("CC(=O)Nc1ccc(O)cc1", None),
    ("CCO", None),
])
def test_acidic_centre_detection(smiles, expected):
    """A phenol and an alcohol are not ionised at pH 7.4; an acid is."""
    found = acidic_centres(smiles)
    assert (expected in found) if expected else (found == [])


def test_basic_pka_finds_the_strongest_centre_and_names_it():
    pka, source = basic_pka_estimate("CCCCCCCCc1ccc(CCC(N)(CO)CO)cc1")
    assert source == "primary_aliphatic_amine"
    assert pka > 9.0


def test_amide_nitrogen_is_not_treated_as_basic():
    """An amide NH is not a basic centre; scoring it as one wrecks the MPO pKa term."""
    pka, source = basic_pka_estimate("CC(=O)Nc1ccc(O)cc1")
    assert source == "none"
    assert pka == 0.0


def test_fumarate_esters_are_flagged_as_michael_acceptors():
    """Dimethyl fumarate's electrophile is its mechanism and must be visible."""
    alerts = [a["alert"] for a in structural_alerts("COC(=O)C=CC(=O)OC")]
    assert "michael_acceptor" in alerts


def test_benign_molecule_carries_no_alerts():
    assert structural_alerts("CC(=O)Nc1ccc(O)cc1") == []


def test_fingerprint_is_stable_across_processes():
    """Novelty claims must be reproducible, so hashing cannot be salted."""
    assert circular_fingerprint("CC(=O)Oc1ccccc1C(=O)O") == circular_fingerprint(
        "CC(=O)Oc1ccccc1C(=O)O"
    )


def test_identical_molecules_are_maximally_similar_and_distinct_ones_are_not():
    assert similarity("CCO", "CCO") == 1.0
    assert similarity("CCO", "c1ccc2ccccc2c1") < 0.1


def test_novelty_recognises_a_compound_already_in_the_reference_set():
    library = [(e["name"], e["smiles"]) for e in WITH_STRUCTURE]
    report = assess_novelty("CCCCCCCCc1ccc(CCC(N)(CO)CO)cc1", library)
    assert report.nearest_name == "Fingolimod"
    assert report.max_similarity == 1.0
    assert not report.is_novel


def test_tanimoto_of_empty_fingerprints_is_defined():
    assert tanimoto(set(), set()) == 1.0
    assert tanimoto({1, 2}, set()) == 0.0

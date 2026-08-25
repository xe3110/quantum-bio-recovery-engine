"""The disease-agnostic claim, tested against both registered diseases.

The design stack previously claimed to be disease-agnostic on the strength of
one registry entry, while its own modules imported MS-specific scoring. These
tests hold the claim to something checkable: the same code path must derive a
sensible profile, Hamiltonian, and molecule for two diseases whose pathways,
therapeutic axes, and safety vocabularies have nothing in common.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from core.design.denovo import DesignWeights, run_design
from core.design.evidence import ASSESSMENT_AXES, evidence_status
from core.design.pharmacophores import load_library, unreachable_requirements
from core.design.quantum_assembly import build_selection_problem, solve_enumeration
from core.design.target_profile import build_target_profile
from core.models.disease import available_diseases, load_disease

ROOT = Path(__file__).resolve().parents[1]
DISEASES = ["multiple_sclerosis", "parkinsons"]


@pytest.fixture(scope="module")
def library():
    return load_library()


@pytest.mark.parametrize("identifier", DISEASES)
def test_registry_entry_loads_all_of_its_data(identifier):
    disease = load_disease(identifier)
    assert disease.signature().genes
    assert disease.network().number_of_nodes() > 0
    assert disease.panel()
    assert disease.druggability()


def test_registry_holds_more_than_one_disease():
    """The generality claim needs at least two instances to be checkable."""
    assert set(DISEASES) <= set(available_diseases())


def test_the_two_diseases_have_largely_different_safety_vocabularies():
    """Safety burden is a property of the disease, not of the scoring code.

    MS therapy is constrained by infection, malignancy, and autoimmunity;
    Parkinson's by dyskinesia, impulse-control disorders, and psychosis. The
    two overlap only on generic organ toxicity -- cardiac and hepatic -- which
    genuinely applies to both. What the registry has to carry is the majority
    that does not transfer.
    """
    ms = set(load_disease("multiple_sclerosis").risk_domains)
    pd = set(load_disease("parkinsons").risk_domains)
    assert ms - pd and pd - ms
    shared = ms & pd
    assert shared <= {"cardiac", "hepatic"}, f"unexpected shared domains: {shared}"
    assert len(shared) < min(len(ms), len(pd)) / 2


def test_the_two_diseases_share_no_therapeutic_axes():
    ms = set(load_disease("multiple_sclerosis").therapeutic_axes)
    pd = set(load_disease("parkinsons").therapeutic_axes)
    assert not (ms & pd)


@pytest.mark.parametrize("identifier", DISEASES)
def test_risk_weights_come_from_the_registry(identifier):
    disease = load_disease(identifier)
    assert set(disease.risk_weights) == set(disease.risk_domains)
    assert all(0.0 <= v <= 1.0 for v in disease.risk_weights.values())


def test_design_modules_do_not_import_disease_specific_scoring():
    """The import graph is part of the claim.

    A layer that says it is disease-agnostic while importing ``ms_scoring`` is
    not, whatever its registry contains.
    """
    for path in (ROOT / "core/design").glob("*.py"):
        assert "ms_scoring" not in path.read_text(), path
    assert "ms_scoring" not in (ROOT / "core/models/disease.py").read_text()


@pytest.mark.parametrize("identifier", DISEASES)
def test_profile_derives_and_ranks_for_each_disease(identifier, library):
    profile = build_target_profile(load_disease(identifier), top_n=12)
    assert profile.requirements
    priorities = [r.priority for r in profile.requirements]
    assert priorities == sorted(priorities, reverse=True)
    assert profile.property_window.rationale


@pytest.mark.parametrize("identifier", DISEASES)
def test_each_disease_has_a_fully_served_and_an_unserved_axis(identifier):
    """Both diseases have effective symptomatic therapy and no disease modification.

    The gap analysis should find that shape without being told it.
    """
    profile = build_target_profile(load_disease(identifier), top_n=12)
    gaps = profile.axis_gaps
    assert min(gaps.values()) == pytest.approx(0.0)
    assert max(gaps.values()) == pytest.approx(1.0)


def test_parkinsons_routes_around_its_least_tractable_central_target():
    """Alpha-synuclein is the central protein in PD and is not small-molecule tractable.

    The profile must not put a designed arm on it just because its leverage is
    high -- the same guard that keeps MS designs off CD20.
    """
    disease = load_disease("parkinsons")
    assert disease.druggability()["SNCA"]["small_molecule_tractability"] < 0.3
    profile = build_target_profile(disease, top_n=8)
    assert "SNCA" not in profile.genes


def test_gene_aliases_recover_a_target_the_interactome_names_differently():
    """STRING still calls glucocerebrosidase GBA; HGNC calls it GBA1.

    Unmapped, the most common genetic risk factor in PD would have zero
    network leverage and drop out of the profile silently.
    """
    disease = load_disease("parkinsons")
    assert disease.gene_aliases.get("GBA") == "GBA1"
    assert "GBA1" in disease.network()
    assert "GBA1" in build_target_profile(disease, top_n=14).genes


@pytest.mark.parametrize("identifier", DISEASES)
def test_hamiltonian_builds_and_solves_for_each_disease(identifier, library):
    profile = build_target_profile(load_disease(identifier), top_n=12)
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    best = solve_enumeration(problem, top=1)[0]
    assert best.feasible
    assert len(best.identifiers) == 2


@pytest.mark.parametrize("identifier", DISEASES)
def test_a_molecule_is_produced_for_each_disease(identifier, library):
    """The end-to-end claim: registry entry in, novel structure out."""
    from core.chemistry.molecule import parse_smiles

    disease = load_disease(identifier)
    profile = build_target_profile(disease, top_n=12)
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    arms = [s.identifiers for s in solve_enumeration(problem, top=1)]
    result = run_design(
        disease=disease, profile=profile, library=library, arm_sets=arms,
        weights=DesignWeights(), seed=7, top=3, iterations=20,
    )
    assert result["candidates"]
    for candidate in result["candidates"]:
        parse_smiles(candidate["smiles"])
        assert candidate["molecular_formula"]


@pytest.mark.parametrize("identifier", DISEASES)
def test_unreachable_requirements_are_reported_not_hidden(identifier, library):
    """A profile target with no chemical matter is the most useful thing a run says."""
    profile = build_target_profile(load_disease(identifier), top_n=14)
    unreachable = unreachable_requirements(library, profile.genes)
    assert isinstance(unreachable, list)
    assert not set(unreachable) - set(profile.genes)


# ---------------------------------------------------------------------------
# Evidence and provenance
# ---------------------------------------------------------------------------

def test_every_pharmacophore_declares_its_evidence_and_primary_targets(library):
    for fragment in library.pharmacophores:
        assert fragment.evidence_tier in (
            "approved_drug", "clinical_candidate", "published_chemotype", "speculative"
        )
        assert fragment.primary_targets
        assert not set(fragment.primary_targets) - set(fragment.engages)


def test_a_binding_claim_outranks_a_downstream_claim_from_the_same_fragment(library):
    """Otherwise a scaffold of speculative inference scores like real pharmacology."""
    fragment = library.get("nrf2_fumarate")
    primary = fragment.claim_confidence("NFE2L2")
    downstream = fragment.claim_confidence("NQO1")
    assert primary > downstream > 0


def test_approved_chemotypes_outrank_speculative_ones(library):
    assert (
        library.get("maob_propargylamine").mean_confidence
        > library.get("sirt_polyphenol").mean_confidence
    )


def test_confidence_weighting_reduces_coverage_for_weak_evidence(library):
    """The optimiser must have a reason to prefer well-evidenced arms."""
    from core.design.quantum_assembly import combined_confidence, profile_coverage

    profile = build_target_profile(load_disease("multiple_sclerosis"), top_n=14)
    fragment = library.get("gpr17_indole_acid")
    unweighted = profile_coverage(fragment.engages, profile)
    weighted = profile_coverage(fragment.engages, profile, combined_confidence([fragment]))
    assert 0 < weighted < unweighted


def test_designs_carry_a_machine_readable_evidence_status(library):
    """Caveats in a protocol do not travel with the data; this field does."""
    status = evidence_status().as_dict()
    assert status["readiness"] == "hypothesis_only"
    assert set(status["unassessed"]) == set(ASSESSMENT_AXES)
    for axis in ("target_binding", "efflux_liability", "cardiac_safety", "synthetic_route"):
        assert axis in status["detail"]


def test_provenance_records_input_digests_and_a_dirty_tree_flag():
    """A revision hash from a modified tree describes code that never existed."""
    from core.provenance import run_provenance

    payload = run_provenance(
        [ROOT / "data/ms_expression_v3.csv"], command=["pytest"]
    )
    assert payload["artifact_class"] == "run_snapshot"
    assert payload["inputs"][0]["sha256"]
    assert "dirty" in payload["git"]
    assert payload["environment"]["chemistry_backend"]["active"] in ("rdkit", "local")


def test_input_digest_changes_when_a_file_changes(tmp_path):
    from core.provenance import file_digest

    path = tmp_path / "signature.csv"
    path.write_text("gene,logFC,desired_direction\nA,1.0,-1\n")
    before = file_digest(path)["sha256"]
    path.write_text("gene,logFC,desired_direction\nA,1.1,-1\n")
    assert file_digest(path)["sha256"] != before

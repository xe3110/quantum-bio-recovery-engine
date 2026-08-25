"""Behavioural tests for the de novo design stack.

These cover the claims the pipeline makes about itself: that the target
profile refuses to point a small molecule at an antibody-only target, that the
Hamiltonian's truncated objective stays honest about what it approximates,
that assembly produces chemically valid structures, and that a seeded campaign
is reproducible.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import pytest

from core.chemistry.molecule import SmilesError
from core.design.denovo import (
    DesignWeights, Recipe, assemble, run_design, window_violations,
)
from core.design.pharmacophores import load_library, unreachable_requirements
from core.design.quantum_assembly import (
    SelectionWeights, arm_heavy_atom_budget, build_selection_problem,
    combined_effects, profile_coverage, redundancy, solve_enumeration,
)
from core.design.target_profile import build_target_profile
from core.models.disease import available_diseases, load_disease

ROOT = Path(__file__).resolve().parents[1]
DISEASE = "multiple_sclerosis"


@pytest.fixture(scope="module")
def disease():
    return load_disease(DISEASE)


@pytest.fixture(scope="module")
def library():
    return load_library()


@pytest.fixture(scope="module")
def profile(disease):
    return build_target_profile(disease, top_n=14)


# ---------------------------------------------------------------------------
# Disease model
# ---------------------------------------------------------------------------

def test_registry_lists_the_disease_and_loads_its_data(disease):
    assert DISEASE in available_diseases()
    assert len(disease.signature().genes) == 112
    assert disease.network().number_of_nodes() > 0
    assert len(disease.panel()) == 74


def test_every_panel_target_carries_a_druggability_annotation(disease):
    """An unannotated target would be silently dropped from the profile."""
    annotated = set(disease.druggability())
    panel_targets = {g for drug in disease.panel() for g in drug["target_effects"]}
    assert not panel_targets - annotated


def test_missing_disease_names_what_is_available():
    with pytest.raises(FileNotFoundError, match="Available"):
        load_disease("a_disease_that_does_not_exist")


def test_known_structures_exclude_biologics(disease):
    """Antibodies have no SMILES and must not enter the novelty reference set."""
    names = {name for name, _ in disease.known_structures()}
    assert "Ocrelizumab" not in names
    assert "Fingolimod" in names


# ---------------------------------------------------------------------------
# Target product profile
# ---------------------------------------------------------------------------

def test_profile_never_targets_an_antibody_only_antigen(profile):
    """CD20 has enormous leverage in MS and no small-molecule route.

    This is the single most important guard in the profile: without the
    tractability floor the engine would confidently specify a molecule to bind
    the target of ocrelizumab.
    """
    assert "MS4A1" not in profile.genes
    assert "CD52" not in profile.genes
    assert "MS4A1" in profile.readouts


def test_profile_requirements_respect_the_tractability_floor(disease):
    profile = build_target_profile(disease, top_n=20, tractability_floor=0.5)
    druggability = disease.druggability()
    for requirement in profile.requirements:
        assert druggability[requirement.gene]["small_molecule_tractability"] >= 0.5


def test_profile_directions_follow_the_signature_where_it_speaks(disease, profile):
    signature = disease.signature()
    for requirement in profile.requirements:
        if requirement.gene in signature.desired:
            assert requirement.desired_direction == signature.desired[requirement.gene]


def test_requirements_are_ranked_by_priority(profile):
    priorities = [r.priority for r in profile.requirements]
    assert priorities == sorted(priorities, reverse=True)


def test_axis_gaps_identify_remyelination_as_unserved(profile):
    """No approved MS therapy targets remyelination; the gap analysis must say so."""
    assert profile.axis_gaps["remyelination"] == pytest.approx(1.0)
    assert profile.axis_gaps["immunomodulation"] == pytest.approx(0.0)


def test_cns_disease_gets_a_tighter_envelope_than_a_peripheral_one(profile):
    window = profile.property_window
    assert window.molecular_weight[1] <= 420
    assert window.cns_mpo_floor >= 4.0
    assert window.rationale


# ---------------------------------------------------------------------------
# Fragment library
# ---------------------------------------------------------------------------

def test_every_fragment_parses_with_the_attachment_points_its_role_requires(library):
    """load_library validates on read, so loading at all is the assertion."""
    assert len(library) >= 40
    for fragment in library.pharmacophores:
        assert fragment.attachment_count == 1
    for linker in library.linkers:
        assert linker.attachment_count == 2
    for scaffold in library.scaffolds:
        assert scaffold.attachment_count >= 2


def test_every_profile_requirement_has_chemical_matter(library, profile):
    assert unreachable_requirements(library, profile.genes) == []


def test_engagement_values_stay_within_the_signed_convention(library):
    for fragment in library.fragments:
        for value in fragment.engages.values():
            assert -1.0 <= value <= 1.0


# ---------------------------------------------------------------------------
# The Hamiltonian
# ---------------------------------------------------------------------------

def test_every_pair_has_a_stored_coupling(profile, library):
    """Couplings are keyed on the sorted pair, not on ranking order.

    Storing them in ranked order and reading them back sorted made every
    coupling whose ranking disagreed with the alphabet return zero, so the
    optimiser silently solved a different Hamiltonian than the one reported.
    """
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    for pair in combinations(problem.identifiers, 2):
        assert tuple(sorted(pair)) in problem.quadratic
        problem.qubo_objective(pair)  # must not raise


def test_size_linearisation_never_under_charges_a_design(profile, library):
    """The QUBO may reject a feasible design; it must never admit an infeasible one.

    Whole-molecule budgets cannot be written as a quadratic, so each fragment
    is charged against a 1/k share. Subadditivity makes that conservative,
    and this pins the direction of the error.
    """
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    for subset in combinations(problem.identifiers, 2):
        assert problem.qubo_objective(subset) <= problem.exact_objective(subset) + 1e-5


def test_coverage_component_is_exact_for_two_arms(profile, library):
    """Setting the budget costs to zero, the truncation is exact at k=2."""
    free = SelectionWeights(mass_weight=0.0, envelope_weight=0.0, liability_weight=0.0)
    problem = build_selection_problem(profile, library, k=2, max_variables=8, weights=free)
    for subset in combinations(problem.identifiers, 2):
        assert problem.qubo_objective(subset) == pytest.approx(
            problem.exact_objective(subset), abs=1e-5
        )


def test_enumeration_returns_the_true_optimum(profile, library):
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    best = solve_enumeration(problem, top=1)[0]
    every = [
        problem.qubo_objective(s) for s in combinations(problem.identifiers, 2)
    ]
    assert best.qubo_objective == pytest.approx(max(every))
    assert best.feasible
    assert len(best.identifiers) == 2


def test_arm_budget_reserves_mass_for_the_scaffold_and_linkers(profile):
    """The arms cannot be handed the whole molecular-weight ceiling.

    Doing so produced 550 Da molecules against a 420 Da envelope, because
    assembly adds a scaffold and a linker per arm afterwards.
    """
    ceiling_atoms = profile.property_window.molecular_weight[1] / 13.5
    assert arm_heavy_atom_budget(profile, 3) < ceiling_atoms


def test_redundancy_penalises_two_arms_from_one_chemotype(library):
    weights = SelectionWeights()
    same_family = redundancy(
        library.get("s1p_azetidine_acid"), library.get("s1p_aminodiol"), weights
    )
    different = redundancy(
        library.get("s1p_azetidine_acid"), library.get("nmda_adamantyl"), weights
    )
    assert same_family > different


def test_coverage_is_negative_when_effects_oppose_the_profile(profile, library):
    """Pushing a target the wrong way must cost, not merely fail to help."""
    fragment = library.get("nlrp3_sulfonylurea")
    inverted = {gene: -value for gene, value in fragment.engages.items()}
    assert profile_coverage(fragment.engages, profile) > 0
    assert profile_coverage(inverted, profile) < 0


def test_bliss_combination_saturates_rather_than_summing(library):
    a = library.get("csf1r_picolinamide")
    combined = combined_effects([a, a])
    for gene, value in combined.items():
        assert abs(value) <= 1.0
        assert abs(value) >= abs(a.engages[gene])


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def test_assembly_produces_a_capped_valid_molecule(library):
    recipe = Recipe(
        "scaffold_benzene_14",
        ("btk_acrylamide_piperidine", "csf1r_picolinamide"),
        ("linker_direct", "linker_amide"),
        (),
    )
    molecule = assemble(library, recipe)
    assert not molecule.attachment_points
    assert molecule.molecular_weight() > 200


def test_assembly_refuses_more_arms_than_the_scaffold_can_hold(library):
    recipe = Recipe(
        "scaffold_benzene_14",
        ("btk_acrylamide_piperidine", "csf1r_picolinamide", "nmda_adamantyl"),
        ("linker_direct",) * 3,
        (),
    )
    with pytest.raises(SmilesError, match="attachment points"):
        assemble(library, recipe)


def test_window_violations_name_the_property_and_the_bound(profile):
    from core.chemistry.descriptors import compute_descriptors

    violations = window_violations(
        compute_descriptors("CCCCCCCCCCCCCCCCCCCCCCCCCC"), profile.property_window
    )
    assert any("clogp" in v for v in violations)


# ---------------------------------------------------------------------------
# The campaign
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def campaign(disease, profile, library):
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    arms = [s.identifiers for s in solve_enumeration(problem, top=1)]
    return run_design(
        disease=disease, profile=profile, library=library, arm_sets=arms,
        weights=DesignWeights(), seed=7, top=5, iterations=30,
    )


def test_campaign_produces_valid_novel_structures(campaign):
    from core.chemistry.molecule import parse_smiles

    assert campaign["candidates"]
    for candidate in campaign["candidates"]:
        parse_smiles(candidate["smiles"])  # must round-trip
        assert candidate["molecular_formula"]
        assert candidate["novelty"]["max_tanimoto"] < 1.0


def test_campaign_is_reproducible_under_a_fixed_seed(disease, profile, library):
    problem = build_selection_problem(profile, library, k=2, max_variables=8)
    arms = [s.identifiers for s in solve_enumeration(problem, top=1)]
    kwargs = dict(
        disease=disease, profile=profile, library=library, arm_sets=arms,
        weights=DesignWeights(), top=5, iterations=30,
    )
    first = run_design(seed=11, **kwargs)
    second = run_design(seed=11, **kwargs)
    assert [c["smiles"] for c in first["candidates"]] == [
        c["smiles"] for c in second["candidates"]
    ]


def test_campaign_reports_rather_than_hides_assembly_failures(campaign):
    """Invalid recipes are expected; silently dropping them would misreport search effort."""
    search = campaign["search"]
    assert search["recipes_attempted"] > 0
    assert search["assembly_failures"] >= 0
    assert search["unique_structures"] <= search["recipes_attempted"]


def test_designs_carry_provenance_back_to_their_chemotypes(campaign):
    for candidate in campaign["candidates"]:
        assert candidate["provenance"]["arm_chemotypes"]
        assert candidate["recipe"]["arms"]


def test_acidic_designs_fail_the_cns_delivery_gate(disease, profile, library):
    """An anionic molecule cannot reach the brain, whatever its MPO says."""
    from core.chemistry.druglikeness import profile_molecule
    from core.design.denovo import score_candidate

    recipe = Recipe(
        "scaffold_benzene_14", ("hmgcr_dihydroxyacid",), ("linker_direct",), ("cap_methyl",)
    )
    molecule = assemble(library, recipe)
    candidate = score_candidate(
        molecule, recipe, library, profile, disease.signature(),
        DesignWeights(), require_cns=True,
    )
    assert profile_molecule(molecule).anionic_at_physiological_ph
    assert not candidate.meets_delivery_gate

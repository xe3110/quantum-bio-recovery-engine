"""De novo assembly: turning a fragment selection into a proposed molecule.

This is where the engine stops ranking existing drugs and starts producing
structures.  A design is a *recipe* -- a scaffold, an arm per attachment
point, a linker for each arm, and a cap for whatever is left over -- and the
recipe is assembled into a real molecular graph, capped, and measured.

The search runs in two stages.  Enumeration walks every scaffold, linker, and
cap combination for a given set of arms, which is exhaustive at this library
size.  Refinement then hill-climbs from the best recipes, mutating one
component at a time, which matters because the arm set is fixed by the
Hamiltonian while the assembly around it is not.

What a designed molecule here *is*: a structure whose predicted
polypharmacology follows from its parts, whose physicochemical properties are
computed from its graph, and whose novelty is measured against known
chemical matter.

What it is not: docked, simulated, synthesised, or assayed.  There is no
binding-affinity model anywhere in this pipeline.  A high-scoring molecule is
a hypothesis worth a synthesis and a panel of assays -- the first step of a
programme, not a result.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Any, Iterable, Sequence

import numpy as np

from core.biology.signature import Signature, alignment_metrics
from core.chemistry.descriptors import Descriptors, compute_descriptors
from core.chemistry.druglikeness import DevelopabilityProfile, profile_molecule
from core.chemistry.fingerprint import NoveltyReport, assess_novelty, diversity_select
from core.chemistry.molecule import Molecule, SmilesError, attach, cap_attachments, parse_smiles
from core.design.evidence import evidence_status
from core.design.pharmacophores import Fragment, FragmentLibrary
from core.design.quantum_assembly import combined_confidence, combined_effects, profile_coverage
from core.design.target_profile import PropertyWindow, TargetProductProfile
from core.models.disease import DiseaseContext


@dataclass(frozen=True)
class DesignWeights:
    """Declared weights for the design fitness function.

    Coverage of the target profile leads; everything else is a developability
    term that decides between molecules of comparable coverage.  A design that
    hits every target and cannot reach the tissue is not a design.
    """

    profile_coverage: float = 1.00
    signature_alignment: float = 0.55
    cns_exposure: float = 0.50
    drug_likeness: float = 0.35
    tractability: float = 0.30
    ligand_efficiency: float = 0.25
    counter_therapeutic_penalty: float = 0.60
    alert_penalty: float = 0.10
    window_penalty: float = 0.45
    delivery_gate_penalty: float = 0.50

    def as_dict(self) -> dict[str, float]:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


# A design achieving this much profile coverage per heavy atom is treated as
# fully efficient. Set from the library's own best solo fragments rather than
# from a literature ligand-efficiency value, which is defined on binding free
# energy and is not what this quantity measures.
LIGAND_EFFICIENCY_REFERENCE = 0.012


@dataclass(frozen=True)
class Recipe:
    """The assembly instructions for one candidate."""

    scaffold: str
    arms: tuple[str, ...]
    linkers: tuple[str, ...]
    caps: tuple[str, ...]

    def as_dict(self) -> dict[str, Any]:
        return {
            "scaffold": self.scaffold,
            "arms": list(self.arms),
            "linkers": list(self.linkers),
            "caps": list(self.caps),
        }

    def key(self) -> tuple:
        return (self.scaffold, self.arms, self.linkers, self.caps)


@dataclass(frozen=True)
class DesignCandidate:
    """A proposed molecule and everything computed about it."""

    smiles: str
    formula: str
    recipe: Recipe
    predicted_effects: dict[str, float]
    descriptors: Descriptors
    developability: DevelopabilityProfile
    profile_coverage: float
    signature_metrics: dict[str, float]
    ligand_efficiency: float
    window_violations: tuple[str, ...]
    meets_delivery_gate: bool
    fitness: float
    novelty: NoveltyReport | None = None
    provenance: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "smiles": self.smiles,
            "molecular_formula": self.formula,
            "fitness": self.fitness,
            "recipe": self.recipe.as_dict(),
            "profile_coverage": self.profile_coverage,
            "signature_metrics": self.signature_metrics,
            "ligand_efficiency": self.ligand_efficiency,
            "meets_delivery_gate": self.meets_delivery_gate,
            "acidic_centres": self.developability.acidic_centres,
            "window_violations": list(self.window_violations),
            "descriptors": self.descriptors.as_dict(),
            "cns_mpo": self.developability.cns_mpo,
            "lipinski": self.developability.lipinski,
            "veber": self.developability.veber,
            "drug_likeness": self.developability.drug_likeness,
            "synthetic_tractability": self.developability.tractability,
            "structural_alerts": self.developability.alerts,
            "predicted_target_effects": {
                k: round(v, 4) for k, v in sorted(self.predicted_effects.items())
            },
            "novelty": self.novelty.as_dict() if self.novelty else None,
            "evidence_status": evidence_status().as_dict(),
            "provenance": self.provenance,
        }


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def build_arm(library: FragmentLibrary, arm_id: str, linker_id: str) -> Molecule:
    """Join one pharmacophore to its linker, leaving a single free point.

    The linker's first attachment point takes the arm; its second survives to
    connect to the scaffold.  ``linker_direct`` is the degenerate case: its
    two attachment atoms are bonded to each other, so consuming one leaves the
    arm carrying a bare point, which is exactly a direct bond.
    """
    arm = library.get(arm_id).molecule()
    linker = library.get(linker_id).molecule()
    return attach(linker, arm)


def assemble(library: FragmentLibrary, recipe: Recipe) -> Molecule:
    """Assemble a recipe into a capped molecular graph.

    Raises :class:`SmilesError` if the recipe cannot be built -- too many arms
    for the scaffold, or a valence that cannot be satisfied.  Failures are
    expected during search and are caught by the caller.
    """
    scaffold_fragment = library.get(recipe.scaffold)
    molecule = scaffold_fragment.molecule()
    points = len(molecule.attachment_points)
    if len(recipe.arms) > points:
        raise SmilesError(
            f"Scaffold {recipe.scaffold!r} offers {points} attachment points but the recipe "
            f"supplies {len(recipe.arms)} arms"
        )
    for arm_id, linker_id in zip(recipe.arms, recipe.linkers):
        molecule = attach(molecule, build_arm(library, arm_id, linker_id))
    for cap_id in recipe.caps:
        if not molecule.attachment_points:
            break
        molecule = attach(molecule, library.get(cap_id).molecule())
    return cap_attachments(molecule)


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def window_violations(descriptors: Descriptors, window: PropertyWindow) -> tuple[str, ...]:
    """Which of the profile's property limits a molecule breaks, and by how much."""
    checks = (
        ("molecular_weight", descriptors.molecular_weight, window.molecular_weight),
        ("clogp", descriptors.clogp, window.clogp),
        ("tpsa", descriptors.tpsa, window.tpsa),
        ("hbd", float(descriptors.hbd), window.hbd),
        ("rotatable_bonds", float(descriptors.rotatable_bonds), window.rotatable_bonds),
    )
    breaks = []
    for name, value, (lower, upper) in checks:
        if value < lower:
            breaks.append(f"{name}={value:g} below {lower:g}")
        elif value > upper:
            breaks.append(f"{name}={value:g} above {upper:g}")
    return tuple(breaks)


def score_candidate(
    molecule: Molecule,
    recipe: Recipe,
    library: FragmentLibrary,
    profile: TargetProductProfile,
    signature: Signature,
    weights: DesignWeights,
    require_cns: bool,
) -> DesignCandidate:
    """Compute every quantity for one assembled molecule and its fitness."""
    arms = [library.get(a) for a in recipe.arms]
    effects = combined_effects(arms)
    confidence = combined_confidence(arms)
    coverage = profile_coverage(effects, profile, confidence)
    alignment = alignment_metrics(effects, signature)

    developability = profile_molecule(molecule)
    descriptors = developability.descriptors
    violations = window_violations(descriptors, profile.property_window)

    mpo = developability.cns_mpo
    # Two independent ways to fail delivery: too little of the multi-parameter
    # score, or an acid group that MPO cannot see but that leaves the molecule
    # anionic at pH 7.4 and effectively excluded from the brain.
    meets_gate = (not require_cns) or (
        mpo["total"] >= profile.property_window.cns_mpo_floor
        and not developability.anionic_at_physiological_ph
    )

    efficiency = coverage / max(descriptors.heavy_atoms, 1)
    efficiency_norm = float(np.clip(efficiency / LIGAND_EFFICIENCY_REFERENCE, 0.0, 1.0))

    fitness = (
        weights.profile_coverage * coverage
        + weights.signature_alignment * alignment["reversal_efficiency"]
        + weights.cns_exposure * mpo["normalised"]
        + weights.drug_likeness * developability.drug_likeness["score"]
        + weights.tractability * developability.tractability["score"]
        + weights.ligand_efficiency * efficiency_norm
        - weights.counter_therapeutic_penalty * alignment["counter_therapeutic"]
        - weights.alert_penalty * len(developability.alerts)
        - weights.window_penalty * (len(violations) / 5.0)
        - (0.0 if meets_gate else weights.delivery_gate_penalty)
    )

    return DesignCandidate(
        smiles=molecule.to_smiles(),
        formula=molecule.molecular_formula(),
        recipe=recipe,
        predicted_effects=effects,
        descriptors=descriptors,
        developability=developability,
        profile_coverage=round(coverage, 5),
        signature_metrics={k: round(float(v), 5) for k, v in alignment.items()},
        ligand_efficiency=round(efficiency, 6),
        window_violations=violations,
        meets_delivery_gate=meets_gate,
        fitness=round(float(fitness), 5),
        provenance={
            "arm_chemotypes": [f.parent_chemotype for f in arms],
            "arm_families": [f.target_family for f in arms],
            "arm_evidence_tiers": [f.evidence_tier for f in arms],
            "therapeutic_axes": sorted({ax for f in arms for ax in f.therapeutic_axes}),
            "engagement_confidence": {g: round(c, 4) for g, c in sorted(confidence.items())},
            "mean_engagement_confidence": round(
                sum(confidence.values()) / len(confidence), 4
            ) if confidence else 0.0,
            "weakest_claim": min(confidence.items(), key=lambda kv: kv[1])[0] if confidence else None,
        },
    )


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

@dataclass
class DesignCampaign:
    """A seeded search over assemblies for one fixed set of arms."""

    library: FragmentLibrary
    profile: TargetProductProfile
    signature: Signature
    require_cns: bool
    weights: DesignWeights = field(default_factory=DesignWeights)
    seed: int = 7
    max_enumeration: int = 20000

    def __post_init__(self) -> None:
        self._rng = random.Random(self.seed)
        self._cache: dict[tuple, DesignCandidate | None] = {}
        self.attempted = 0
        self.failed = 0

    # -- evaluation -----------------------------------------------------
    def evaluate(self, recipe: Recipe) -> DesignCandidate | None:
        key = recipe.key()
        if key in self._cache:
            return self._cache[key]
        self.attempted += 1
        try:
            molecule = assemble(self.library, recipe)
            candidate = score_candidate(
                molecule, recipe, self.library, self.profile,
                self.signature, self.weights, self.require_cns,
            )
        except (SmilesError, KeyError, ValueError):
            self.failed += 1
            candidate = None
        self._cache[key] = candidate
        return candidate

    # -- stage one: enumeration -----------------------------------------
    def enumerate_assemblies(self, arms: Sequence[str]) -> list[DesignCandidate]:
        """Every scaffold, linker, and cap combination for a fixed arm set."""
        scaffolds = [
            s for s in self.library.scaffolds
            if s.attachment_count >= len(arms)
        ]
        linkers = [l.identifier for l in self.library.linkers]
        caps = [c.identifier for c in self.library.caps]

        results: list[DesignCandidate] = []
        for scaffold in scaffolds:
            spare = scaffold.attachment_count - len(arms)
            cap_options = list(product(caps, repeat=spare)) if spare else [()]
            for linker_choice in product(linkers, repeat=len(arms)):
                for cap_choice in cap_options:
                    if self.attempted >= self.max_enumeration:
                        break
                    recipe = Recipe(scaffold.identifier, tuple(arms), linker_choice, cap_choice)
                    candidate = self.evaluate(recipe)
                    if candidate is not None:
                        results.append(candidate)
        results.sort(key=lambda c: (-c.fitness, c.smiles))
        return results

    # -- stage two: refinement ------------------------------------------
    def _mutate(self, recipe: Recipe, allowed_arms: Sequence[str]) -> Recipe:
        """One random single-component change to a recipe."""
        rng = self._rng
        choice = rng.choice(("scaffold", "linker", "cap", "arm"))
        scaffolds = [s.identifier for s in self.library.scaffolds
                     if s.attachment_count >= len(recipe.arms)]
        if choice == "scaffold" and scaffolds:
            new_scaffold = rng.choice(scaffolds)
            spare = self.library.get(new_scaffold).attachment_count - len(recipe.arms)
            caps = tuple(
                rng.choice([c.identifier for c in self.library.caps]) for _ in range(spare)
            )
            return Recipe(new_scaffold, recipe.arms, recipe.linkers, caps)
        if choice == "linker" and recipe.linkers:
            index = rng.randrange(len(recipe.linkers))
            linkers = list(recipe.linkers)
            linkers[index] = rng.choice([l.identifier for l in self.library.linkers])
            return Recipe(recipe.scaffold, recipe.arms, tuple(linkers), recipe.caps)
        if choice == "cap" and recipe.caps:
            index = rng.randrange(len(recipe.caps))
            caps = list(recipe.caps)
            caps[index] = rng.choice([c.identifier for c in self.library.caps])
            return Recipe(recipe.scaffold, recipe.arms, recipe.linkers, tuple(caps))
        if choice == "arm" and len(allowed_arms) > len(recipe.arms):
            index = rng.randrange(len(recipe.arms))
            arms = list(recipe.arms)
            replacements = [a for a in allowed_arms if a not in arms]
            if replacements:
                arms[index] = rng.choice(replacements)
                return Recipe(recipe.scaffold, tuple(arms), recipe.linkers, recipe.caps)
        return recipe

    def refine(
        self,
        seeds: Sequence[DesignCandidate],
        allowed_arms: Sequence[str],
        iterations: int = 400,
        restarts: int = 6,
    ) -> list[DesignCandidate]:
        """Hill-climb from the best enumerated recipes.

        Accepts only strict improvements, so the walk is a local search rather
        than annealing: enumeration has already covered the assembly space
        broadly, and this stage is for the arm substitutions and rare
        combinations that enumeration holds fixed.
        """
        found: dict[str, DesignCandidate] = {}
        for seed_candidate in seeds[:restarts]:
            current = seed_candidate
            found[current.smiles] = current
            for _ in range(iterations):
                proposal = self._mutate(current.recipe, allowed_arms)
                if proposal.key() == current.recipe.key():
                    continue
                candidate = self.evaluate(proposal)
                if candidate is None:
                    continue
                found.setdefault(candidate.smiles, candidate)
                if candidate.fitness > current.fitness:
                    current = candidate
        return sorted(found.values(), key=lambda c: (-c.fitness, c.smiles))


def run_design(
    disease: DiseaseContext,
    profile: TargetProductProfile,
    library: FragmentLibrary,
    arm_sets: Sequence[Sequence[str]],
    weights: DesignWeights | None = None,
    seed: int = 7,
    top: int = 10,
    iterations: int = 400,
    novelty_threshold: float = 0.6,
) -> dict[str, Any]:
    """Run assembly, refinement, novelty assessment, and diversity selection.

    ``arm_sets`` are the fragment selections handed over by the Hamiltonian --
    typically the optimum plus the next few, so the campaign explores more
    than one mechanistic hypothesis.
    """
    weights = weights or DesignWeights()
    signature = disease.signature()
    campaign = DesignCampaign(
        library=library,
        profile=profile,
        signature=signature,
        require_cns=disease.delivery.requires_cns_exposure,
        weights=weights,
        seed=seed,
    )

    everything: dict[str, DesignCandidate] = {}
    per_arm_set = []
    for arms in arm_sets:
        enumerated = campaign.enumerate_assemblies(list(arms))
        refined = campaign.refine(enumerated[:6], list(arms), iterations=iterations)
        for candidate in (*enumerated, *refined):
            existing = everything.get(candidate.smiles)
            if existing is None or candidate.fitness > existing.fitness:
                everything[candidate.smiles] = candidate
        per_arm_set.append({
            "arms": list(arms),
            "n_assemblies": len(enumerated),
            "best_fitness": enumerated[0].fitness if enumerated else None,
            "best_smiles": enumerated[0].smiles if enumerated else None,
        })

    ranked = sorted(everything.values(), key=lambda c: (-c.fitness, c.smiles))

    reference = [*disease.known_structures(), *library.parent_structures()]
    scored: list[DesignCandidate] = []
    for candidate in ranked[: max(top * 6, 60)]:
        novelty = assess_novelty(candidate.smiles, reference, threshold=novelty_threshold)
        scored.append(
            DesignCandidate(**{**candidate.__dict__, "novelty": novelty})
        )

    # Advancing near-identical molecules tests one hypothesis, not several.
    indices = diversity_select([(c.formula, c.smiles) for c in scored], k=top)
    selected = [scored[i] for i in indices]

    return {
        "candidates": [c.as_dict() for c in selected],
        "all_ranked": [c.as_dict() for c in scored[:top]],
        "per_arm_set": per_arm_set,
        "search": {
            "recipes_attempted": campaign.attempted,
            "assembly_failures": campaign.failed,
            "unique_structures": len(everything),
            "seed": seed,
            "refinement_iterations": iterations,
            "weights": weights.as_dict(),
            "novelty_threshold": novelty_threshold,
            "ligand_efficiency_reference": LIGAND_EFFICIENCY_REFERENCE,
        },
    }

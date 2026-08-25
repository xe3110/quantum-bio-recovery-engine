"""Fragment selection as a Hamiltonian, solved exactly and on QAOA.

Choosing which pharmacophores to fuse into one molecule is a constrained
binary optimisation, and not a trivial one.  Each fragment carries a benefit
of its own; every pair carries an interaction -- Bliss saturation when two
arms push the same target, redundancy when they come from the same chemotype,
and a shared mass budget that a CNS-restricted profile cannot overrun.  That
is exactly a QUBO:

    maximise    sum_f b_f x_f  +  sum_{f<g} c_fg x_f x_g
    subject to  sum_f x_f = k,        x_f in {0,1}

The pairwise form is a truncation in two places.  ``b_f`` is a fragment's
solo coverage of the target profile minus its costs, and ``c_fg`` is the
*interaction* term ``cov(f,g) - cov(f) - cov(g)`` minus redundancy.  The
coverage part is exact at ``k = 2`` and a second-order approximation beyond
it, slightly over-crediting three-way combinations of a saturating function.

The size and envelope costs are a second, separate approximation.  Mass,
polar surface, donors, and lipophilicity are budgets on the *whole* molecule,
which a quadratic form cannot express, so each fragment is charged against a
``1/k`` share instead.  That linearisation is deliberately conservative:
because ``sum_f max(0, x_f - s) >= max(0, sum_f x_f - k*s)``, the QUBO never
under-charges a design, and a molecule the Hamiltonian rejects on size may
still have fitted.

:func:`exact_objective` recomputes the true value -- joint budgets and all --
for any selection, and every solver's result is reported against it.

Three backends solve the same Hamiltonian:

``enumeration``
    Brute force over all C(n,k) feasible subsets.  Cheap at this size and the
    ground truth everything else is checked against.
``eigensolver``
    Qiskit's ``MinimumEigenOptimizer`` with an exact classical eigensolver on
    the Ising Hamiltonian.  This is the standard quantum-algorithms baseline:
    the same operator a device would receive, diagonalised classically.
``qaoa``
    QAOA on the Aer simulator, seeded, swept across circuit depths because a
    single depth is not reliable -- see :func:`solve_qaoa_sweep`.

Running QAOA against exact enumeration on a problem small enough to solve
both ways is the honest use of quantum optimisation here.  It demonstrates
the formulation transfers to a device; it does not demonstrate an advantage,
and at this problem size there is none to demonstrate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Any, Iterable, Sequence

import numpy as np

from core.biology.signature import bliss_combine
from core.design.pharmacophores import Fragment, FragmentLibrary
from core.design.target_profile import TargetProductProfile

try:  # pragma: no cover - environment dependent
    from qiskit_optimization import QuadraticProgram
    from qiskit_optimization.algorithms import MinimumEigenOptimizer

    HAS_QISKIT = True
except Exception:  # pragma: no cover
    QuadraticProgram = None
    MinimumEigenOptimizer = None
    HAS_QISKIT = False


@dataclass(frozen=True)
class SelectionWeights:
    """Declared costs in the fragment-selection objective.

    Calibrated so profile coverage drives the choice and the costs modulate
    it.  An earlier calibration charged mass linearly from the first atom,
    which made a four-atom zinc binder beat every real pharmacophore in the
    library -- the optimiser was minimising molecular weight, not designing a
    drug.  Size is a budget, so it is charged only where a fragment exceeds
    its share of one.
    """

    mass_weight: float = 0.35
    liability_weight: float = 0.03
    redundancy_weight: float = 0.20
    family_clash_penalty: float = 0.30
    envelope_weight: float = 0.60

    def as_dict(self) -> dict[str, float]:
        return {
            "mass_weight": self.mass_weight,
            "liability_weight": self.liability_weight,
            "redundancy_weight": self.redundancy_weight,
            "family_clash_penalty": self.family_clash_penalty,
            "envelope_weight": self.envelope_weight,
        }


def profile_coverage(
    effects: dict[str, float],
    profile: TargetProductProfile,
    confidence: dict[str, float] | None = None,
) -> float:
    """Priority-weighted, direction-aware coverage of the target profile.

    Positive where the effects push each requirement the way the profile asks,
    negative where they push against it. Normalised by total profile priority
    so the value sits on [-1, 1] and is comparable across profiles.

    ``confidence`` discounts each contribution by how well-evidenced the claim
    behind it is. Without it, an arm built on speculative downstream inference
    buys exactly as much coverage as one built on approved-drug pharmacology,
    and the optimiser has no reason to prefer the evidence.
    """
    total = sum(r.priority for r in profile.requirements)
    if total <= 0:
        return 0.0
    score = 0.0
    for requirement in profile.requirements:
        effect = effects.get(requirement.gene)
        if not effect:
            continue
        aligned = float(np.clip(requirement.desired_direction * effect, -1.0, 1.0))
        weight = 1.0 if confidence is None else confidence.get(requirement.gene, 0.0)
        score += requirement.priority * aligned * weight
    return score / total


def combined_confidence(fragments: Sequence[Fragment]) -> dict[str, float]:
    """Per-gene confidence for a set of arms: the best claim behind each gene.

    Taking the maximum rather than a mean means a well-evidenced arm is not
    penalised for sharing a molecule with a speculative one, while a gene
    supported only by speculation stays discounted.
    """
    best: dict[str, float] = {}
    for fragment in fragments:
        for gene in fragment.engages:
            value = fragment.claim_confidence(gene)
            if value > best.get(gene, 0.0):
                best[gene] = value
    return best


def combined_effects(fragments: Sequence[Fragment]) -> dict[str, float]:
    """Union of engagement vectors, combined under Bliss independence."""
    combined: dict[str, float] = {}
    for fragment in fragments:
        for gene, value in fragment.engages.items():
            combined[gene] = bliss_combine(combined.get(gene, 0.0), value)
    return combined


def redundancy(a: Fragment, b: Fragment, weights: SelectionWeights) -> float:
    """How much two fragments duplicate each other, in [0, ~1].

    Target-set overlap plus an explicit penalty for sharing a chemotype
    family.  Two arms from the same family are one hypothesis wearing two
    hats: they cost mass twice and buy coverage once.
    """
    union = a.targets | b.targets
    overlap = len(a.targets & b.targets) / len(union) if union else 0.0
    clash = weights.family_clash_penalty if (
        a.target_family and a.target_family == b.target_family
    ) else 0.0
    return overlap + clash


@dataclass(frozen=True)
class FragmentSelectionProblem:
    """A QUBO over fragment choice, with everything needed to audit it."""

    fragments: tuple[Fragment, ...]
    linear: dict[str, float]
    quadratic: dict[tuple[str, str], float]
    k: int
    profile: TargetProductProfile
    weights: SelectionWeights
    heavy_atom_budget: int

    @property
    def identifiers(self) -> list[str]:
        return [f.identifier for f in self.fragments]

    def fragment(self, identifier: str) -> Fragment:
        for fragment in self.fragments:
            if fragment.identifier == identifier:
                return fragment
        raise KeyError(identifier)

    def qubo_objective(self, selection: Iterable[str]) -> float:
        """The truncated pairwise objective the solvers actually optimise."""
        chosen = list(selection)
        value = sum(self.linear[i] for i in chosen)
        for pair in combinations(sorted(chosen), 2):
            if pair not in self.quadratic:
                raise KeyError(f"No coupling stored for the pair {pair}")
            value += self.quadratic[pair]
        return round(value, 6)

    def exact_objective(self, selection: Iterable[str]) -> float:
        """The true objective, recomputed without the pairwise truncation."""
        fragments = [self.fragment(i) for i in selection]
        if not fragments:
            return 0.0
        coverage = profile_coverage(
            combined_effects(fragments), self.profile, combined_confidence(fragments)
        )
        used = sum(f.heavy_atoms for f in fragments)
        mass = max(0.0, used - self.heavy_atom_budget) / self.heavy_atom_budget
        window = self.profile.property_window
        tpsa_budget = max(window.tpsa[1], 1.0)
        hbd_budget = max(window.hbd[1], 1.0)
        logp_budget = max(window.clogp[1], 1.0)
        descriptors = [f.capped_descriptors() for f in fragments]
        envelope = (
            max(0.0, sum(d.tpsa for d in descriptors) - tpsa_budget) / tpsa_budget
            + max(0.0, sum(d.hbd for d in descriptors) - hbd_budget) / hbd_budget
            + max(0.0, sum(d.clogp for d in descriptors) - logp_budget) / logp_budget
        )
        liabilities = sum(len(f.liabilities) for f in fragments)
        penalty = 0.0
        for a, b in combinations(fragments, 2):
            penalty += self.weights.redundancy_weight * redundancy(a, b, self.weights)
        return round(
            coverage
            - self.weights.mass_weight * mass
            - self.weights.envelope_weight * envelope
            - self.weights.liability_weight * liabilities
            - penalty,
            6,
        )

    def hamiltonian_summary(self) -> dict[str, Any]:
        return {
            "n_variables": len(self.fragments),
            "n_linear_terms": len(self.linear),
            "n_quadratic_terms": len(self.quadratic),
            "cardinality_constraint": self.k,
            "feasible_states": _binomial(len(self.fragments), self.k),
            "search_space": 2 ** len(self.fragments),
        }


def _binomial(n: int, k: int) -> int:
    from math import comb

    return comb(n, k) if 0 <= k <= n else 0


# Mean mass of a heavy atom in drug-like space, for converting a molecular
# weight ceiling into an atom count.
MEAN_HEAVY_ATOM_MASS = 13.5
# Heavy atoms a typical scaffold plus one linker per arm consumes. The arms
# only get what is left after the assembly around them is paid for.
SCAFFOLD_RESERVE = 6
LINKER_RESERVE_PER_ARM = 1.5


def arm_heavy_atom_budget(profile: TargetProductProfile, k: int) -> int:
    """Heavy atoms available to the arms, after the assembly takes its share.

    A molecular-weight ceiling is a whole-molecule constraint, but the
    Hamiltonian only chooses arms.  Handing the arms the entire budget lets
    the optimiser pick three ten-atom pharmacophores and then discover, after
    assembly adds a scaffold and three linkers, that the molecule is 550 Da
    against a 420 Da ceiling.  The reserve is subtracted up front instead.
    """
    ceiling = profile.property_window.molecular_weight[1] / MEAN_HEAVY_ATOM_MASS
    available = ceiling - SCAFFOLD_RESERVE - LINKER_RESERVE_PER_ARM * k
    return max(int(available), 4 * k)


def build_selection_problem(
    profile: TargetProductProfile,
    library: FragmentLibrary,
    k: int = 3,
    heavy_atom_budget: int | None = None,
    weights: SelectionWeights | None = None,
    max_variables: int = 12,
) -> FragmentSelectionProblem:
    """Build the QUBO for choosing ``k`` pharmacophores against a profile.

    The library is pre-filtered to the ``max_variables`` fragments with the
    highest solo benefit.  That keeps the Hamiltonian inside the range where
    QAOA can be simulated and exact enumeration stays instant -- and it is a
    real restriction: a fragment that is weak alone but excellent in
    combination can be filtered out before the optimiser ever sees it.
    """
    weights = weights or SelectionWeights()
    if heavy_atom_budget is None:
        heavy_atom_budget = arm_heavy_atom_budget(profile, k)
    candidates = list(library.pharmacophores)

    share = max(heavy_atom_budget / max(k, 1), 1.0)
    window = profile.property_window
    tpsa_budget = max(window.tpsa[1], 1.0)
    hbd_budget = max(window.hbd[1], 1.0)
    logp_budget = max(window.clogp[1], 1.0)
    tpsa_share = tpsa_budget / max(k, 1)
    hbd_share = hbd_budget / max(k, 1)
    logp_share = logp_budget / max(k, 1)

    def envelope_excess(fragment: Fragment) -> float:
        """Cost of a fragment overrunning its share of the property envelope.

        Charged on all three near-additive properties, because charging on
        only some of them just moves the failure.  Penalising polarity alone
        made the optimiser swap a statin acid and a sulfonylurea -- too polar
        to cross the barrier -- for an adamantane and a chloroarene, landing
        at cLogP above 6, which fails the same gate from the other side.
        Polar surface area, donor count, and lipophilicity are each charged
        against a ``1/k`` share of the budget.
        """
        descriptors = fragment.capped_descriptors()
        over_tpsa = max(0.0, descriptors.tpsa - tpsa_share) / tpsa_budget
        over_hbd = max(0.0, descriptors.hbd - hbd_share) / hbd_budget
        over_logp = max(0.0, descriptors.clogp - logp_share) / logp_budget
        return over_tpsa + over_hbd + over_logp

    def mass_excess(fragment: Fragment) -> float:
        """Cost of a fragment overrunning its share of the size budget.

        Zero for anything at or under its share: a design is not rewarded for
        being small, only penalised for being too large to fit the envelope.
        """
        return max(0.0, fragment.heavy_atoms - share) / heavy_atom_budget

    def solo_benefit(fragment: Fragment) -> float:
        return (
            profile_coverage(fragment.engages, profile, combined_confidence([fragment]))
            - weights.mass_weight * mass_excess(fragment)
            - weights.envelope_weight * envelope_excess(fragment)
            - weights.liability_weight * len(fragment.liabilities)
        )

    ranked = sorted(candidates, key=lambda f: (-solo_benefit(f), f.identifier))
    selected = tuple(ranked[:max_variables])

    linear = {f.identifier: round(solo_benefit(f), 6) for f in selected}
    quadratic: dict[tuple[str, str], float] = {}
    for a, b in combinations(selected, 2):
        joint = profile_coverage(combined_effects([a, b]), profile, combined_confidence([a, b]))
        interaction = (
            joint
            - profile_coverage(a.engages, profile, combined_confidence([a]))
            - profile_coverage(b.engages, profile, combined_confidence([b]))
        )
        value = interaction - weights.redundancy_weight * redundancy(a, b, weights)
        # Keyed on the sorted pair. The fragments arrive in ranked order, so
        # storing them in that order and looking them up sorted silently
        # misses -- every coupling whose ranking disagrees with the alphabet
        # reads back as zero, and the optimiser solves a different problem
        # than the one reported.
        quadratic[tuple(sorted((a.identifier, b.identifier)))] = round(value, 6)

    return FragmentSelectionProblem(
        fragments=selected,
        linear=linear,
        quadratic=quadratic,
        k=k,
        profile=profile,
        weights=weights,
        heavy_atom_budget=heavy_atom_budget,
    )


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FragmentSelection:
    """One solved selection, with the backend that produced it.

    ``feasible`` is not decoration.  A heuristic solver can return a state
    that breaks the cardinality constraint, and such a state scores *better*
    on the objective precisely because it is buying coverage with arms it was
    not allowed to spend.  Reading the objective without reading this flag
    inverts the result.
    """

    identifiers: tuple[str, ...]
    qubo_objective: float
    exact_objective: float
    backend: str
    feasible: bool = True
    detail: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "fragments": list(self.identifiers),
            "n_fragments": len(self.identifiers),
            "qubo_objective": self.qubo_objective,
            "exact_objective": self.exact_objective,
            "backend": self.backend,
            "feasible": self.feasible,
            "detail": self.detail,
        }


def _selection(problem: FragmentSelectionProblem, ids: Iterable[str], backend: str,
               detail: dict[str, Any] | None = None) -> FragmentSelection:
    ordered = tuple(sorted(ids))
    return FragmentSelection(
        identifiers=ordered,
        qubo_objective=problem.qubo_objective(ordered),
        exact_objective=problem.exact_objective(ordered),
        backend=backend,
        feasible=len(ordered) == problem.k,
        detail=detail or {},
    )


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

def solve_enumeration(problem: FragmentSelectionProblem, top: int = 1) -> list[FragmentSelection]:
    """Exhaustive search over feasible subsets: the ground truth."""
    results = [
        _selection(problem, subset, "enumeration")
        for subset in combinations(problem.identifiers, problem.k)
    ]
    results.sort(key=lambda s: (-s.qubo_objective, s.identifiers))
    return results[:top]


def to_quadratic_program(problem: FragmentSelectionProblem):
    """Express the QUBO as a Qiskit ``QuadraticProgram``.

    This is the handover point to quantum hardware: the returned program is
    what a device-facing stack converts to an Ising Hamiltonian.
    """
    if not HAS_QISKIT:
        raise ImportError(
            "Qiskit is not installed in this environment. The enumeration backend solves "
            "the same problem; install the quantum stack from requirements.txt to use "
            "the eigensolver and QAOA backends."
        )
    program = QuadraticProgram(name="pharmacophore_selection")
    for identifier in problem.identifiers:
        program.binary_var(name=identifier)
    program.maximize(
        linear=dict(problem.linear),
        quadratic={(a, b): v for (a, b), v in problem.quadratic.items()},
    )
    program.linear_constraint(
        linear={identifier: 1 for identifier in problem.identifiers},
        sense="==",
        rhs=problem.k,
        name="arm_count",
    )
    return program


def _read_solution(problem: FragmentSelectionProblem, result) -> list[str]:
    return [
        identifier
        for identifier, value in zip(problem.identifiers, result.x)
        if round(float(value)) == 1
    ]


def solve_eigensolver(problem: FragmentSelectionProblem) -> FragmentSelection:
    """Exact classical diagonalisation of the Ising Hamiltonian, via Qiskit."""
    from qiskit_algorithms.minimum_eigensolvers import NumPyMinimumEigensolver

    program = to_quadratic_program(problem)
    optimizer = MinimumEigenOptimizer(NumPyMinimumEigensolver())
    result = optimizer.solve(program)
    return _selection(
        problem, _read_solution(problem, result), "eigensolver",
        {"status": str(result.status), "reported_value": float(result.fval)},
    )


def build_qaoa_ansatz(problem: FragmentSelectionProblem, reps: int = 2):
    """Construct the QAOA circuit for this problem's Ising Hamiltonian.

    The cardinality constraint is folded into the cost operator as a penalty
    by ``QuadraticProgramToQubo``, so the circuit optimises an unconstrained
    Ising model whose ground state is the constrained optimum.

    Two mechanical details matter and are easy to get wrong.  The ansatz is
    decomposed to elementary gates, because a sampler primitive rejects the
    composite ``QAOA`` instruction it is otherwise built from; and it is
    given explicit measurements, because the V2 sampler interface samples
    registers rather than statevectors.
    """
    from qiskit.circuit.library import QAOAAnsatz
    from qiskit_optimization.converters import QuadraticProgramToQubo

    program = to_quadratic_program(problem)
    qubo = QuadraticProgramToQubo().convert(program)
    cost_operator, offset = qubo.to_ising()
    ansatz = QAOAAnsatz(cost_operator=cost_operator, reps=reps).decompose(reps=5)
    ansatz.measure_all()
    return ansatz, cost_operator, offset


def _make_sampler(seed: int, shots: int) -> tuple[Any, str]:
    """The best available V2 sampler primitive.

    Qiskit's algorithm layer emits V2 primitive pubs, so a V1 sampler fails
    here with a type error that says nothing about the cause; the fallback
    chain exists to keep that failure from looking like a problem with the
    chemistry.
    """
    try:
        from qiskit_aer.primitives import SamplerV2

        return SamplerV2(seed=seed, default_shots=shots), "qiskit_aer.primitives.SamplerV2"
    except Exception:  # pragma: no cover - Aer availability
        from qiskit.primitives import StatevectorSampler

        return StatevectorSampler(seed=seed, default_shots=shots), "qiskit.primitives.StatevectorSampler"


def solve_qaoa(
    problem: FragmentSelectionProblem,
    reps: int = 3,
    seed: int = 7,
    maxiter: int = 400,
    shots: int = 4096,
) -> FragmentSelection:
    """Solve the same Hamiltonian with QAOA on a simulator.

    Seeded end to end -- sampler and classical optimiser -- so a run is
    reproducible.  QAOA is a heuristic: it can return a feasible but
    sub-optimal selection, which is why the result is always reported next to
    the enumeration optimum rather than in place of it.
    """
    from qiskit_algorithms.minimum_eigensolvers import SamplingVQE
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit_algorithms.utils import algorithm_globals

    algorithm_globals.random_seed = seed
    ansatz, cost_operator, offset = build_qaoa_ansatz(problem, reps=reps)
    sampler, sampler_name = _make_sampler(seed, shots)

    solver = SamplingVQE(sampler=sampler, ansatz=ansatz, optimizer=COBYLA(maxiter=maxiter))
    result = MinimumEigenOptimizer(solver).solve(to_quadratic_program(problem))
    return _selection(
        problem, _read_solution(problem, result), "qaoa",
        {
            "reps": reps, "seed": seed, "shots": shots, "maxiter": maxiter,
            "sampler": sampler_name,
            "status": str(result.status),
            "reported_value": float(result.fval),
            "qubits": cost_operator.num_qubits,
            "ansatz_depth": ansatz.depth(),
            "variational_parameters": ansatz.num_parameters,
            "ising_offset": round(float(offset), 6),
        },
    )


def solve_qaoa_sweep(
    problem: FragmentSelectionProblem,
    reps_range: Sequence[int] = (1, 2, 3, 4, 5),
    seed: int = 7,
    maxiter: int = 400,
    shots: int = 4096,
) -> FragmentSelection:
    """Run QAOA across several circuit depths and keep the best feasible result.

    A single fixed depth is not a defensible default. Measured across the two
    registered diseases, QAOA's success is **non-monotonic in depth and
    specific to the problem instance**: depth 3 solved one campaign and failed
    both after the fragment library grew, while depths 2, 4, and 5 each
    succeeded on some instances and not others. Deeper circuits have more
    variational parameters and a harder classical landscape, so they do not
    reliably improve on shallower ones at a fixed iteration budget.

    Reporting one lucky depth would present a heuristic as though it were a
    solver. The sweep keeps the best *feasible* selection and records what
    every depth returned, so the variance is visible in the output.
    """
    attempts: list[dict[str, Any]] = []
    best: FragmentSelection | None = None
    for reps in reps_range:
        try:
            candidate = solve_qaoa(problem, reps=reps, seed=seed, maxiter=maxiter, shots=shots)
        except Exception as error:  # pragma: no cover - backend availability
            attempts.append({"reps": reps, "error": f"{type(error).__name__}: {error}"})
            continue
        attempts.append({
            "reps": reps,
            "feasible": candidate.feasible,
            "qubo_objective": candidate.qubo_objective,
            "fragments": list(candidate.identifiers),
        })
        if not candidate.feasible:
            continue
        if best is None or candidate.qubo_objective > best.qubo_objective:
            best = candidate
    if best is None:
        raise RuntimeError(
            f"QAOA returned no feasible selection at any depth in {list(reps_range)}"
        )
    detail = dict(best.detail)
    detail["depth_sweep"] = attempts
    detail["depths_tried"] = list(reps_range)
    detail["depths_feasible"] = [a["reps"] for a in attempts if a.get("feasible")]
    return FragmentSelection(
        identifiers=best.identifiers,
        qubo_objective=best.qubo_objective,
        exact_objective=best.exact_objective,
        backend="qaoa",
        feasible=best.feasible,
        detail=detail,
    )


def solve_selection(
    problem: FragmentSelectionProblem,
    backend: str = "enumeration",
    **kwargs: Any,
) -> FragmentSelection:
    """Solve with the named backend: ``enumeration``, ``eigensolver``, or ``qaoa``."""
    if backend == "enumeration":
        return solve_enumeration(problem, top=1)[0]
    if backend == "eigensolver":
        return solve_eigensolver(problem)
    if backend == "qaoa":
        return solve_qaoa_sweep(problem, **kwargs)
    raise ValueError(f"Unknown backend {backend!r}; expected enumeration, eigensolver, or qaoa")


def benchmark_backends(
    problem: FragmentSelectionProblem,
    backends: Sequence[str] = ("enumeration", "eigensolver", "qaoa"),
    **kwargs: Any,
) -> dict[str, Any]:
    """Run several backends on one Hamiltonian and report whether they agree.

    Agreement against the enumeration optimum is the only claim being made:
    that the QUBO formulation is faithful and transfers to a quantum
    algorithm.  A disagreement is a QAOA convergence failure, not evidence
    about the chemistry.
    """
    import time

    reference = solve_enumeration(problem, top=1)[0]
    outcomes: dict[str, Any] = {}
    for backend in backends:
        if backend != "enumeration" and not HAS_QISKIT:
            outcomes[backend] = {"skipped": "qiskit not installed in this environment"}
            continue
        start = time.perf_counter()
        try:
            selection = solve_selection(problem, backend, **(kwargs if backend == "qaoa" else {}))
        except Exception as error:  # pragma: no cover - backend availability
            outcomes[backend] = {"error": f"{type(error).__name__}: {error}"}
            continue
        elapsed = time.perf_counter() - start
        outcomes[backend] = {
            **selection.as_dict(),
            "seconds": round(elapsed, 3),
            "matches_enumeration_optimum": (
                selection.feasible and selection.identifiers == reference.identifiers
            ),
            "objective_gap": round(reference.qubo_objective - selection.qubo_objective, 6),
            "note": (
                ""
                if selection.feasible
                else f"Constraint violated: {len(selection.identifiers)} arms selected where "
                     f"{problem.k} were allowed. The objective is higher than the true optimum "
                     f"only because the extra arm was not paid for; this is a QAOA convergence "
                     f"failure, not a better design."
            ),
        }
    return {
        "hamiltonian": problem.hamiltonian_summary(),
        "weights": problem.weights.as_dict(),
        "reference_optimum": reference.as_dict(),
        "backends": outcomes,
    }

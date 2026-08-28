"""Disease-agnostic scoring of drug *combinations of any order*, k = 1, 2, 3, ...

Two things separate this module from ``core.biology.ms_scoring``, which it
generalises.

**It takes its vocabulary from the disease, not from itself.** ``ms_scoring``
holds MS's pathways, therapeutic axes, risk domains, and risk weights as module
constants, so adding a second disease meant editing scoring logic. Everything
here is read from a :class:`~core.models.disease.DiseaseContext`, which is a
registry entry under ``data/diseases/``. Adding a third disease is a data task.

**It scores subsets, not pairs.** ``pair_metrics`` is the k = 2 case of
:func:`combination_metrics`. Monotherapy (k = 1) is scored through the same
code path as a triple, which is the only way a monotherapy-versus-combination
comparison means anything: if singles and pairs were scored by different
functions, any difference between them would be partly an artefact of the two
implementations.

Comparing across k, carefully
-----------------------------
``priority_score`` is a weighted composite that includes terms only a
combination can earn -- target complementarity, compartment complementarity,
network separation. A single agent scores zero on those by construction, so
**ranking a single against a pair by ``priority_score`` is not a fair
comparison and this module does not ask anyone to do it.** Use it to rank
within a fixed k.

Across k, compare on the efficacy block, which is defined identically at every
order: ``signed_reversal``, ``reversal_efficiency``, ``gene_coverage``,
``pathway_coverage``, ``axis_coverage``, ``counter_therapeutic``. The
comparison the screen actually exists to make is carried by three fields:

* ``reversal_gain_over_best_single`` -- what the combination adds over its own
  best member;
* ``score_gain_over_best_subset`` -- for k >= 3, what the k-th agent adds over
  the best (k-1)-subset it contains, i.e. whether it earns its place;
* ``additivity_ratio`` -- combined reversal divided by the sum of the members'
  solo reversals. Below 1 the members are covering the same signal; near 1 they
  are independent.

Nothing here is a dose-response measurement, a synergy estimate, or a
co-administration recommendation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations
from typing import Any, Iterable, Sequence

import numpy as np

from core.biology.signature import (
    EVIDENCE_RANK, Signature, alignment_metrics, combine_effects,
)
from core.models.disease import DiseaseContext

# A combination touching this fraction of total signature weight is treated as
# broad multi-pathway coverage. Declared explicitly so the coverage term lands
# on [0, 1] rather than being crushed by a 90-to-112-gene denominator.
BROAD_COVERAGE_REFERENCE = 0.25

# Ordinal distance between half-life classes; washout planning and the
# interpretation of an adverse event both depend on the gap.
_HALF_LIFE_ORDER = {"short": 0, "intermediate": 1, "long": 2, "reconstitution": 3}

# How much each route costs a chronic regimen. Unknown routes are treated as
# mid-burden rather than silently free.
_ROUTE_BURDEN = {
    "oral": 0.0, "transdermal": 0.15, "inhaled": 0.2, "subcutaneous": 0.4,
    "intramuscular": 0.5, "infusion": 0.7, "intrathecal": 0.9, "intraputaminal": 1.0,
}
_DEFAULT_ROUTE_BURDEN = 0.5

# Where an agent acts, as a scalar. Used for compartment complementarity.
_COMPARTMENT_SCALE = {"peripheral": 0.0, "both": 0.5, "cns": 1.0}


@dataclass(frozen=True)
class CombinationConfig:
    """Declared weights for the composite prioritisation score.

    Every weight is named, defaulted, and serialised into the result file. The
    accompanying Pareto analysis reports which combinations are non-dominated
    irrespective of how these are set, which is the check on the fact that they
    were chosen by hand.
    """

    reversal_efficiency_weight: float = 1.00
    coverage_weight: float = 0.45
    pathway_weight: float = 0.20
    axis_weight: float = 0.25
    complementarity_weight: float = 0.15
    compartment_weight: float = 0.20
    cns_weight: float = 0.15
    separation_weight: float = 0.20
    evidence_weight: float = 0.15
    counter_therapeutic_penalty: float = 0.50
    safety_penalty: float = 0.35
    safety_overlap_penalty: float = 0.30
    uncertainty_penalty: float = 0.20
    regimen_penalty: float = 0.15
    minimum_evidence: str = "phase_2"
    # Two agents can carry different mechanism_class labels and still hit the
    # same biology -- an anti-alpha-synuclein antibody and an aggregation
    # inhibitor, say. Target-set overlap above this fraction is treated as
    # pharmacodynamic redundancy regardless of how the classes are named.
    max_target_overlap: float = 0.6
    broad_coverage_reference: float = BROAD_COVERAGE_REFERENCE

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Set-level helpers, each reducing to the familiar pairwise form at k = 2
# ---------------------------------------------------------------------------

def _mean_pairwise(members: Sequence[Any], fn) -> float:
    """Mean of ``fn`` over every unordered pair; 0.0 for a single member.

    Every k-ary term in this module is defined as the mean of its pairwise
    definition. That keeps k = 2 numerically identical to the pairwise screen
    while giving k >= 3 a definition that does not privilege any member.
    """
    pairs = list(combinations(members, 2))
    if not pairs:
        return 0.0
    return float(np.mean([fn(a, b) for a, b in pairs]))


def combine_many(effect_maps: Iterable[dict[str, float]]) -> dict[str, float]:
    """Fold any number of signed effect vectors together under Bliss."""
    combined: dict[str, float] = {}
    for effects in effect_maps:
        combined = combine_effects(combined, effects)
    return combined


def _target_complementarity(drugs: Sequence[dict]) -> float:
    """Mean pairwise Jaccard distance between target sets."""
    def jaccard_distance(a: dict, b: dict) -> float:
        ga, gb = set(a["target_effects"]), set(b["target_effects"])
        return 1.0 - len(ga & gb) / max(1, len(ga | gb))

    return _mean_pairwise(drugs, jaccard_distance)


def _compartment_complementarity(drugs: Sequence[dict], requires_cns: bool) -> float:
    """Reward covering the compartments where the disease actually lives.

    For a disease whose pathology sits behind the blood-brain barrier, a
    regimen that never reaches the CNS cannot modify it however well it scores
    on transcriptional reversal. Where the registry declares
    ``requires_cns_exposure``, this term is CNS coverage; where it does not, it
    is compartment *spread*, which is what the MS screen rewards for the
    peripheral/central split of that disease.
    """
    values = [_COMPARTMENT_SCALE[d["compartment"]] for d in drugs]
    if requires_cns:
        # Best achievable is at least one agent fully central; a second central
        # agent adds nothing, a peripheral partner is not penalised here (it is
        # priced by the safety and regimen terms instead).
        return float(max(values))
    return _mean_pairwise(values, lambda a, b: abs(a - b) + min(a, b) * (1.0 - abs(a - b)))


def _axis_coverage(drugs: Sequence[dict], axes: Sequence[str]) -> float:
    covered = {a for d in drugs for a in d.get("therapeutic_axes", [])}
    return len(covered & set(axes)) / len(axes)


def _pathway_coverage(effects: dict[str, float], signature: Signature,
                      pathways: Sequence[str]) -> float:
    hit = {signature.pathway[g] for g in effects
           if g in signature.pathway and signature.pathway[g]}
    return len(hit & set(pathways)) / len(pathways)


def _regimen_burden(drugs: Sequence[dict]) -> float:
    """Practical cost of running this regimen: route burden and half-life spread.

    Reported as a burden on [0, 1] -- higher is worse -- because it enters the
    composite as a penalty. A four-drug all-oral regimen is cheap here; adding
    an infused or intraputaminal agent is not.
    """
    routes = [_ROUTE_BURDEN.get(d["route"], _DEFAULT_ROUTE_BURDEN) for d in drugs]
    order = [_HALF_LIFE_ORDER[d["half_life_class"]] for d in drugs]
    spread = (max(order) - min(order)) / 3.0 if len(order) > 1 else 0.0
    # Route burden is the max, not the mean: one infusion sets the ceiling on
    # how convenient the whole regimen can be.
    return float(0.6 * max(routes) + 0.4 * spread)


def combined_safety(drugs: Sequence[dict], disease: DiseaseContext) -> dict[str, Any]:
    """Aggregate safety burden across a combination, per the disease's domains.

    Per domain the union risk is ``1 - prod(1 - a_i)``; ``overlap`` is the mean
    pairwise product, which isolates the domains where more than one agent
    contributes -- the case where combination use is most likely to be
    additive in a clinically meaningful way. Domain weights come from the
    registry entry, because how much a domain constrains a regimen is a
    property of the disease and its population, not of a scoring function.
    """
    weights = disease.risk_weights
    weight_sum = sum(weights.values()) or 1.0
    union_total, overlap_total, worst, worst_domain = 0.0, 0.0, 0.0, ""
    per_domain: dict[str, float] = {}
    for domain in disease.risk_domains:
        values = [float(d.get("safety_burden", {}).get(domain, 0.0)) for d in drugs]
        w = weights.get(domain, 1.0)
        union = 1.0 - float(np.prod([1.0 - v for v in values]))
        overlap = _mean_pairwise(values, lambda a, b: a * b)
        per_domain[domain] = round(union, 4)
        union_total += w * union
        overlap_total += w * overlap
        if union * w > worst:
            worst, worst_domain = union * w, domain
    return {
        "safety_union": union_total / weight_sum,
        "safety_overlap": overlap_total / weight_sum,
        "worst_risk_domain": worst_domain,
        "safety_by_domain": per_domain,
    }


# ---------------------------------------------------------------------------
# Redundancy
# ---------------------------------------------------------------------------

def redundancy_flags(drugs: Sequence[dict], config: CombinationConfig) -> dict[str, bool]:
    """Flag a combination whose members duplicate one another.

    A combination is redundant if *any* pair inside it is redundant, which is
    what makes triples buildable only from mutually distinct agents. Duplicated
    members are computed and written to the full CSV but held out of the
    primary ranking, so the exclusion is auditable rather than invisible.
    """
    same_mechanism = shared_safety = redundant_targets = False
    for a, b in combinations(drugs, 2):
        same_mechanism |= a["mechanism_class"] == b["mechanism_class"]
        shared_safety |= bool(set(a["safety_classes"]) & set(b["safety_classes"]))
        family_a = a.get("target_family", a["mechanism_class"])
        family_b = b.get("target_family", b["mechanism_class"])
        ga, gb = set(a["target_effects"]), set(b["target_effects"])
        overlap = len(ga & gb) / max(1, len(ga | gb))
        redundant_targets |= (family_a == family_b) or overlap > config.max_target_overlap
    return {
        "same_mechanism": same_mechanism,
        "shared_safety_class": shared_safety,
        "redundant_targets": redundant_targets,
        "excluded_from_primary_ranking": same_mechanism or shared_safety or redundant_targets,
    }


# ---------------------------------------------------------------------------
# The scorer
# ---------------------------------------------------------------------------

def combination_metrics(
    drugs: Sequence[dict[str, Any]],
    signature: Signature,
    disease: DiseaseContext,
    config: CombinationConfig,
    exposure: Any = None,
) -> dict[str, Any]:
    """Score one combination of any order across every declared parameter.

    ``exposure`` is an optional :class:`~core.biology.network_proximity.
    ComplementaryExposure`; when omitted the topology terms are reported as NaN
    and contribute nothing to the composite, rather than being imputed.
    """
    if not drugs:
        raise ValueError("combination_metrics requires at least one agent")
    names = [d["name"] for d in drugs]
    order = len(drugs)

    combined = combine_many(d["target_effects"] for d in drugs)
    combo = alignment_metrics(combined, signature)
    solo = [alignment_metrics(d["target_effects"], signature) for d in drugs]
    solo_reversals = [s["reversal"] for s in solo]
    best_single = max(solo_reversals)
    solo_sum = sum(solo_reversals)

    complementarity = _target_complementarity(drugs)
    compartment = _compartment_complementarity(drugs, disease.delivery.requires_cns_exposure)
    axis = _axis_coverage(drugs, disease.therapeutic_axes)
    pathway = _pathway_coverage(combined, signature, disease.pathways)
    regimen_burden = _regimen_burden(drugs)
    safety = combined_safety(drugs, disease)

    evidence = min(EVIDENCE_RANK[d["evidence_tier"]] for d in drugs) / 3.0
    uncertainty = float(np.mean([float(d["target_uncertainty"]) for d in drugs]))
    cns_reach = max(float(d["cns_penetration"]) for d in drugs)

    separation = float(exposure.separation) if exposure is not None else float("nan")
    # Only positive separation earns credit; overlap earns nothing, not a penalty.
    separation_term = (
        0.0 if (exposure is None or np.isnan(separation))
        else float(np.clip(separation, 0.0, 2.0) / 2.0)
    )

    coverage_scaled = float(min(1.0, combo["gene_coverage"] / config.broad_coverage_reference))
    counter_rate = (
        combo["counter_therapeutic"] / combo["gene_coverage"]
        if combo["gene_coverage"] > 0 else 0.0
    )

    priority = (
        config.reversal_efficiency_weight * combo["reversal_efficiency"]
        + config.coverage_weight * coverage_scaled
        + config.pathway_weight * pathway
        + config.axis_weight * axis
        + config.complementarity_weight * complementarity
        + config.compartment_weight * compartment
        + config.cns_weight * cns_reach
        + config.separation_weight * separation_term
        + config.evidence_weight * evidence
        - config.counter_therapeutic_penalty * counter_rate
        - config.safety_penalty * safety["safety_union"]
        - config.safety_overlap_penalty * safety["safety_overlap"]
        - config.uncertainty_penalty * uncertainty
        - config.regimen_penalty * regimen_burden
    )

    row = {
        "combination": " + ".join(names),
        "members": names,
        "order": order,
        "priority_score": round(float(priority), 4),
        # --- efficacy block: defined identically at every order, and the only
        #     block that may be compared across k -----------------------------
        "signed_reversal": round(combo["reversal"], 4),
        "reversal_efficiency": round(combo["reversal_efficiency"], 4),
        "gene_coverage": round(combo["gene_coverage"], 4),
        "pathway_coverage": round(pathway, 4),
        "axis_coverage": round(axis, 4),
        "counter_therapeutic": round(combo["counter_therapeutic"], 4),
        # --- efficacy comparison ------------------------------------------
        "best_single_reversal": round(best_single, 4),
        "reversal_gain_over_best_single": round(combo["reversal"] - best_single, 4),
        # A single agent is trivially additive with itself. Without this, an
        # agent whose every effect is counter-therapeutic (carbidopa, whose only
        # modelled action is DDC inhibition) has reversal 0, and 0/0 would be
        # reported as total redundancy rather than as no signal.
        "additivity_ratio": (
            1.0 if order == 1
            else round(float(combo["reversal"] / solo_sum), 4) if solo_sum > 0 else 0.0
        ),
        # --- pharmacology / delivery --------------------------------------
        "target_complementarity": round(complementarity, 4),
        "compartment_complementarity": round(compartment, 4),
        "cns_reach": round(cns_reach, 4),
        "regimen_burden": round(regimen_burden, 4),
        # --- network topology ---------------------------------------------
        "network_separation": separation,
        "proximity_z_min": float(exposure.z_a) if exposure is not None else float("nan"),
        "proximity_z_max": float(exposure.z_b) if exposure is not None else float("nan"),
        "complementary_exposure": bool(exposure.is_complementary) if exposure is not None else False,
        # --- evidence / risk ----------------------------------------------
        "evidence_score": round(evidence, 4),
        "uncertainty": round(uncertainty, 4),
        "safety_union": round(safety["safety_union"], 4),
        "safety_overlap": round(safety["safety_overlap"], 4),
        "worst_risk_domain": safety["worst_risk_domain"],
    }
    row.update(redundancy_flags(drugs, config))
    return row


# Objectives for the weight-free Pareto analysis: (key, maximise?)
PARETO_OBJECTIVES: tuple[tuple[str, bool], ...] = (
    ("reversal_efficiency", True),
    ("gene_coverage", True),
    ("axis_coverage", True),
    ("compartment_complementarity", True),
    ("evidence_score", True),
    ("counter_therapeutic", False),
    ("safety_union", False),
    ("uncertainty", False),
    ("regimen_burden", False),
)


def eligible_drugs(drugs: Sequence[dict], config: CombinationConfig) -> list[dict]:
    floor = EVIDENCE_RANK[config.minimum_evidence]
    return [d for d in drugs if EVIDENCE_RANK[d["evidence_tier"]] >= floor]


def combination_key(names: Iterable[str]) -> tuple[str, ...]:
    """Order-independent identity for a combination."""
    return tuple(sorted(names))


def rank_combinations(
    drugs: Sequence[dict[str, Any]],
    signature: Signature,
    disease: DiseaseContext,
    config: CombinationConfig | None = None,
    order: int = 2,
    exposures: dict[tuple[str, ...], Any] | None = None,
    prefilter: set[tuple[str, ...]] | None = None,
) -> list[dict[str, Any]]:
    """Score every combination of the given order and sort deterministically.

    ``prefilter``, when given, restricts which combinations are built: only
    subsets all of whose (order-1) faces appear in it are scored. At order 3
    the unrestricted space is large and mostly uninteresting -- the vast
    majority of triples contain a pair that was already excluded as redundant
    -- so the runner passes the surviving pairs and the triple enumeration
    inherits their exclusions instead of rediscovering them.
    """
    config = config or CombinationConfig()
    exposures = exposures or {}
    pool = eligible_drugs(drugs, config)
    rows = []
    for members in combinations(pool, order):
        key = combination_key(d["name"] for d in members)
        if prefilter is not None and order > 1:
            faces = [combination_key(f) for f in combinations(key, order - 1)]
            if not all(face in prefilter for face in faces):
                continue
        rows.append(combination_metrics(members, signature, disease, config, exposures.get(key)))
    return sort_rows(rows)


def sort_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda r: (r["excluded_from_primary_ranking"], -r["priority_score"], r["combination"]),
    )


def attach_subset_gain(
    rows: Sequence[dict[str, Any]], lower_order: Sequence[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Record what the k-th agent adds over the best (k-1)-subset it contains.

    This is the question a triple has to answer that a pair does not: adding a
    third chronic medication to an elderly patient's regimen has to buy
    something, and ``score_gain_over_best_subset <= 0`` says it does not.
    """
    lookup = {combination_key(r["members"]): r["priority_score"] for r in lower_order}
    for row in rows:
        faces = [combination_key(f) for f in combinations(row["members"], row["order"] - 1)]
        scores = [lookup[f] for f in faces if f in lookup]
        if scores:
            best = max(scores)
            row["best_subset_score"] = round(float(best), 4)
            row["score_gain_over_best_subset"] = round(float(row["priority_score"] - best), 4)
        else:
            row["best_subset_score"] = float("nan")
            row["score_gain_over_best_subset"] = float("nan")
    return list(rows)

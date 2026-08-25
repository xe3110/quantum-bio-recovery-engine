"""Extended multi-parameter prioritisation of MS drug-combination candidates.

Design commitments
------------------
1. **Direction-aware.** The disease signature carries an explicit
   ``desired_direction`` per gene.  Raw log fold-change is the wrong target for
   compensatory transcripts (e.g. ``HMOX1`` rises in MS as an antioxidant
   response), so reversing it would penalise Nrf2 activators.  Scoring uses the
   curated therapeutic direction instead.
2. **Bliss independence, not clipping.** Two agents acting on one gene in the
   same direction combine as ``a + b - ab`` on the fractional scale, which
   saturates at 1 without the discontinuity that hard clipping introduces.
   Opposing effects add, so mutual antagonism is represented rather than hidden.
3. **Penalties are first-class.** Counter-therapeutic movement, overlapping
   safety-risk domains, and target uncertainty each reduce the score directly.
4. **Weights are declared and auditable.** ``ScoringConfig`` holds every weight,
   and the accompanying Pareto analysis reports which pairs are non-dominated
   irrespective of weighting.

Nothing here is a dose-response measurement, a synergy estimate, or a
co-administration recommendation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field, asdict
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

# The generic contracts now live in core.biology.signature so that campaigns
# for other diseases can use them without importing MS scoring. Re-exported
# here so existing callers keep working unchanged.
from core.biology.signature import (  # noqa: F401
    EVIDENCE_RANK, Signature, alignment_metrics, bliss_combine, combine_effects,
    load_signature,
)

MS_DISEASE_PATHWAYS = (
    "adaptive_immunity", "B_cell_immunity", "innate_immunity",
    "microglial_activation", "lymphocyte_trafficking", "blood_brain_barrier",
    "oxidative_stress", "remyelination", "neuroprotection", "metabolism",
)

THERAPEUTIC_AXES = (
    "immunomodulation", "cns_innate", "remyelination", "neuroprotection", "metabolic_repair",
)

RISK_DOMAINS = (
    "infection", "malignancy", "cardiac", "hepatic", "teratogenicity", "autoimmunity", "ocular",
)

# A pair touching this fraction of total signature weight is treated as broad
# multi-pathway coverage. Declared explicitly so the coverage term is on [0,1]
# rather than being crushed by the 112-gene denominator.
BROAD_COVERAGE_REFERENCE = 0.25

# Half-life mismatch matters for washout planning and for interpreting an
# adverse event, so an ordinal distance between classes feeds the score.
_HALF_LIFE_ORDER = {"short": 0, "intermediate": 1, "long": 2, "reconstitution": 3}


@dataclass(frozen=True)
class ScoringConfig:
    """Declared weights for the composite prioritisation score."""

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
    minimum_evidence: str = "phase_2"
    # Two agents can carry different mechanism_class labels yet hit the same
    # biological target (an anti-alpha4 antibody and an oral alpha4 antagonist,
    # say). Target-set overlap above this fraction is treated as
    # pharmacodynamic redundancy regardless of how the classes are named.
    max_target_overlap: float = 0.6
    # Risk domains weighted by how strongly they constrain combination use.
    risk_weights: dict[str, float] = field(default_factory=lambda: {
        "infection": 1.0, "malignancy": 0.9, "cardiac": 0.7,
        "hepatic": 0.7, "teratogenicity": 0.5, "autoimmunity": 0.9, "ocular": 0.4,
    })

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _pathway_coverage(effects: dict[str, float], signature: Signature) -> float:
    hit = {signature.pathway[g] for g in effects if g in signature.pathway and signature.pathway[g]}
    return len(hit & set(MS_DISEASE_PATHWAYS)) / len(MS_DISEASE_PATHWAYS)


def combined_safety(drug_a: dict, drug_b: dict, config: ScoringConfig) -> dict[str, float]:
    """Aggregate pairwise safety burden and the overlap between risk domains.

    Per domain the union risk is ``1-(1-a)(1-b)``; ``overlap`` captures domains
    where *both* agents contribute, which is where combination use is most
    likely to be additive in a clinically meaningful way.
    """
    a_risk = drug_a.get("safety_burden", {})
    b_risk = drug_b.get("safety_burden", {})
    weights = config.risk_weights
    weight_sum = sum(weights.values())
    union, overlap, worst = 0.0, 0.0, 0.0
    worst_domain = ""
    for domain in RISK_DOMAINS:
        a = float(a_risk.get(domain, 0.0))
        b = float(b_risk.get(domain, 0.0))
        w = weights.get(domain, 1.0)
        combined = 1.0 - (1.0 - a) * (1.0 - b)
        union += w * combined
        overlap += w * (a * b)
        if combined * w > worst:
            worst, worst_domain = combined * w, domain
    return {
        "safety_union": union / weight_sum,
        "safety_overlap": overlap / weight_sum,
        "worst_risk_domain": worst_domain,
    }


def _compartment_complementarity(drug_a: dict, drug_b: dict) -> float:
    """Reward pairing a peripheral agent with one acting inside the CNS.

    Compartmentalised inflammation behind an intact blood-brain barrier is a
    leading explanation for why peripheral immunosuppression alone plateaus in
    progressive MS, so covering both compartments is scored explicitly.
    """
    scale = {"peripheral": 0.0, "both": 0.5, "cns": 1.0}
    a = scale[drug_a["compartment"]]
    b = scale[drug_b["compartment"]]
    return float(abs(a - b) + 0.5 * min(a, b) * 2 * (1 - abs(a - b)))


def _axis_coverage(drug_a: dict, drug_b: dict) -> float:
    axes = set(drug_a.get("therapeutic_axes", [])) | set(drug_b.get("therapeutic_axes", []))
    return len(axes & set(THERAPEUTIC_AXES)) / len(THERAPEUTIC_AXES)


def _regimen_compatibility(drug_a: dict, drug_b: dict) -> float:
    """Score practical co-administration: route burden and half-life mismatch."""
    oral = {drug_a["route"], drug_b["route"]} <= {"oral"}
    both_infusion = drug_a["route"] == drug_b["route"] == "infusion"
    route_score = 1.0 if oral else (0.4 if both_infusion else 0.7)
    gap = abs(_HALF_LIFE_ORDER[drug_a["half_life_class"]] - _HALF_LIFE_ORDER[drug_b["half_life_class"]])
    return round(0.6 * route_score + 0.4 * (1.0 - gap / 3.0), 4)


def pair_metrics(
    drug_a: dict[str, Any],
    drug_b: dict[str, Any],
    signature: Signature,
    config: ScoringConfig,
    exposure: Any = None,
) -> dict[str, Any]:
    """Score one candidate pair across every declared parameter.

    ``exposure`` is an optional :class:`ComplementaryExposure` from
    ``network_proximity``; when omitted the topology terms are reported as NaN
    and contribute nothing to the composite score.
    """
    combined = combine_effects(drug_a["target_effects"], drug_b["target_effects"])
    combo = alignment_metrics(combined, signature)
    solo_a = alignment_metrics(drug_a["target_effects"], signature)
    solo_b = alignment_metrics(drug_b["target_effects"], signature)
    best_single = max(solo_a["reversal"], solo_b["reversal"])

    genes_a, genes_b = set(drug_a["target_effects"]), set(drug_b["target_effects"])
    complementarity = 1.0 - len(genes_a & genes_b) / max(1, len(genes_a | genes_b))

    evidence = min(EVIDENCE_RANK[drug_a["evidence_tier"]], EVIDENCE_RANK[drug_b["evidence_tier"]]) / 3.0
    uncertainty = (float(drug_a["target_uncertainty"]) + float(drug_b["target_uncertainty"])) / 2.0
    safety = combined_safety(drug_a, drug_b, config)

    cns_reach = max(float(drug_a["cns_penetration"]), float(drug_b["cns_penetration"]))
    compartment = _compartment_complementarity(drug_a, drug_b)
    axis = _axis_coverage(drug_a, drug_b)
    pathway = _pathway_coverage(combined, signature)
    regimen = _regimen_compatibility(drug_a, drug_b)

    separation = float(exposure.separation) if exposure is not None else float("nan")
    # Only positive separation earns credit; overlap earns nothing, not a penalty.
    separation_term = 0.0 if (exposure is None or np.isnan(separation)) else float(np.clip(separation, 0.0, 2.0) / 2.0)

    same_mechanism = drug_a["mechanism_class"] == drug_b["mechanism_class"]
    shared_safety_class = bool(set(drug_a["safety_classes"]) & set(drug_b["safety_classes"]))
    same_family = drug_a.get("target_family", drug_a["mechanism_class"]) == drug_b.get(
        "target_family", drug_b["mechanism_class"]
    )
    redundant_targets = same_family or (1.0 - complementarity) > config.max_target_overlap
    excluded = same_mechanism or shared_safety_class or redundant_targets

    coverage_scaled = float(min(1.0, combo["gene_coverage"] / BROAD_COVERAGE_REFERENCE))
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
        - config.counter_therapeutic_penalty * (combo["counter_therapeutic"] / combo["gene_coverage"] if combo["gene_coverage"] > 0 else 0.0)
        - config.safety_penalty * safety["safety_union"]
        - config.safety_overlap_penalty * safety["safety_overlap"]
        - config.uncertainty_penalty * uncertainty
    )

    return {
        "drug_a": drug_a["name"],
        "drug_b": drug_b["name"],
        "priority_score": round(float(priority), 4),
        # --- efficacy-side parameters -------------------------------------
        "signed_reversal": round(combo["reversal"], 4),
        "reversal_efficiency": round(combo["reversal_efficiency"], 4),
        "coverage_scaled": round(coverage_scaled, 4),
        "best_single_reversal": round(best_single, 4),
        "reversal_gain_over_best_single": round(combo["reversal"] - best_single, 4),
        "counter_therapeutic": round(combo["counter_therapeutic"], 4),
        "gene_coverage": round(combo["gene_coverage"], 4),
        "pathway_coverage": round(pathway, 4),
        "axis_coverage": round(axis, 4),
        "target_complementarity": round(complementarity, 4),
        # --- pharmacology / delivery --------------------------------------
        "compartment_complementarity": round(compartment, 4),
        "cns_reach": round(cns_reach, 4),
        "regimen_compatibility": regimen,
        # --- network topology ---------------------------------------------
        "network_separation": separation,
        "proximity_z_a": float(exposure.z_a) if exposure is not None else float("nan"),
        "proximity_z_b": float(exposure.z_b) if exposure is not None else float("nan"),
        "complementary_exposure": bool(exposure.is_complementary) if exposure is not None else False,
        # --- evidence / risk ----------------------------------------------
        "evidence_score": round(evidence, 4),
        "uncertainty": round(uncertainty, 4),
        "safety_union": round(safety["safety_union"], 4),
        "safety_overlap": round(safety["safety_overlap"], 4),
        "worst_risk_domain": safety["worst_risk_domain"],
        # --- redundancy flags ----------------------------------------------
        "same_mechanism": same_mechanism,
        "shared_safety_class": shared_safety_class,
        "redundant_targets": redundant_targets,
        "excluded_from_primary_ranking": excluded,
    }


# Objectives used for the weight-free Pareto analysis: (key, maximise?)
PARETO_OBJECTIVES: tuple[tuple[str, bool], ...] = (
    ("reversal_efficiency", True),
    ("gene_coverage", True),
    ("axis_coverage", True),
    ("compartment_complementarity", True),
    ("evidence_score", True),
    ("counter_therapeutic", False),
    ("safety_union", False),
    ("uncertainty", False),
)


def eligible_drugs(drugs: Sequence[dict], config: ScoringConfig) -> list[dict]:
    floor = EVIDENCE_RANK[config.minimum_evidence]
    return [d for d in drugs if EVIDENCE_RANK[d["evidence_tier"]] >= floor]


def rank_pairs(
    drugs: Sequence[dict[str, Any]],
    signature: Signature,
    config: ScoringConfig | None = None,
    exposures: dict[tuple[str, str], Any] | None = None,
) -> list[dict[str, Any]]:
    """Score every eligible pair and sort deterministically."""
    config = config or ScoringConfig()
    exposures = exposures or {}
    rows = []
    for a, b in combinations(eligible_drugs(drugs, config), 2):
        key = (min(a["name"], b["name"]), max(a["name"], b["name"]))
        rows.append(pair_metrics(a, b, signature, config, exposures.get(key)))
    return sorted(
        rows,
        key=lambda r: (r["excluded_from_primary_ranking"], -r["priority_score"], r["drug_a"], r["drug_b"]),
    )

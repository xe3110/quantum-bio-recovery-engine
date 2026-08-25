"""Deterministic, research-use-only prioritisation for MS drug combinations.

This module deliberately separates *candidate generation* from claims of
clinical efficacy.  Scores are interpretable, reproducible heuristics derived
from a signed disease signature and a curated drug-target/pathway panel; they
are not dose-response measurements and must be validated experimentally.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, Iterable

import numpy as np


@dataclass(frozen=True)
class PrioritizationConfig:
    pathway_weight: float = 0.25
    complementarity_weight: float = 0.20
    evidence_weight: float = 0.15
    uncertainty_penalty: float = 0.20
    minimum_evidence: str = "preclinical"


EVIDENCE_RANK = {"approved": 3, "phase_3": 2, "phase_2": 1, "preclinical": 0}
MS_DISEASE_PATHWAYS = frozenset({
    "adaptive_immunity", "B_cell_immunity", "innate_immunity",
    "microglial_activation", "lymphocyte_trafficking", "blood_brain_barrier",
    "oxidative_stress", "remyelination", "neuroprotection", "metabolism",
})


def _clamp(value: float) -> float:
    return float(np.clip(value, 0.0, 1.0))


def load_signature(rows: Iterable[dict[str, Any]]) -> dict[str, float]:
    """Return signed gene effects, validating the compact CSV/JSON schema."""
    signature = {}
    for row in rows:
        gene, effect = row["gene"].strip(), float(row["logFC"])
        if not gene:
            raise ValueError("Disease signature contains an empty gene")
        signature[gene] = effect
    if not signature:
        raise ValueError("Disease signature is empty")
    return signature


def single_agent_metrics(drug: dict[str, Any], signature: dict[str, float]) -> dict[str, float]:
    """Measure signed reversal and target coverage for one candidate.

    ``target_effects`` values are signed [-1, 1] activity changes.  A positive
    disease logFC is improved by a negative drug effect, and conversely.
    """
    effects = drug["target_effects"]
    weights = {gene: abs(value) for gene, value in signature.items()}
    total = sum(weights.values())
    overlap = [gene for gene in effects if gene in signature]
    if not overlap or total == 0:
        return {"reversal": 0.0, "gene_coverage": 0.0, "pathway_coverage": 0.0}

    reversal = sum(
        weights[gene] * max(0.0, -np.sign(signature[gene]) * effects[gene])
        for gene in overlap
    ) / total
    gene_coverage = sum(weights[gene] for gene in overlap) / total
    affected_pathways = set(drug.get("pathways", []))
    pathway_coverage = len(affected_pathways & MS_DISEASE_PATHWAYS) / len(MS_DISEASE_PATHWAYS)
    return {
        "reversal": _clamp(reversal),
        "gene_coverage": _clamp(gene_coverage),
        "pathway_coverage": _clamp(pathway_coverage),
    }


def pair_metrics(
    drug_a: dict[str, Any], drug_b: dict[str, Any], signature: dict[str, float], config: PrioritizationConfig
) -> dict[str, Any]:
    """Score a pair using a conservative Bliss-like target effect model.

    The pair score is not a biological synergy estimate.  It is a transparent
    *prioritisation score* that rewards nonredundant signed reversal while
    penalising safety-class conflicts and low-evidence candidates.
    """
    a, b = single_agent_metrics(drug_a, signature), single_agent_metrics(drug_b, signature)
    genes = set(drug_a["target_effects"]) | set(drug_b["target_effects"])
    combined = dict(drug_a["target_effects"])
    for gene, effect in drug_b["target_effects"].items():
        # Saturating signed perturbation prevents unbounded sequential effects.
        combined[gene] = float(np.clip(combined.get(gene, 0.0) + effect, -1.0, 1.0))
    combined_drug = {"target_effects": combined, "pathways": list(set(drug_a["pathways"]) | set(drug_b["pathways"]))}
    combo = single_agent_metrics(combined_drug, signature)
    overlap = set(drug_a["target_effects"]) & set(drug_b["target_effects"])
    complementarity = 1.0 - len(overlap) / max(1, len(genes))
    evidence = min(EVIDENCE_RANK[drug_a["evidence_tier"]], EVIDENCE_RANK[drug_b["evidence_tier"]]) / 3.0
    uncertainty = (float(drug_a["target_uncertainty"]) + float(drug_b["target_uncertainty"])) / 2.0
    same_safety_class = bool(set(drug_a["safety_classes"]) & set(drug_b["safety_classes"]))
    same_mechanism = drug_a["mechanism_class"] == drug_b["mechanism_class"]
    exclusion = same_safety_class or same_mechanism
    priority = (
        combo["reversal"]
        + config.pathway_weight * combo["pathway_coverage"]
        + config.complementarity_weight * complementarity
        + config.evidence_weight * evidence
        - config.uncertainty_penalty * uncertainty
    )
    return {
        "drug_a": drug_a["name"], "drug_b": drug_b["name"], "priority_score": round(float(priority), 4),
        "signed_reversal": round(combo["reversal"], 4), "gene_coverage": round(combo["gene_coverage"], 4),
        "pathway_coverage": round(combo["pathway_coverage"], 4), "target_complementarity": round(complementarity, 4),
        "evidence_score": round(evidence, 4), "uncertainty": round(uncertainty, 4),
        "same_mechanism": same_mechanism, "shared_safety_class": same_safety_class,
        "excluded_from_primary_ranking": exclusion,
    }


def rank_pairs(drugs: list[dict[str, Any]], signature: dict[str, float], config: PrioritizationConfig | None = None) -> list[dict[str, Any]]:
    config = config or PrioritizationConfig()
    minimum_rank = EVIDENCE_RANK[config.minimum_evidence]
    eligible = [d for d in drugs if EVIDENCE_RANK[d["evidence_tier"]] >= minimum_rank]
    rows = [pair_metrics(a, b, signature, config) for a, b in combinations(eligible, 2)]
    return sorted(rows, key=lambda row: (row["excluded_from_primary_ranking"], -row["priority_score"], row["drug_a"], row["drug_b"]))


def sensitivity_interval(drug_a: dict[str, Any], drug_b: dict[str, Any], signature: dict[str, float], draws: int = 500, seed: int = 7) -> tuple[float, float]:
    """Bootstrap target-effect uncertainty; deterministic under a stated seed."""
    rng = np.random.default_rng(seed)
    config = PrioritizationConfig()
    scores = []
    for _ in range(draws):
        perturbed = []
        for drug in (drug_a, drug_b):
            clone = dict(drug)
            sd = float(drug["target_uncertainty"])
            clone["target_effects"] = {g: float(np.clip(v + rng.normal(0, sd), -1, 1)) for g, v in drug["target_effects"].items()}
            perturbed.append(clone)
        scores.append(pair_metrics(*perturbed, signature, config)["priority_score"])
    low, high = np.quantile(scores, [0.025, 0.975])
    return round(float(low), 4), round(float(high), 4)

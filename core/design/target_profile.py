"""Deriving a Target Product Profile: what one molecule would have to do.

The combination screen answers "which two existing agents pair best".  This
module answers the prior question -- *what should a molecule do at all* --
and answers it as a ranked, signed, weighted list of proteins to engage plus
the physicochemical envelope the molecule has to live inside.

Each candidate target earns a priority from four independent quantities:

``leverage``
    How much of the disease signature this protein commands.  Half of it is
    direct (the protein is itself a weighted signature gene), half is
    topological: influence over signature genes through the interactome,
    decayed by shortest-path distance.  A hub two steps upstream of twenty
    disease genes outranks a terminal gene that is merely dysregulated.

``tractability``
    The per-target small-molecule prior from the druggability annotation.
    This is the term that separates design from wishful thinking: CD20 has
    enormous leverage in MS and a tractability of 0.05, so no arm of a
    designed molecule will ever be pointed at it.

``liability``
    Safety burden inherited from the agents already known to engage the
    target, weighted by the disease context's risk domains.  A target reached
    only by agents carrying heavy infection risk starts with that debt.

``direction``
    The sign the molecule must push.  Taken from the signature where the
    target is itself a signature gene, and otherwise inferred from the panel:
    the direction applied by agents that actually move the signature the
    right way.

The output is a specification, not a prediction.  It says which proteins are
worth engaging given curated directional hypotheses; whether any molecule can
engage them together, at the same exposure, without toxicity, is exactly what
the rest of the pipeline cannot tell you.
"""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from typing import Any, Iterable, Sequence

import networkx as nx
import numpy as np

from core.biology.signature import Signature, alignment_metrics
from core.models.disease import DiseaseContext

# Influence decays by this factor per interactome step. At three steps a
# target retains an eighth of its weight, which is where the contribution
# stops being distinguishable from network background.
NETWORK_DECAY = 0.5
MAX_NETWORK_DISTANCE = 3

# Split between owning a signature gene outright and commanding it through
# the network. Declared rather than tuned: an even-ish split keeps a
# well-connected hub competitive with a heavily dysregulated terminal gene
# without letting topology dominate the transcriptional evidence.
DIRECT_LEVERAGE_SHARE = 0.6


@dataclass(frozen=True)
class TargetRequirement:
    """One protein the designed molecule should engage, and why."""

    gene: str
    target_class: str
    desired_direction: int
    priority: float
    leverage: float
    direct_leverage: float
    network_leverage: float
    tractability: float
    safety_liability: float
    therapeutic_axes: tuple[str, ...]
    pathway: str
    engaged_by: tuple[str, ...]
    rationale: str

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["therapeutic_axes"] = list(self.therapeutic_axes)
        payload["engaged_by"] = list(self.engaged_by)
        return payload


@dataclass(frozen=True)
class PropertyWindow:
    """The physicochemical envelope a candidate has to land inside."""

    molecular_weight: tuple[float, float]
    clogp: tuple[float, float]
    tpsa: tuple[float, float]
    hbd: tuple[float, float]
    rotatable_bonds: tuple[float, float]
    cns_mpo_floor: float
    rationale: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "molecular_weight": list(self.molecular_weight),
            "clogp": list(self.clogp),
            "tpsa": list(self.tpsa),
            "hbd": list(self.hbd),
            "rotatable_bonds": list(self.rotatable_bonds),
            "cns_mpo_floor": self.cns_mpo_floor,
            "rationale": self.rationale,
        }


@dataclass(frozen=True)
class TargetProductProfile:
    """The full specification handed to the generator."""

    disease: str
    requirements: tuple[TargetRequirement, ...]
    readouts: tuple[str, ...]
    property_window: PropertyWindow
    axis_gaps: dict[str, float]
    method: dict[str, Any]

    @property
    def genes(self) -> list[str]:
        return [r.gene for r in self.requirements]

    def direction(self, gene: str) -> int:
        for requirement in self.requirements:
            if requirement.gene == gene:
                return requirement.desired_direction
        return 0

    def weight(self, gene: str) -> float:
        for requirement in self.requirements:
            if requirement.gene == gene:
                return requirement.priority
        return 0.0

    def as_dict(self) -> dict[str, Any]:
        return {
            "disease": self.disease,
            "requirements": [r.as_dict() for r in self.requirements],
            "readouts_not_targeted": list(self.readouts),
            "property_window": self.property_window.as_dict(),
            "therapeutic_axis_gaps": self.axis_gaps,
            "method": self.method,
        }


# ---------------------------------------------------------------------------
# Component quantities
# ---------------------------------------------------------------------------

def network_leverage(
    graph: nx.Graph,
    signature: Signature,
    genes: Iterable[str],
    decay: float = NETWORK_DECAY,
    max_distance: int = MAX_NETWORK_DISTANCE,
) -> dict[str, float]:
    """Distance-decayed influence of each gene over the disease signature.

    A single-source shortest-path search per candidate, cut off at
    ``max_distance``; the score is the signature weight reachable within that
    horizon, discounted geometrically by distance.
    """
    scores: dict[str, float] = {}
    signature_weights = signature.weight
    for gene in genes:
        if gene not in graph:
            scores[gene] = 0.0
            continue
        total = 0.0
        distances = nx.single_source_shortest_path_length(graph, gene, cutoff=max_distance)
        for other, distance in distances.items():
            weight = signature_weights.get(other)
            if weight is None or distance == 0:
                continue
            total += weight * (decay ** distance)
        scores[gene] = total
    return scores


def infer_direction(gene: str, signature: Signature, panel: Sequence[dict]) -> tuple[int, str]:
    """The sign a molecule should push this target, and where that came from.

    Signature membership settles it outright.  Otherwise the panel votes:
    each agent's applied effect on the target is weighted by how well that
    agent moves the signature overall, so agents that get the biology right
    carry the vote.
    """
    if gene in signature.desired:
        return signature.desired[gene], "disease signature"
    vote = 0.0
    for drug in panel:
        effect = drug["target_effects"].get(gene)
        if not effect:
            continue
        alignment = alignment_metrics(drug["target_effects"], signature)
        vote += effect * alignment["reversal_efficiency"]
    if vote == 0.0:
        return 0, "undetermined"
    return (1 if vote > 0 else -1), "panel consensus weighted by signature alignment"


def target_safety_liability(gene: str, panel: Sequence[dict], weights: dict[str, float]) -> float:
    """Risk-weighted safety burden inherited from agents engaging this target.

    ``weights`` come from the disease registry entry, so a disease whose
    population tolerates infection risk poorly weights it accordingly.
    """
    weight_sum = sum(weights.values())
    numerator, denominator = 0.0, 0.0
    for drug in panel:
        effect = drug["target_effects"].get(gene)
        if not effect:
            continue
        burden = drug.get("safety_burden", {})
        weighted = sum(weights.get(d, 1.0) * float(burden.get(d, 0.0)) for d in weights) / weight_sum
        numerator += abs(effect) * weighted
        denominator += abs(effect)
    return round(numerator / denominator, 4) if denominator else 0.0


def target_axes(gene: str, panel: Sequence[dict]) -> tuple[str, ...]:
    axes: set[str] = set()
    for drug in panel:
        if gene in drug["target_effects"]:
            axes.update(drug.get("therapeutic_axes", []))
    return tuple(sorted(axes))


def axis_coverage_gaps(disease: DiseaseContext, panel: Sequence[dict]) -> dict[str, float]:
    """How poorly each therapeutic axis is served by *approved* agents.

    This is the argument for designing anything at all.  An axis covered by
    six approved drugs does not need a new molecule; an axis with none is
    where a new molecule has somewhere to go.  Scored as an unmet fraction in
    [0,1] so it can weight the profile directly.
    """
    approved = [d for d in panel if d.get("evidence_tier") == "approved"]
    counts = {axis: 0 for axis in disease.therapeutic_axes}
    for drug in approved:
        for axis in drug.get("therapeutic_axes", []):
            if axis in counts:
                counts[axis] += 1
    highest = max(counts.values()) or 1
    return {axis: round(1.0 - count / highest, 4) for axis, count in counts.items()}


# ---------------------------------------------------------------------------
# Profile construction
# ---------------------------------------------------------------------------

def _property_window(disease: DiseaseContext) -> PropertyWindow:
    """The envelope, tightened when the disease demands CNS exposure.

    The CNS numbers are the conventional central-nervous-system drug space
    (Wager's MPO desirability windows and Pajouhesh's guidelines): smaller,
    less polar, fewer donors than general oral space, because those are the
    properties that correlate with brain penetration.
    """
    if disease.delivery.requires_cns_exposure:
        return PropertyWindow(
            molecular_weight=(220.0, 420.0),
            clogp=(1.0, 3.5),
            tpsa=(40.0, 80.0),
            hbd=(0.0, 2.0),
            rotatable_bonds=(0.0, 7.0),
            cns_mpo_floor=disease.delivery.cns_mpo_floor,
            rationale=(
                "CNS-restricted envelope. "
                + disease.delivery.sanctuary_rationale
            ),
        )
    return PropertyWindow(
        molecular_weight=(250.0, 500.0),
        clogp=(0.5, 5.0),
        tpsa=(40.0, 130.0),
        hbd=(0.0, 4.0),
        rotatable_bonds=(0.0, 10.0),
        cns_mpo_floor=0.0,
        rationale="General oral small-molecule envelope; no sanctuary compartment to reach.",
    )


def build_target_profile(
    disease: DiseaseContext,
    top_n: int = 12,
    tractability_floor: float = 0.2,
    axis_gap_weight: float = 0.35,
) -> TargetProductProfile:
    """Derive the Target Product Profile for a disease context.

    ``top_n`` requirements are returned, ranked by priority.  Targets below
    ``tractability_floor`` are reported separately as readouts: the profile
    still wants their expression to move, but no designed arm will bind them.
    """
    signature = disease.signature()
    graph = disease.network()
    panel = disease.panel()
    druggability = disease.druggability()

    candidates = sorted({g for drug in panel for g in drug["target_effects"]} | set(signature.genes))
    candidates = [g for g in candidates if g in druggability]

    direct = {
        gene: signature.weight.get(gene, 0.0) for gene in candidates
    }
    networked = network_leverage(graph, signature, candidates)
    direct_max = max(direct.values()) or 1.0
    network_max = max(networked.values()) or 1.0

    gaps = axis_coverage_gaps(disease, panel)

    requirements: list[TargetRequirement] = []
    readouts: list[str] = []
    for gene in candidates:
        entry = druggability[gene]
        tractability = float(entry["small_molecule_tractability"])
        direction, direction_source = infer_direction(gene, signature, panel)
        direct_norm = direct[gene] / direct_max
        network_norm = networked[gene] / network_max
        leverage = (
            DIRECT_LEVERAGE_SHARE * direct_norm
            + (1.0 - DIRECT_LEVERAGE_SHARE) * network_norm
        )
        if tractability < tractability_floor or direction == 0:
            if leverage > 0:
                readouts.append(gene)
            continue

        liability = target_safety_liability(gene, panel, disease.risk_weights)
        axes = target_axes(gene, panel)
        # An axis nobody has an approved drug for is worth more.
        gap_bonus = max((gaps.get(axis, 0.0) for axis in axes), default=0.0)
        priority = (
            leverage
            * tractability
            * (1.0 - 0.5 * liability)
            * (1.0 + axis_gap_weight * gap_bonus)
        )
        engaged_by = tuple(sorted(
            d["name"] for d in panel if gene in d["target_effects"]
        ))[:6]
        arrow = "increase" if direction > 0 else "decrease"
        requirements.append(TargetRequirement(
            gene=gene,
            target_class=entry["target_class"],
            desired_direction=direction,
            priority=round(priority, 5),
            leverage=round(leverage, 5),
            direct_leverage=round(direct_norm, 5),
            network_leverage=round(network_norm, 5),
            tractability=tractability,
            safety_liability=liability,
            therapeutic_axes=axes,
            pathway=signature.pathway.get(gene, ""),
            engaged_by=engaged_by,
            rationale=(
                f"{arrow} {gene} ({entry['target_class']}); direction from {direction_source}; "
                f"{entry['note']}"
            ),
        ))

    requirements.sort(key=lambda r: (-r.priority, r.gene))
    selected = tuple(requirements[:top_n])

    return TargetProductProfile(
        disease=disease.identifier,
        requirements=selected,
        readouts=tuple(sorted(readouts)),
        property_window=_property_window(disease),
        axis_gaps=gaps,
        method={
            "n_candidates_considered": len(candidates),
            "n_requirements": len(selected),
            "n_readouts": len(readouts),
            "tractability_floor": tractability_floor,
            "network_decay": NETWORK_DECAY,
            "max_network_distance": MAX_NETWORK_DISTANCE,
            "direct_leverage_share": DIRECT_LEVERAGE_SHARE,
            "axis_gap_weight": axis_gap_weight,
            "priority_formula": (
                "leverage * tractability * (1 - 0.5*safety_liability) * "
                "(1 + axis_gap_weight * unmet_axis_fraction)"
            ),
            "caveat": (
                "A specification derived from curated directional hypotheses and a cached "
                "interactome. It states what would be worth engaging, not that any molecule "
                "can engage it."
            ),
        },
    )

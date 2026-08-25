"""The disease model: one object holding everything a campaign needs to know.

Before this module the MS campaign's vocabulary -- its pathways, therapeutic
axes, and safety domains -- lived as module constants inside the MS scoring
code.  That is fine for one disease and wrong for two: adding Parkinson's or
lupus would mean editing scoring logic rather than adding data.

A :class:`DiseaseContext` is loaded from a registry entry under
``data/diseases/`` and carries:

* the pointers to that disease's signature, interactome, and candidate panel;
* its controlled vocabularies (pathways, therapeutic axes, risk domains);
* its delivery constraint -- whether a molecule must cross into a sanctuary
  compartment to work, which for MS is the whole argument for the progressive
  phase and for a peripheral autoimmune disease is irrelevant;
* the known chemical matter to measure novelty against.

Everything downstream -- target profiling, fragment scoring, the Hamiltonian,
the generator -- reads the context and never names a disease.  Adding a
disease is therefore a data task.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[2]
REGISTRY_DIR = ROOT / "data/diseases"


@dataclass(frozen=True)
class DeliveryConstraint:
    """Where a molecule has to reach, and what that costs it.

    ``requires_cns_exposure`` is the flag that turns the CNS multi-parameter
    score from a nice-to-have into a gate.  ``sanctuary_rationale`` records
    *why*, because a constraint without a stated reason gets relaxed by the
    next person who finds it inconvenient.
    """

    requires_cns_exposure: bool = False
    cns_mpo_floor: float = 4.0
    preferred_routes: tuple[str, ...] = ("oral",)
    sanctuary_rationale: str = ""

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "DeliveryConstraint":
        return cls(
            requires_cns_exposure=bool(payload.get("requires_cns_exposure", False)),
            cns_mpo_floor=float(payload.get("cns_mpo_floor", 4.0)),
            preferred_routes=tuple(payload.get("preferred_routes", ("oral",))),
            sanctuary_rationale=str(payload.get("sanctuary_rationale", "")),
        )


@dataclass(frozen=True)
class DiseaseContext:
    """A disease as the engine sees it: vocabularies, data, and constraints."""

    identifier: str
    name: str
    description: str
    pathways: tuple[str, ...]
    therapeutic_axes: tuple[str, ...]
    risk_domains: tuple[str, ...]
    risk_weights: dict[str, float]
    gene_aliases: dict[str, str]
    delivery: DeliveryConstraint
    signature_path: Path
    network_path: Path
    panel_path: Path
    druggability_path: Path
    structures_path: Path | None = None
    provenance: dict[str, Any] = field(default_factory=dict)

    # -- data access ----------------------------------------------------
    def signature(self):
        """The signed disease signature (see ``core.biology.signature``)."""
        from core.biology.signature import load_signature

        return load_signature(self.signature_path)

    def network(self, min_score: float = 0.4):
        """The cached interactome, with node symbols mapped onto the signature's.

        Gene symbols drift, and interaction databases lag HGNC. STRING v12
        still calls glucocerebrosidase ``GBA`` where the current symbol is
        ``GBA1`` -- and in Parkinson's that is the single most common genetic
        risk factor, so leaving it unmapped silently zeroes the network
        leverage of the most important target in the disease. Aliases are
        declared per disease in the registry rather than patched into the
        cached data, so the cache stays a faithful record of what the source
        returned.
        """
        import networkx as nx

        from core.biology.network_proximity import load_network

        graph = load_network(self.network_path, min_score=min_score)
        if self.gene_aliases:
            applicable = {k: v for k, v in self.gene_aliases.items() if k in graph}
            if applicable:
                graph = nx.relabel_nodes(graph, applicable, copy=True)
        return graph

    def panel(self) -> list[dict[str, Any]]:
        payload = json.loads(Path(self.panel_path).read_text())
        return payload["drugs"] if isinstance(payload, dict) else payload

    def panel_metadata(self) -> dict[str, Any]:
        payload = json.loads(Path(self.panel_path).read_text())
        return payload.get("metadata", {}) if isinstance(payload, dict) else {}

    def druggability(self) -> dict[str, dict[str, Any]]:
        return json.loads(Path(self.druggability_path).read_text())["targets"]

    def known_structures(self) -> list[tuple[str, str]]:
        """``(name, SMILES)`` for every known agent that has a structure.

        Entries recorded without a structure -- biologics, and a few complex
        natural products -- are dropped here, which means novelty is measured
        against small-molecule chemical matter only.
        """
        if self.structures_path is None:
            return []
        payload = json.loads(Path(self.structures_path).read_text())
        return [
            (entry["name"], entry["smiles"])
            for entry in payload["structures"]
            if entry.get("smiles")
        ]

    # -- derived views --------------------------------------------------
    def tractable_targets(self, floor: float = 0.2) -> dict[str, float]:
        """Targets a small molecule may be assigned to bind, and their priors.

        Below ``floor`` a protein is a transcriptional readout: the profile
        can still require its expression to move, but no arm of a designed
        molecule will be pointed at it.
        """
        return {
            gene: entry["small_molecule_tractability"]
            for gene, entry in self.druggability().items()
            if entry["small_molecule_tractability"] >= floor
        }

    def target_classes(self) -> dict[str, str]:
        return {gene: entry["target_class"] for gene, entry in self.druggability().items()}

    def as_dict(self) -> dict[str, Any]:
        return {
            "identifier": self.identifier,
            "name": self.name,
            "description": self.description,
            "pathways": list(self.pathways),
            "therapeutic_axes": list(self.therapeutic_axes),
            "risk_domains": list(self.risk_domains),
            "risk_weights": self.risk_weights,
            "gene_aliases": self.gene_aliases,
            "delivery": {
                "requires_cns_exposure": self.delivery.requires_cns_exposure,
                "cns_mpo_floor": self.delivery.cns_mpo_floor,
                "preferred_routes": list(self.delivery.preferred_routes),
                "sanctuary_rationale": self.delivery.sanctuary_rationale,
            },
            "provenance": self.provenance,
        }


def _risk_weights(vocab: dict[str, Any]) -> dict[str, float]:
    """Per-domain risk weighting, defaulting to equal weight.

    How much each safety domain constrains use is a property of the disease
    and its patient population, not of a scoring function: infection risk
    dominates in a chronically immunosuppressed population and matters far
    less elsewhere. It therefore belongs in the registry entry, where it can
    be argued with, rather than as a default inside a scoring config.
    """
    declared = vocab.get("risk_weights") or {}
    weights = {domain: float(declared.get(domain, 1.0)) for domain in vocab["risk_domains"]}
    unknown = set(declared) - set(weights)
    if unknown:
        raise ValueError(
            f"risk_weights names domains absent from risk_domains: {sorted(unknown)}"
        )
    return weights


def _resolve(path_value: str | None) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value)
    return path if path.is_absolute() else ROOT / path


@lru_cache(maxsize=None)
def load_disease(identifier: str, registry_dir: Path | None = None) -> DiseaseContext:
    """Load a disease context by identifier from the registry directory."""
    directory = registry_dir or REGISTRY_DIR
    path = Path(directory) / f"{identifier}.json"
    if not path.exists():
        available = ", ".join(available_diseases(directory)) or "none"
        raise FileNotFoundError(
            f"No disease registry entry {identifier!r} in {directory}. Available: {available}"
        )
    payload = json.loads(path.read_text())
    data = payload["data"]
    missing = [
        key for key in ("signature", "network", "panel", "druggability")
        if not _resolve(data.get(key)) or not _resolve(data[key]).exists()
    ]
    if missing:
        raise FileNotFoundError(
            f"Disease {identifier!r} references data files that do not exist: {missing}"
        )
    vocab = payload["vocabulary"]
    return DiseaseContext(
        identifier=payload["identifier"],
        name=payload["name"],
        description=payload.get("description", ""),
        pathways=tuple(vocab["pathways"]),
        therapeutic_axes=tuple(vocab["therapeutic_axes"]),
        risk_domains=tuple(vocab["risk_domains"]),
        risk_weights=_risk_weights(vocab),
        gene_aliases=dict(payload.get("gene_aliases", {})),
        delivery=DeliveryConstraint.from_dict(payload.get("delivery", {})),
        signature_path=_resolve(data["signature"]),
        network_path=_resolve(data["network"]),
        panel_path=_resolve(data["panel"]),
        druggability_path=_resolve(data["druggability"]),
        structures_path=_resolve(data.get("structures")),
        provenance=payload.get("provenance", {}),
    )


def available_diseases(registry_dir: Path | None = None) -> list[str]:
    directory = Path(registry_dir or REGISTRY_DIR)
    if not directory.exists():
        return []
    return sorted(p.stem for p in directory.glob("*.json"))

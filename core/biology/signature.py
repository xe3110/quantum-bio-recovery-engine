"""Disease-agnostic contracts: signatures, effect algebra, evidence ranking.

These pieces began inside the MS screen's scoring module, which was the right
place for them while there was one disease. It stopped being the right place
the moment a second campaign needed them: the design layer was importing
``core.biology.ms_scoring`` to score molecules for any disease, which made the
generality claim false in the import graph regardless of what the registry
said.

Nothing here knows about multiple sclerosis, or about any disease. A signature
is a set of genes with a weight and a direction; an effect is a signed
fractional movement; alignment is how well one matches the other.
``core.biology.ms_scoring`` re-exports these names, so existing callers are
unaffected.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

# Ordinal rank over clinical-development stage. Disease-agnostic: the stages
# mean the same thing whatever is being developed.
EVIDENCE_RANK = {"approved": 3, "phase_3": 2, "phase_2": 1, "preclinical": 0}


@dataclass(frozen=True)
class Signature:
    """A signed disease signature with an explicit therapeutic direction.

    ``desired`` is the direction a therapy should move each gene, which is not
    the same as the direction the disease moved it. Compensatory transcripts
    rise as part of a protective response, and optimising to reverse them
    penalises agents that work. Carrying the therapeutic direction as its own
    field is what keeps that distinction from being lost.
    """

    logfc: dict[str, float]
    desired: dict[str, int]
    weight: dict[str, float]
    pathway: dict[str, str]

    @property
    def genes(self) -> list[str]:
        return list(self.logfc)

    @property
    def total_weight(self) -> float:
        return sum(self.weight.values())


def load_signature(path: Path | str) -> Signature:
    """Read a signature CSV, weighting genes by ``|logFC| x confidence``.

    Expected columns: ``gene``, ``logFC``, ``desired_direction`` (-1 or +1),
    and optionally ``pathway`` and ``confidence``. Duplicate or empty gene
    symbols and out-of-range directions are rejected here rather than being
    discovered as a strange score later.
    """
    logfc, desired, weight, pathway = {}, {}, {}, {}
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle):
            gene = row["gene"].strip()
            if not gene:
                raise ValueError("Disease signature contains an empty gene symbol")
            if gene in logfc:
                raise ValueError(f"Duplicate gene in signature: {gene}")
            direction = int(row["desired_direction"])
            if direction not in (-1, 1):
                raise ValueError(f"desired_direction must be -1 or 1 for {gene}")
            logfc[gene] = float(row["logFC"])
            desired[gene] = direction
            weight[gene] = abs(float(row["logFC"])) * float(row.get("confidence", 1.0))
            pathway[gene] = row.get("pathway", "")
    if not logfc:
        raise ValueError("Disease signature is empty")
    return Signature(logfc, desired, weight, pathway)


def bliss_combine(a: float, b: float) -> float:
    """Combine two signed fractional effects under Bliss independence.

    Same-signed effects saturate toward +/-1; opposing effects partially
    cancel, so mutual antagonism is represented rather than hidden by the
    discontinuity a hard clip would introduce.
    """
    if a == 0.0:
        return b
    if b == 0.0:
        return a
    if (a > 0) == (b > 0):
        sign = 1.0 if a > 0 else -1.0
        return float(sign * (abs(a) + abs(b) - abs(a) * abs(b)))
    return float(np.clip(a + b, -1.0, 1.0))


def combine_effects(effects_a: dict[str, float], effects_b: dict[str, float]) -> dict[str, float]:
    """Union of two effect vectors, combined gene-wise under Bliss."""
    combined = dict(effects_a)
    for gene, value in effects_b.items():
        combined[gene] = bliss_combine(combined.get(gene, 0.0), value)
    return combined


def alignment_metrics(effects: dict[str, float], signature: Signature) -> dict[str, float]:
    """Split an effect vector into therapeutic and counter-therapeutic movement.

    ``reversal`` is absolute and spans a narrow range, because any one agent
    engages a handful of a signature's genes. ``reversal_efficiency`` --
    therapeutic movement divided by the signal actually engaged -- spans [0,1]
    and is what makes agents of different breadth comparable.
    """
    total = signature.total_weight
    overlap = [g for g in effects if g in signature.logfc]
    if not overlap or total == 0:
        return {"reversal": 0.0, "counter_therapeutic": 0.0, "gene_coverage": 0.0,
                "reversal_efficiency": 0.0}

    therapeutic = 0.0
    counter = 0.0
    for gene in overlap:
        aligned = float(np.clip(signature.desired[gene] * effects[gene], -1.0, 1.0))
        w = signature.weight[gene]
        if aligned >= 0:
            therapeutic += w * aligned
        else:
            counter += w * -aligned
    coverage = sum(signature.weight[g] for g in overlap) / total
    return {
        "reversal": therapeutic / total,
        "counter_therapeutic": counter / total,
        "gene_coverage": coverage,
        "reversal_efficiency": (therapeutic / coverage / total) if coverage > 0 else 0.0,
    }

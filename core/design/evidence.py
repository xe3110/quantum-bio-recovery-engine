"""What has and has not been assessed about a designed molecule.

Every caveat in this project's documentation is true and none of it travels
with the data. A result file read by a script, pasted into a slide, or handed
to someone who never opened the protocol carries a ``fitness`` number that
looks like a prediction and is not one.

So the absence of evidence is emitted as **structured data on every
candidate**, not as prose elsewhere. A consumer that wants to filter on "has
any binding evidence" can do so; one that ignores the field is at least
ignoring something explicit.

The list is deliberately long and deliberately unflattering. It is the honest
gap between what this pipeline computes -- a specification match and a
physicochemical property vector -- and what would justify calling a structure
a drug candidate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# axis -> (what would satisfy it, why its absence matters here)
ASSESSMENT_AXES: dict[str, tuple[str, str]] = {
    "target_binding": (
        "docking, free-energy perturbation, or a measured Ki/IC50/EC50",
        "Target engagement is inherited from the parent chemotype of each fragment, which "
        "assumes a motif keeps its activity when transplanted onto a scaffold it did not "
        "evolve on. This is the assumption most likely to be wrong and the first that "
        "should be tested.",
    ),
    "selectivity": (
        "a kinome or receptor off-target panel",
        "Several library arms are deliberately promiscuous. A designed multi-target ligand "
        "and an uncontrolled polypharmacology liability are the same molecule seen from "
        "two sides, and nothing here distinguishes them.",
    ),
    "cell_activity": (
        "a target-engagement or phenotypic assay in a relevant cell type",
        "No compound has been made or tested. The predicted transcriptional effects are "
        "curated directional hypotheses, not measurements.",
    ),
    "in_vivo_efficacy": (
        "a disease model with a functional endpoint",
        "Neither disease in this registry has a model that reliably predicts clinical "
        "benefit, which is why both have long records of translational failure.",
    ),
    "permeability": (
        "PAMPA, Caco-2, or MDCK permeability",
        "Passive permeability is assumed from physicochemical properties alone.",
    ),
    "bbb_transport": (
        "in situ perfusion, brain:plasma ratio, or unbound partition coefficient",
        "The CNS multi-parameter score is a physicochemical prior, not a transport model. "
        "It says a molecule looks like CNS drugs; it does not say it reaches the brain.",
    ),
    "efflux_liability": (
        "a P-glycoprotein and BCRP substrate assay",
        "P-gp efflux decides many CNS programmes on its own and is not modelled anywhere "
        "in this pipeline. A molecule can satisfy every property window and still be "
        "pumped straight back out.",
    ),
    "metabolic_stability": (
        "microsomal or hepatocyte clearance",
        "Metabolic soft spots are noted qualitatively in fragment annotations and are not "
        "computed.",
    ),
    "cyp_liability": (
        "CYP inhibition and induction panels",
        "Drug-drug interaction risk is unassessed, which matters most for the combination "
        "use these molecules are designed to replace.",
    ),
    "cardiac_safety": (
        "hERG patch-clamp",
        "Basic, lipophilic amines are the classic hERG chemotype and the CNS envelope "
        "actively pushes designs toward them.",
    ),
    "solubility": (
        "kinetic and thermodynamic solubility",
        "Not computed. A design can satisfy every property window and be unformulable.",
    ),
    "synthetic_route": (
        "a proposed route and a made compound",
        "The tractability score is a declared size-, stereo-, and topology-based proxy, "
        "explicitly not a retrosynthetic analysis.",
    ),
    "intellectual_property": (
        "a freedom-to-operate and novelty search",
        "The Tanimoto novelty check compares against this project's own reference set of "
        "a few dozen structures. It is not a prior-art search and carries no legal weight.",
    ),
}


@dataclass(frozen=True)
class EvidenceStatus:
    """The assessment state of one candidate. Currently: none of it."""

    assessed: dict[str, Any]
    unassessed: list[str]

    def as_dict(self) -> dict[str, Any]:
        return {
            "assessed": self.assessed,
            "unassessed": self.unassessed,
            "n_unassessed": len(self.unassessed),
            "readiness": (
                "hypothesis_only" if not self.assessed else "partially_assessed"
            ),
            "statement": (
                "No experimental or physics-based evidence supports this structure. What "
                "was computed is a match to a target product profile and a physicochemical "
                "property vector, both derived from curated hypotheses. This is a starting "
                "point for a medicinal chemistry programme, not a drug candidate."
            ),
            "detail": {
                axis: {"would_require": requirement, "why_it_matters": reason}
                for axis, (requirement, reason) in ASSESSMENT_AXES.items()
                if axis in self.unassessed
            },
        }


def evidence_status(assessed: dict[str, Any] | None = None) -> EvidenceStatus:
    """Build the evidence block.

    ``assessed`` is the hook for a future assay-feedback loop: pass in whatever
    real evidence exists and those axes drop out of the unassessed list.
    Nothing in the pipeline populates it today, which is the point.
    """
    assessed = assessed or {}
    return EvidenceStatus(
        assessed=assessed,
        unassessed=[axis for axis in ASSESSMENT_AXES if axis not in assessed],
    )

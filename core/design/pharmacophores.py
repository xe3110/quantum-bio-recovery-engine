"""The fragment library: chemical matter annotated with what it is for.

A :class:`Fragment` pairs a substructure with the target-engagement
hypothesis its chemotype carries.  That pairing is the whole basis of
assembly-based design, and it is the weakest link in the pipeline: it assumes
a motif keeps its activity when transplanted onto a scaffold it did not
evolve on.  Real medicinal chemistry programmes spend most of their effort
discovering exactly where that assumption breaks.

Nothing here is a measurement.  The engagement values are curated directional
hypotheses on the same signed scale as the drug panel, which is what lets an
assembled molecule be scored against a disease signature at all.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Sequence

from core.chemistry.descriptors import Descriptors, compute_descriptors
from core.chemistry.molecule import Molecule, cap_attachments, parse_smiles

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LIBRARY = ROOT / "data/chemistry/pharmacophore_library.json"

ROLE_ATTACHMENTS = {"pharmacophore": 1, "cap": 1, "linker": 2}

# How much to believe a chemotype engages its primary targets. Ordinal, and
# deliberately coarse: the distance between an approved drug and a speculative
# chemotype is what matters, not a second decimal place.
EVIDENCE_CONFIDENCE = {
    "approved_drug": 0.90,
    "clinical_candidate": 0.75,
    "published_chemotype": 0.60,
    "speculative": 0.30,
}

# A downstream transcriptional consequence is a weaker claim than a binding
# claim from the same chemotype, so it is discounted by this factor. An
# approved drug's binding target lands at 0.90; the transcripts it is expected
# to move land at 0.45.
DOWNSTREAM_DISCOUNT = 0.5


@dataclass(frozen=True)
class Fragment:
    """One building block, its structure, and its engagement hypothesis."""

    identifier: str
    name: str
    role: str
    smiles: str
    engages: dict[str, float] = field(default_factory=dict)
    target_family: str = ""
    parent_chemotype: str = ""
    therapeutic_axes: tuple[str, ...] = ()
    liabilities: tuple[str, ...] = ()
    note: str = ""
    geometry: str = ""
    length: int | None = None
    evidence_tier: str = "speculative"
    primary_targets: tuple[str, ...] = ()

    def claim_confidence(self, gene: str) -> float:
        """How much to believe this fragment engages ``gene``.

        Two things set it: the strength of published precedent for the
        chemotype, and whether the gene is a target it is claimed to *bind* or
        a transcript it is expected to move downstream. Without this
        distinction a scaffold assembled from speculative downstream claims
        scores identically to one built on approved-drug pharmacology.
        """
        if gene not in self.engages:
            return 0.0
        base = EVIDENCE_CONFIDENCE.get(self.evidence_tier, 0.3)
        return round(base if gene in self.primary_targets else base * DOWNSTREAM_DISCOUNT, 4)

    @property
    def mean_confidence(self) -> float:
        if not self.engages:
            return 0.0
        return round(
            sum(self.claim_confidence(g) for g in self.engages) / len(self.engages), 4
        )

    def molecule(self) -> Molecule:
        return parse_smiles(self.smiles, name=self.name)

    @property
    def attachment_count(self) -> int:
        return len(self.molecule().attachment_points)

    @property
    def heavy_atoms(self) -> int:
        """Heavy-atom count excluding the attachment markers.

        This is the fragment's cost against the molecule's size budget, which
        for a CNS-restricted profile is the binding constraint on how many
        arms a design can afford.
        """
        return len(self.molecule().heavy_atoms)

    def capped_descriptors(self) -> Descriptors:
        """Descriptors for the fragment with its attachment points hydrogenated.

        Polarity and donor count are close to additive across an assembly, so
        a fragment's capped TPSA and HBD are a usable estimate of what it will
        contribute to the finished molecule.  Computing them here is what lets
        the delivery envelope enter fragment selection rather than only
        judging molecules after they are built.
        """
        return compute_descriptors(cap_attachments(self.molecule()))

    @property
    def targets(self) -> set[str]:
        return set(self.engages)

    def as_dict(self) -> dict[str, Any]:
        return {
            "id": self.identifier,
            "name": self.name,
            "role": self.role,
            "smiles": self.smiles,
            "target_family": self.target_family,
            "parent_chemotype": self.parent_chemotype,
            "therapeutic_axes": list(self.therapeutic_axes),
            "engages": self.engages,
            "evidence_tier": self.evidence_tier,
            "primary_targets": list(self.primary_targets),
            "claim_confidence": {g: self.claim_confidence(g) for g in sorted(self.engages)},
            "liabilities": list(self.liabilities),
            "heavy_atoms": self.heavy_atoms,
            "note": self.note,
        }


@dataclass(frozen=True)
class FragmentLibrary:
    """Every fragment, indexed by role and identifier."""

    fragments: tuple[Fragment, ...]
    metadata: dict[str, Any] = field(default_factory=dict)

    def by_role(self, role: str) -> list[Fragment]:
        return [f for f in self.fragments if f.role == role]

    def get(self, identifier: str) -> Fragment:
        for fragment in self.fragments:
            if fragment.identifier == identifier:
                return fragment
        raise KeyError(f"No fragment {identifier!r} in the library")

    @property
    def pharmacophores(self) -> list[Fragment]:
        return self.by_role("pharmacophore")

    @property
    def scaffolds(self) -> list[Fragment]:
        return self.by_role("scaffold")

    @property
    def linkers(self) -> list[Fragment]:
        return self.by_role("linker")

    @property
    def caps(self) -> list[Fragment]:
        return self.by_role("cap")

    def covering(self, gene: str) -> list[Fragment]:
        """Every pharmacophore claiming to move ``gene``."""
        return [f for f in self.pharmacophores if gene in f.engages]

    def parent_structures(self) -> list[tuple[str, str]]:
        """``(name, SMILES)`` for each fragment, for novelty comparison.

        The attachment marker is left in place: comparing against the bare
        fragment is the check for whether a design is merely its own
        starting material.
        """
        return [(f.name, f.smiles.replace("*", "C")) for f in self.fragments
                if f.role == "pharmacophore"]

    def __len__(self) -> int:
        return len(self.fragments)


def _validate(fragment: Fragment) -> None:
    expected = ROLE_ATTACHMENTS.get(fragment.role)
    actual = fragment.attachment_count
    if expected is not None and actual != expected:
        raise ValueError(
            f"Fragment {fragment.identifier!r} has role {fragment.role!r} and so needs "
            f"{expected} attachment point(s), but its SMILES carries {actual}"
        )
    if fragment.role == "scaffold" and actual < 2:
        raise ValueError(
            f"Scaffold {fragment.identifier!r} needs at least two attachment points, has {actual}"
        )
    for gene, value in fragment.engages.items():
        if not -1.0 <= float(value) <= 1.0:
            raise ValueError(
                f"Fragment {fragment.identifier!r} claims an engagement of {value} for {gene}; "
                "the convention is a signed value in [-1, 1]"
            )
    if fragment.role == "pharmacophore":
        if fragment.evidence_tier not in EVIDENCE_CONFIDENCE:
            raise ValueError(
                f"Fragment {fragment.identifier!r} declares evidence_tier "
                f"{fragment.evidence_tier!r}; expected one of {sorted(EVIDENCE_CONFIDENCE)}"
            )
        if not fragment.primary_targets:
            raise ValueError(
                f"Pharmacophore {fragment.identifier!r} declares no primary_targets. Every "
                "engagement claim would be scored as downstream inference."
            )
        stray = [g for g in fragment.primary_targets if g not in fragment.engages]
        if stray:
            raise ValueError(
                f"Fragment {fragment.identifier!r} lists primary targets absent from its "
                f"engagement map: {stray}"
            )


@lru_cache(maxsize=None)
def load_library(path: Path | str = DEFAULT_LIBRARY) -> FragmentLibrary:
    """Load and validate the fragment library.

    Validation is strict and happens at load: a fragment whose SMILES does not
    carry the attachment points its role requires, or whose engagement values
    fall outside the signed convention, stops the campaign here rather than
    producing a molecule that cannot be assembled or scored.
    """
    payload = json.loads(Path(path).read_text())
    fragments = []
    seen: set[str] = set()
    for entry in payload["fragments"]:
        identifier = entry["id"]
        if identifier in seen:
            raise ValueError(f"Duplicate fragment id {identifier!r} in {path}")
        seen.add(identifier)
        fragment = Fragment(
            identifier=identifier,
            name=entry["name"],
            role=entry["role"],
            smiles=entry["smiles"],
            engages={k: float(v) for k, v in entry.get("engages", {}).items()},
            target_family=entry.get("target_family", ""),
            parent_chemotype=entry.get("parent_chemotype", ""),
            therapeutic_axes=tuple(entry.get("therapeutic_axes", ())),
            liabilities=tuple(entry.get("liabilities", ())),
            note=entry.get("note", ""),
            geometry=entry.get("geometry", ""),
            length=entry.get("length"),
            evidence_tier=entry.get("evidence_tier", "speculative"),
            primary_targets=tuple(entry.get("primary_targets", ())),
        )
        _validate(fragment)
        fragments.append(fragment)
    return FragmentLibrary(tuple(fragments), payload.get("metadata", {}))


def unreachable_requirements(library: FragmentLibrary, genes: Iterable[str]) -> list[str]:
    """Profile targets no fragment in the library claims to engage.

    Reported rather than silently dropped: a high-priority requirement with no
    chemical matter behind it is the most useful thing a design run can tell
    you, because it names what the library is missing.
    """
    covered = {gene for fragment in library.pharmacophores for gene in fragment.engages}
    return [gene for gene in genes if gene not in covered]

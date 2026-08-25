"""Circular fingerprints and similarity, for novelty and redundancy checks.

A generator that rediscovers dimethyl fumarate has not designed anything.  The
only way to say so is to compare every proposal against the known chemical
matter it was assembled from, which needs a structural fingerprint.

:func:`circular_fingerprint` follows the Morgan/ECFP construction (Rogers &
Hahn 2010): each atom is given an initial invariant, then that invariant is
repeatedly re-hashed together with its neighbours' invariants, once per
radius step.  The set of every identifier seen at every step is the
fingerprint.  Radius 2 corresponds to the commonly used ECFP4.

Two deliberate differences from a production toolkit: identifiers are hashed
into a large space and kept as a set rather than folded into a fixed-length
bit vector, so there is no folding collision to reason about; and
stereochemistry is absent, because the molecule model does not carry it.  A
Tanimoto computed here is therefore comparable *within* this project and not
with a published ECFP4 value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from core.chemistry.descriptors import ring_atoms
from core.chemistry.molecule import Molecule, parse_smiles

_HASH_SPACE = 2 ** 32


def _stable_hash(payload: tuple) -> int:
    """A hash that does not change between interpreter runs.

    ``hash()`` is salted per process, which would make every fingerprint --
    and so every novelty claim -- irreproducible.
    """
    value = 2166136261  # FNV-1a offset basis
    for item in payload:
        for byte in str(item).encode("utf-8"):
            value ^= byte
            value = (value * 16777619) % _HASH_SPACE
        value = (value * 16777619) % _HASH_SPACE
    return value


def _initial_invariants(mol: Molecule) -> dict[int, int]:
    """Daylight-style atom invariants: the ECFP starting identifiers."""
    in_ring = ring_atoms(mol)
    invariants = {}
    for atom in mol.atoms:
        invariants[atom.index] = _stable_hash((
            atom.element,
            mol.degree(atom.index),
            mol.implicit_hydrogens(atom.index),
            atom.charge,
            int(atom.aromatic),
            int(atom.index in in_ring),
        ))
    return invariants


def circular_fingerprint(molecule: Molecule | str, radius: int = 2) -> set[int]:
    """Return the set of circular-substructure identifiers up to ``radius``."""
    mol = parse_smiles(molecule) if isinstance(molecule, str) else molecule
    invariants = _initial_invariants(mol)
    features = set(invariants.values())
    for _ in range(radius):
        updated = {}
        for atom in mol.atoms:
            environment = sorted(
                (bond.order, int(bond.aromatic), invariants[bond.other(atom.index)])
                for bond in mol.bonds_of(atom.index)
            )
            updated[atom.index] = _stable_hash((invariants[atom.index], *environment))
        invariants = updated
        features.update(invariants.values())
    return features


def tanimoto(a: set[int], b: set[int]) -> float:
    """Jaccard/Tanimoto coefficient between two fingerprint sets."""
    if not a and not b:
        return 1.0
    union = len(a | b)
    return round(len(a & b) / union, 4) if union else 0.0


def similarity(first: Molecule | str, second: Molecule | str, radius: int = 2) -> float:
    return tanimoto(circular_fingerprint(first, radius), circular_fingerprint(second, radius))


@dataclass(frozen=True)
class NoveltyReport:
    """How structurally close a proposal sits to the nearest known compound."""

    nearest_name: str
    nearest_smiles: str
    max_similarity: float
    is_novel: bool
    threshold: float

    def as_dict(self) -> dict:
        return {
            "nearest_known_compound": self.nearest_name,
            "nearest_known_smiles": self.nearest_smiles,
            "max_tanimoto": self.max_similarity,
            "is_novel": self.is_novel,
            "novelty_threshold": self.threshold,
        }


def assess_novelty(
    candidate: Molecule | str,
    reference_library: Sequence[tuple[str, str]],
    threshold: float = 0.6,
    radius: int = 2,
) -> NoveltyReport:
    """Compare a candidate against ``(name, smiles)`` pairs of known compounds.

    ``threshold`` is the Tanimoto above which the candidate is treated as a
    close analogue rather than new chemical matter.  0.6 on ECFP4 is a
    conventional line for "same series"; it is a convention, not a legal or
    patentability standard, and nothing here constitutes a novelty search.
    """
    query = circular_fingerprint(candidate, radius)
    best_name, best_smiles, best_score = "", "", 0.0
    for name, smiles in reference_library:
        try:
            score = tanimoto(query, circular_fingerprint(smiles, radius))
        except Exception:
            continue
        if score > best_score:
            best_name, best_smiles, best_score = name, smiles, score
    return NoveltyReport(
        nearest_name=best_name,
        nearest_smiles=best_smiles,
        max_similarity=best_score,
        is_novel=best_score < threshold,
        threshold=threshold,
    )


def diversity_select(
    candidates: Sequence[tuple[str, str]],
    k: int,
    radius: int = 2,
    min_distance: float = 0.35,
) -> list[int]:
    """Greedy MaxMin selection of ``k`` structurally distinct candidates.

    Returns indices into ``candidates``, preserving input order as the
    tie-break so the selection is deterministic.  A design campaign that
    advances five near-identical molecules has tested one hypothesis, not
    five; this is what stops that happening.
    """
    if not candidates:
        return []
    fingerprints = [circular_fingerprint(smiles, radius) for _name, smiles in candidates]
    chosen = [0]
    while len(chosen) < min(k, len(candidates)):
        best_index, best_distance = None, -1.0
        for index in range(len(candidates)):
            if index in chosen:
                continue
            distance = min(1.0 - tanimoto(fingerprints[index], fingerprints[c]) for c in chosen)
            if distance > best_distance:
                best_index, best_distance = index, distance
        if best_index is None or best_distance < min_distance:
            break
        chosen.append(best_index)
    return chosen

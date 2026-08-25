"""Physicochemical descriptors computed directly on the molecular graph.

Every descriptor here is an *estimator*, and the estimator is named in its
docstring so a reader can check the original method rather than trusting a
number.  Two in particular are additive-contribution models:

``tpsa``
    Ertl, Rohde & Selzer (2000), *J. Med. Chem.* 43:3714 -- topological polar
    surface area as a sum of per-fragment contributions for nitrogen and
    oxygen environments, with the sulfur/phosphorus extension available but
    off by default, matching the convention most filters assume.

``clogp``
    Wildman & Crippen (1999), *J. Chem. Inf. Comput. Sci.* 39:868 -- an
    atom-typed additive logP.  This is a reduced implementation of the
    published atom types: the common environments in drug-like space are
    typed individually and the rare ones fall back to their class default, so
    expect agreement with a full implementation to roughly a few tenths of a
    log unit and a growing error for exotic chemotypes.

Both are stereochemistry-blind and conformation-blind by construction, and
neither is a substitute for a measurement.  When RDKit is importable the
molecule module's canonicaliser uses it, but these descriptors deliberately
stay self-contained so that scores are reproducible across environments.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import lru_cache
from typing import Any

import networkx as nx

from core.chemistry import backend as chemistry_backend
from core.chemistry.molecule import Molecule, parse_smiles

HALOGENS = {"F", "Cl", "Br", "I"}
HETEROATOMS = {"N", "O", "S", "P"} | HALOGENS
POLAR = {"N", "O"}


# ---------------------------------------------------------------------------
# Ring perception
# ---------------------------------------------------------------------------

def ring_systems(mol: Molecule) -> list[list[int]]:
    """Smallest set of smallest rings, via a minimum-weight cycle basis.

    ``networkx.cycle_basis`` is a *spanning-tree* basis and is not guaranteed
    to return the smallest cycles: on adamantane it reports two six-membered
    rings and one eight-membered ring, where three six-membered rings exist.
    That mattered beyond bookkeeping -- ``synthetic_tractability`` charges a
    macrocycle penalty for any ring larger than seven atoms, so every
    adamantane-containing design was being penalised for a macrocycle it does
    not have. ``minimum_cycle_basis`` returns the genuine smallest set.

    A residual difference from RDKit remains by construction: RDKit's
    ``CalcNumRings`` uses a *symmetrised* SSSR, which reports extra rings on
    symmetric bridged cages to preserve their symmetry, so it counts five
    rings in adamantane where the cycle rank is three. Neither is wrong; they
    answer different questions.
    """
    edges = tuple(sorted((min(b.a, b.b), max(b.a, b.b)) for b in mol.bonds))
    return [list(ring) for ring in _minimum_cycle_basis(edges, len(mol.atoms))]


@lru_cache(maxsize=8192)
def _minimum_cycle_basis(edges: tuple[tuple[int, int], ...], n_atoms: int) -> tuple[tuple[int, ...], ...]:
    """Cached ring perception. Keyed on the edge set, so it is safe to share
    across molecules and stable across runs.

    The cache matters: ring perception is called from the rotatable-bond
    count, the tractability proxy, and several alert rules, and a design
    campaign evaluates thousands of structures.
    """
    graph = nx.Graph()
    graph.add_nodes_from(range(n_atoms))
    graph.add_edges_from(edges)
    return tuple(tuple(sorted(cycle)) for cycle in nx.minimum_cycle_basis(graph))


def ring_atoms(mol: Molecule) -> set[int]:
    return {index for ring in ring_systems(mol) for index in ring}


def aromatic_rings(mol: Molecule) -> list[list[int]]:
    return [r for r in ring_systems(mol) if all(mol.atoms[i].aromatic for i in r)]


def is_in_ring(mol: Molecule, index: int) -> bool:
    return index in ring_atoms(mol)


# ---------------------------------------------------------------------------
# Hydrogen bonding
# ---------------------------------------------------------------------------

def hbd_count(mol: Molecule) -> int:
    """Lipinski donors: hydrogens attached to nitrogen or oxygen."""
    return sum(
        mol.implicit_hydrogens(a.index)
        for a in mol.atoms
        if a.element in POLAR
    )


def hba_count(mol: Molecule) -> int:
    """Lipinski acceptors: the count of nitrogen and oxygen atoms."""
    return sum(1 for a in mol.atoms if a.element in POLAR)


# ---------------------------------------------------------------------------
# Topological polar surface area (Ertl)
# ---------------------------------------------------------------------------

def _oxygen_tpsa(mol: Molecule, index: int) -> float:
    atom = mol.atoms[index]
    hydrogens = mol.implicit_hydrogens(index)
    if atom.charge < 0:
        return 23.06
    if atom.aromatic:
        return 13.14
    doubly_bonded = any(b.order == 2 for b in mol.bonds_of(index))
    if doubly_bonded:
        return 17.07
    if hydrogens == 1:
        return 20.23
    return 9.23  # ether / ester oxygen


def _nitrogen_tpsa(mol: Molecule, index: int) -> float:
    atom = mol.atoms[index]
    hydrogens = mol.implicit_hydrogens(index)
    bonds = mol.bonds_of(index)
    connections = len(bonds)
    max_order = max((b.order for b in bonds), default=1)

    if atom.aromatic:
        if atom.charge > 0:
            if any(
                mol.atoms[b.other(index)].element == "O" and mol.atoms[b.other(index)].charge < 0
                for b in bonds
            ):
                return 4.10  # aromatic N-oxide
            return 14.14 if hydrogens else 4.10
        if hydrogens:
            return 15.79
        return 4.41 if connections >= 3 else 12.89

    if atom.charge > 0:
        return {0: 0.00, 1: 4.44, 2: 16.61, 3: 27.64}.get(hydrogens, 0.00)
    if max_order == 3:
        return 23.79  # nitrile
    if max_order == 2:
        return 12.36 if hydrogens == 0 else 23.85
    return {0: 3.24, 1: 12.03, 2: 26.02}.get(hydrogens, 3.24)


_SULFUR_TPSA = {"thiol": 38.80, "sulfide": 25.30, "aromatic": 28.24,
                "sulfoxide": 19.21, "sulfone": 8.38}


def tpsa(mol: Molecule, include_s_p: bool = False) -> float:
    """Topological polar surface area in square angstroms (Ertl 2000).

    ``include_s_p`` adds the sulfur and phosphorus contributions from the same
    paper.  It defaults to off because the cut-offs quoted in the literature
    (90 A^2 for oral, 60-70 A^2 for CNS) were derived on the N/O-only form.
    """
    total = 0.0
    for atom in mol.atoms:
        if atom.element == "O":
            total += _oxygen_tpsa(mol, atom.index)
        elif atom.element == "N":
            total += _nitrogen_tpsa(mol, atom.index)
        elif include_s_p and atom.element == "S":
            oxides = sum(1 for b in mol.bonds_of(atom.index)
                         if b.order == 2 and mol.atoms[b.other(atom.index)].element == "O")
            if atom.aromatic:
                total += _SULFUR_TPSA["aromatic"]
            elif oxides >= 2:
                total += _SULFUR_TPSA["sulfone"]
            elif oxides == 1:
                total += _SULFUR_TPSA["sulfoxide"]
            elif mol.implicit_hydrogens(atom.index):
                total += _SULFUR_TPSA["thiol"]
            else:
                total += _SULFUR_TPSA["sulfide"]
        elif include_s_p and atom.element == "P":
            total += 13.59 if mol.implicit_hydrogens(atom.index) == 0 else 23.47
    return round(total, 2)


# ---------------------------------------------------------------------------
# Crippen logP
# ---------------------------------------------------------------------------

_H_CONTRIBUTION = {"C": 0.1230, "O": -0.2677, "N": 0.2142, "S": 0.2980, "P": 0.2980}


def _carbon_logp(mol: Molecule, index: int) -> float:
    atom = mol.atoms[index]
    bonds = mol.bonds_of(index)
    neighbours = [mol.atoms[b.other(index)] for b in bonds]
    hydrogens = mol.implicit_hydrogens(index)
    max_order = max((b.order for b in bonds), default=1)

    # A carbon double-bonded to a heteroatom is a carbonyl or imine carbon
    # whatever ring it sits in. This test has to come before the aromatic
    # branch: the exocyclic carbonyls of a purinone or pyrimidinone are written
    # as aromatic ring atoms, and typing them as aromatic ethers instead
    # over-estimated caffeine and uracil by roughly three log units.
    doubled = [mol.atoms[b.other(index)] for b in bonds if b.order == 2]
    if any(n.element in HETEROATOMS for n in doubled):
        return -0.2783

    if atom.aromatic:
        elements = {n.element for n in neighbours}
        if elements & HALOGENS:
            return -0.8186
        if "O" in elements:
            return 0.5437
        if "N" in elements:
            return 0.4619
        if "S" in elements:
            return 0.1893
        if hydrogens:
            return 0.1581
        if sum(1 for n in neighbours if n.aromatic) >= 3:
            return 0.2955  # ring-fusion carbon
        return 0.1360

    if max_order == 3:
        return 0.0017
    if max_order == 2:
        return 0.1551
    if any(n.element in HETEROATOMS for n in neighbours):
        return -0.2035 if hydrogens else -0.2051
    return 0.1441 if hydrogens else 0.0000


def _nitrogen_logp(mol: Molecule, index: int) -> float:
    atom = mol.atoms[index]
    hydrogens = mol.implicit_hydrogens(index)
    bonds = mol.bonds_of(index)
    max_order = max((b.order for b in bonds), default=1)

    if atom.charge > 0:
        return -1.9500
    if max_order == 3:
        return -0.3239  # nitrile

    # Amide typing precedes the aromatic branch for the same reason as carbon:
    # a ring nitrogen flanking a lactam carbonyl is amide-like, and the
    # heteroaromatic type is the wrong one for it.
    for bond in bonds:
        other = bond.other(index)
        if mol.atoms[other].element == "C" and any(
            b.order == 2 and mol.atoms[b.other(other)].element == "O" for b in mol.bonds_of(other)
        ):
            return -0.4045  # amide / lactam / carbamate nitrogen

    if atom.aromatic:
        return -0.3187
    if max_order == 2:
        return -0.4458
    return {2: -1.0190, 1: -0.7096, 0: -1.0270}.get(hydrogens, -1.0270)


def _oxygen_logp(mol: Molecule, index: int) -> float:
    atom = mol.atoms[index]
    hydrogens = mol.implicit_hydrogens(index)
    bonds = mol.bonds_of(index)
    neighbours = [mol.atoms[b.other(index)] for b in bonds]

    if atom.charge < 0:
        return -1.0000
    if atom.aromatic:
        return 0.1129
    if any(b.order == 2 for b in bonds):
        carbonyl_carbon = mol.atoms[bonds[0].other(index)]
        return 0.1129 if carbonyl_carbon.aromatic else -0.1526
    if hydrogens:
        for neighbour in neighbours:
            if neighbour.element == "C" and any(
                b.order == 2 and mol.atoms[b.other(neighbour.index)].element == "O"
                for b in mol.bonds_of(neighbour.index)
            ):
                return -0.3339  # carboxylic acid hydroxyl
        return -0.2893
    return 0.1350 if any(n.aromatic for n in neighbours) else 0.1129


def clogp(mol: Molecule) -> float:
    """Wildman-Crippen atom-contribution logP (reduced atom typing)."""
    total = 0.0
    for atom in mol.atoms:
        element = atom.element
        if element == "C":
            total += _carbon_logp(mol, atom.index)
        elif element == "N":
            total += _nitrogen_logp(mol, atom.index)
        elif element == "O":
            total += _oxygen_logp(mol, atom.index)
        elif element == "F":
            total += 0.4202
        elif element == "Cl":
            total += 0.6895
        elif element == "Br":
            total += 0.8456
        elif element == "I":
            total += 0.8857
        elif element == "S":
            oxides = sum(1 for b in mol.bonds_of(atom.index)
                         if b.order == 2 and mol.atoms[b.other(atom.index)].element == "O")
            total += -0.0024 if oxides else (0.6237 if atom.aromatic else 0.6482)
        elif element == "P":
            total += 0.8612
        elif element == "B":
            total += -0.3080
        total += mol.implicit_hydrogens(atom.index) * _H_CONTRIBUTION.get(element, 0.1230)
    return round(total, 3)


# ---------------------------------------------------------------------------
# Shape and flexibility
# ---------------------------------------------------------------------------

def rotatable_bonds(mol: Molecule) -> int:
    """Rotatable-bond count under the strict (Veber/RDKit) definition.

    An acyclic single bond between two non-terminal heavy atoms, excluding two
    classes whose rotation is not conformationally meaningful:

    * **conjugated carbonyl bonds** -- a trigonal carbon bearing a double bond
      to N, O, or S, joined to a non-terminal N, O, or S. This covers amides,
      esters, carbamates, anhydrides, and thioesters, all of which have
      restricted rotation from delocalisation. Excluding only amides
      over-counted generated structures by nearly three bonds each, because
      fragment assembly produces carbamates and anhydrides freely.
    * **symmetric terminal rotors** -- CF3, CCl3, CBr3, and tert-butyl, where
      rotation permutes indistinguishable substituents.
    """
    ring_bond = _ring_bonds(mol)
    count = 0
    for bond in mol.bonds:
        if bond.order != 1 or bond.aromatic:
            continue
        if frozenset((bond.a, bond.b)) in ring_bond:
            continue
        a, b = mol.atoms[bond.a], mol.atoms[bond.b]
        if a.is_attachment or b.is_attachment:
            continue
        if mol.degree(bond.a) < 2 or mol.degree(bond.b) < 2:
            continue
        if _is_conjugated_carbonyl_bond(mol, bond.a, bond.b):
            continue
        if _is_symmetric_rotor(mol, bond.a) or _is_symmetric_rotor(mol, bond.b):
            continue
        count += 1
    return count


def _ring_bonds(mol: Molecule) -> set[frozenset[int]]:
    """Bonds lying on a ring, taken from the cycle basis."""
    bonds: set[frozenset[int]] = set()
    for ring in ring_systems(mol):
        members = set(ring)
        for bond in mol.bonds:
            if bond.a in members and bond.b in members:
                bonds.add(frozenset((bond.a, bond.b)))
    return bonds


def _is_conjugated_carbonyl_bond(mol: Molecule, first: int, second: int) -> bool:
    """Amide, ester, carbamate, anhydride, urea, or thioester bond.

    A trigonal carbon double-bonded to N, O, or S, joined to a non-terminal N,
    O, or S. Delocalisation restricts rotation about that bond.

    RDKit's strict definition is stricter still, and *atom*-centric: it also
    disqualifies the conjugated carbon from ending any other bond, which
    removes the central C-C of an oxamide. Matching that exactly was tried and
    over-corrected in the other direction, so this stops at the bond-level
    rule. The residual disagreement is measured rather than assumed -- see
    docs/chemistry_validation.json -- and it is why RDKit is the preferred
    backend wherever it can be installed.
    """
    for carbon, hetero in ((first, second), (second, first)):
        if mol.atoms[carbon].element != "C" or mol.degree(carbon) != 3:
            continue
        double_bonded = any(
            b.order == 2 and mol.atoms[b.other(carbon)].element in ("N", "O", "S")
            for b in mol.bonds_of(carbon)
        )
        if not double_bonded:
            continue
        if mol.atoms[hetero].element in ("N", "O", "S") and mol.degree(hetero) > 1:
            return True
    return False


def _is_symmetric_rotor(mol: Molecule, index: int) -> bool:
    """A carbon carrying three identical terminal substituents."""
    if mol.atoms[index].element != "C":
        return False
    neighbours = mol.neighbors(index)
    elements = [mol.atoms[n].element for n in neighbours]
    if any(elements.count(halogen) == 3 for halogen in ("F", "Cl", "Br")):
        return True
    methyls = sum(
        1 for n in neighbours
        if mol.atoms[n].element == "C" and mol.degree(n) == 1
    )
    return methyls == 3


def fraction_sp3(mol: Molecule) -> float:
    """Fraction of carbons that are sp3 -- a saturation measure that tracks
    with clinical success (Lovering's 'escape from flatland')."""
    carbons = [a for a in mol.atoms if a.element == "C"]
    if not carbons:
        return 0.0
    sp3 = sum(
        1 for a in carbons
        if not a.aromatic and all(b.order == 1 for b in mol.bonds_of(a.index))
    )
    return round(sp3 / len(carbons), 4)


def heteroatom_fraction(mol: Molecule) -> float:
    heavy = mol.heavy_atoms
    if not heavy:
        return 0.0
    return round(sum(1 for a in heavy if a.element in HETEROATOMS) / len(heavy), 4)


def stereocentre_estimate(mol: Molecule) -> int:
    """Count sp3 carbons with four distinct neighbour environments.

    Neighbour identity is approximated by a two-bond environment hash rather
    than full CIP perception, so this over-counts symmetric cases.  It feeds
    the synthetic-accessibility proxy only, never a structural claim.
    """
    count = 0
    for atom in mol.atoms:
        if atom.element != "C" or atom.aromatic:
            continue
        bonds = mol.bonds_of(atom.index)
        if any(b.order != 1 for b in bonds):
            continue
        substituents = [b.other(atom.index) for b in bonds]
        hydrogens = mol.implicit_hydrogens(atom.index)
        if len(substituents) + hydrogens != 4 or hydrogens > 1:
            continue
        environments = {_environment_hash(mol, s, atom.index) for s in substituents}
        if len(environments) + hydrogens == 4:
            count += 1
    return count


def _environment_hash(mol: Molecule, index: int, exclude: int) -> str:
    atom = mol.atoms[index]
    neighbours = sorted(
        mol.atoms[n].element for n in mol.neighbors(index) if n != exclude
    )
    return f"{atom.element}{'a' if atom.aromatic else ''}:{''.join(neighbours)}"


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Descriptors:
    """A molecule's computed property vector."""

    smiles: str
    formula: str
    molecular_weight: float
    heavy_atoms: int
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    aromatic_rings: int
    rings: int
    fraction_sp3: float
    heteroatom_fraction: float
    formal_charge: int
    stereocentres: int
    backend: str = "local"

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def compute_descriptors(molecule: Molecule | str, backend: str | None = None) -> Descriptors:
    """Compute the descriptor vector, preferring RDKit when it is available.

    ``backend`` overrides the automatic choice; see
    ``core.chemistry.backend``. The chosen backend is recorded on the returned
    vector, because the two do not agree exactly and a number without its
    provenance cannot be interpreted.
    """
    chosen = chemistry_backend.active_backend(backend)
    if chosen == "rdkit":
        mol = parse_smiles(molecule) if isinstance(molecule, str) else molecule
        try:
            values = chemistry_backend.rdkit_descriptor_values(mol.to_smiles())
        except Exception:
            # RDKit rejecting a structure the local model accepts is itself a
            # finding, but it must not abort a campaign mid-search.
            return _local_descriptors(mol)
        values["heteroatom_fraction"] = heteroatom_fraction(mol)
        return Descriptors(backend="rdkit", **values)
    return _local_descriptors(molecule)


def _local_descriptors(molecule: Molecule | str) -> Descriptors:
    """The dependency-free implementation. See docs/chemistry_validation.json
    for its measured deviation from RDKit."""
    mol = parse_smiles(molecule) if isinstance(molecule, str) else molecule
    return Descriptors(
        backend="local",
        smiles=mol.to_smiles(),
        formula=mol.molecular_formula(),
        molecular_weight=mol.molecular_weight(),
        heavy_atoms=len(mol.heavy_atoms),
        clogp=clogp(mol),
        tpsa=tpsa(mol),
        hbd=hbd_count(mol),
        hba=hba_count(mol),
        rotatable_bonds=rotatable_bonds(mol),
        aromatic_rings=len(aromatic_rings(mol)),
        rings=len(ring_systems(mol)),
        fraction_sp3=fraction_sp3(mol),
        heteroatom_fraction=heteroatom_fraction(mol),
        formal_charge=mol.net_charge(),
        stereocentres=stereocentre_estimate(mol),
    )

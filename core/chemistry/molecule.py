"""A dependency-free molecular graph: SMILES in, molecule out, SMILES back.

Why this module exists
----------------------
The rest of the engine needs to *emit* structures that do not exist yet, not
merely read curated ones.  That requires a writer as well as a parser, plus a
join operation that fuses two fragments at marked attachment points.

RDKit would be the obvious choice and is strictly better where it is
available.  It has no wheel for the interpreter this project pins, so the
graph, valence model, and SMILES round-trip are implemented here against the
OpenSMILES specification, restricted to the subset a small-molecule design
campaign actually needs:

  * organic-subset atoms and bracket atoms (charge, explicit H, any element)
  * single / double / triple / aromatic bonds, branches, ring-bond closures
  * lowercase aromaticity as written -- this module trusts the input's
    aromatic perception rather than re-deriving Huckel rings
  * ``*`` attachment points, used by the fragment assembler

Not modelled: stereochemistry (``@``/``@@`` and ``/``\\`` are parsed and
discarded), isotopes, and reaction SMILES.  Every downstream descriptor is
therefore stereochemistry-blind, which is stated wherever it matters.

If RDKit *is* importable, :func:`canonical_smiles` and :func:`is_valid`
delegate to it, so a richer environment automatically gets stricter
validation without any other module changing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Iterator, Sequence

try:  # pragma: no cover - environment dependent
    from rdkit import Chem as _rdkit_chem
    from rdkit import RDLogger as _rdkit_logger

    _rdkit_logger.DisableLog("rdApp.*")
except Exception:  # pragma: no cover
    _rdkit_chem = None

HAS_RDKIT = _rdkit_chem is not None

ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}
AROMATIC_SYMBOLS = {"b": "B", "c": "C", "n": "N", "o": "O", "p": "P", "s": "S"}

# Average atomic masses (IUPAC 2021), sufficient for molecular-weight filters.
ATOMIC_MASS = {
    "H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Na": 22.990, "Mg": 24.305, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Zn": 65.38,
    "Br": 79.904, "I": 126.904, "*": 0.0,
}

# Allowed valences, smallest first; the parser picks the smallest that covers
# the bonds actually drawn, which is what gives N(=O)=O five bonds while
# leaving an amine at three.
VALENCES = {
    "B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5), "S": (2, 4, 6),
    "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,), "Si": (4,), "*": (1,),
}

# Elements whose valence rises with positive charge (pnictogens/chalcogens)
# rather than falling with any charge (carbon and friends).
_LONE_PAIR_DONORS = {"N", "P", "O", "S"}

BOND_SYMBOL = {1: "", 2: "=", 3: "#"}


class SmilesError(ValueError):
    """Raised when a SMILES string cannot be parsed under this subset."""


@dataclass
class Atom:
    element: str
    aromatic: bool = False
    charge: int = 0
    explicit_h: int | None = None  # from a bracket atom; None means "infer"
    index: int = 0

    @property
    def is_attachment(self) -> bool:
        return self.element == "*"

    @property
    def mass(self) -> float:
        return ATOMIC_MASS.get(self.element, 0.0)


@dataclass
class Bond:
    a: int
    b: int
    order: int = 1
    aromatic: bool = False

    def other(self, index: int) -> int:
        return self.b if index == self.a else self.a


@dataclass
class Molecule:
    """An undirected molecular graph with an inferred hydrogen count."""

    atoms: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)
    name: str = ""

    # -- construction ---------------------------------------------------
    def add_atom(self, atom: Atom) -> int:
        atom.index = len(self.atoms)
        self.atoms.append(atom)
        return atom.index

    def add_bond(self, a: int, b: int, order: int = 1, aromatic: bool = False) -> None:
        if a == b:
            raise SmilesError("An atom cannot bond to itself")
        if self.bond_between(a, b) is not None:
            raise SmilesError(f"Duplicate bond between atoms {a} and {b}")
        self.bonds.append(Bond(a, b, order, aromatic))

    def copy(self) -> "Molecule":
        return Molecule(
            atoms=[Atom(a.element, a.aromatic, a.charge, a.explicit_h, a.index) for a in self.atoms],
            bonds=[Bond(b.a, b.b, b.order, b.aromatic) for b in self.bonds],
            name=self.name,
        )

    # -- topology -------------------------------------------------------
    def neighbors(self, index: int) -> list[int]:
        return [b.other(index) for b in self.bonds if index in (b.a, b.b)]

    def bonds_of(self, index: int) -> list[Bond]:
        return [b for b in self.bonds if index in (b.a, b.b)]

    def bond_between(self, a: int, b: int) -> Bond | None:
        for bond in self.bonds:
            if {bond.a, bond.b} == {a, b}:
                return bond
        return None

    def degree(self, index: int) -> int:
        return len(self.bonds_of(index))

    @property
    def heavy_atoms(self) -> list[Atom]:
        return [a for a in self.atoms if not a.is_attachment]

    @property
    def attachment_points(self) -> list[int]:
        return [a.index for a in self.atoms if a.is_attachment]

    def __len__(self) -> int:
        return len(self.heavy_atoms)

    # -- hydrogens ------------------------------------------------------
    def bond_order_sum(self, index: int) -> int:
        """Summed bond order, counting an aromatic ring bond as one.

        An aromatic atom is credited one extra unit for the pi bond it carries
        in the delocalised system -- that credit is what gives benzene's
        carbons exactly one hydrogen each and pyridine's nitrogen none.  The
        credit is withheld in the two cases where no such pi bond exists:

        * an atom already carrying an explicit double or triple bond, as in
          the exocyclic carbonyls of caffeine's purine ring;
        * a pyrrole-type nitrogen or phosphorus (three connections) and the
          chalcogens, which donate a lone pair to the ring instead.
        """
        atom = self.atoms[index]
        bonds = self.bonds_of(index)
        total = sum(b.order for b in bonds)
        return total + (1 if self._has_delocalised_pi_bond(index) else 0)

    def _has_delocalised_pi_bond(self, index: int) -> bool:
        atom = self.atoms[index]
        if not atom.aromatic:
            return False
        bonds = self.bonds_of(index)
        if any(b.order > 1 for b in bonds):
            return False
        if atom.element in ("O", "S"):
            return False
        if atom.element in ("N", "P") and len(bonds) >= 3:
            return False
        return True

    def implicit_hydrogens(self, index: int) -> int:
        atom = self.atoms[index]
        if atom.explicit_h is not None:
            return atom.explicit_h
        options = VALENCES.get(atom.element)
        if options is None:  # an element outside the model carries no implicit H
            return 0
        if atom.element in _LONE_PAIR_DONORS:
            options = tuple(v + atom.charge for v in options)
        else:
            options = tuple(v - abs(atom.charge) for v in options)
        used = self.bond_order_sum(index)
        for value in options:
            if value >= used:
                return max(0, value - used)
        return 0

    def total_hydrogens(self) -> int:
        return sum(self.implicit_hydrogens(a.index) for a in self.atoms)

    # -- identity -------------------------------------------------------
    def molecular_formula(self) -> str:
        """Formula in Hill order: carbon, hydrogen, then everything alphabetical."""
        counts: dict[str, int] = {}
        for atom in self.atoms:
            if atom.is_attachment:
                continue
            counts[atom.element] = counts.get(atom.element, 0) + 1
        hydrogens = self.total_hydrogens()
        if hydrogens:
            counts["H"] = counts.get("H", 0) + hydrogens
        ordered: list[str] = []
        for symbol in ("C", "H"):
            if counts.pop(symbol, 0):
                pass
        # recount because pop above discarded the values
        counts = {}
        for atom in self.atoms:
            if atom.is_attachment:
                continue
            counts[atom.element] = counts.get(atom.element, 0) + 1
        if hydrogens:
            counts["H"] = counts.get("H", 0) + hydrogens
        head = [s for s in ("C", "H") if s in counts] if "C" in counts else []
        tail = sorted(s for s in counts if s not in head)
        for symbol in [*head, *tail]:
            n = counts[symbol]
            ordered.append(symbol if n == 1 else f"{symbol}{n}")
        return "".join(ordered)

    def molecular_weight(self) -> float:
        heavy = sum(a.mass for a in self.atoms if not a.is_attachment)
        return round(heavy + self.total_hydrogens() * ATOMIC_MASS["H"], 3)

    def net_charge(self) -> int:
        return sum(a.charge for a in self.atoms)

    def to_smiles(self) -> str:
        return write_smiles(self)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Molecule({self.molecular_formula()!r}, {len(self.atoms)} atoms)"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_BOND_ORDERS = {"-": 1, "=": 2, "#": 3, "/": 1, "\\": 1}


def parse_smiles(text: str, name: str = "") -> Molecule:
    """Parse a SMILES string into a :class:`Molecule`.

    Raises :class:`SmilesError` on anything this subset cannot represent, so a
    generator that assembles a nonsensical string fails loudly rather than
    scoring a molecule that does not exist.
    """
    smiles = text.strip()
    if not smiles:
        raise SmilesError("Empty SMILES string")

    mol = Molecule(name=name)
    previous: int | None = None
    branch_stack: list[int] = []
    ring_open: dict[int, tuple[int, int | None]] = {}
    pending_order: int | None = None
    pending_aromatic = False
    i = 0
    n = len(smiles)

    while i < n:
        char = smiles[i]

        if char in "-=#/\\":
            pending_order = _BOND_ORDERS[char]
            i += 1
            continue
        if char == ":":
            pending_order, pending_aromatic = 1, True
            i += 1
            continue
        if char == ".":
            previous, pending_order, pending_aromatic = None, None, False
            i += 1
            continue
        if char == "(":
            if previous is None:
                raise SmilesError("Branch opened before any atom")
            branch_stack.append(previous)
            i += 1
            continue
        if char == ")":
            if not branch_stack:
                raise SmilesError("Unbalanced ')' in SMILES")
            previous = branch_stack.pop()
            i += 1
            continue
        if char.isdigit() or char == "%":
            if char == "%":
                if i + 2 >= n or not smiles[i + 1:i + 3].isdigit():
                    raise SmilesError("Malformed '%nn' ring-bond number")
                number, i = int(smiles[i + 1:i + 3]), i + 3
            else:
                number, i = int(char), i + 1
            if previous is None:
                raise SmilesError("Ring-bond number before any atom")
            if number in ring_open:
                partner, opening_order = ring_open.pop(number)
                order = pending_order or opening_order or 1
                aromatic = mol.atoms[partner].aromatic and mol.atoms[previous].aromatic
                mol.add_bond(partner, previous, order, aromatic and order == 1)
            else:
                ring_open[number] = (previous, pending_order)
            pending_order, pending_aromatic = None, False
            continue

        atom, i = _read_atom(smiles, i)
        index = mol.add_atom(atom)
        if previous is not None:
            both_aromatic = mol.atoms[previous].aromatic and atom.aromatic
            order = pending_order or 1
            aromatic = pending_aromatic or (both_aromatic and pending_order is None)
            mol.add_bond(previous, index, order, aromatic)
        previous = index
        pending_order, pending_aromatic = None, False

    if branch_stack:
        raise SmilesError("Unbalanced '(' in SMILES")
    if ring_open:
        raise SmilesError(f"Unclosed ring-bond number(s): {sorted(ring_open)}")
    if not mol.atoms:
        raise SmilesError("SMILES contained no atoms")
    _check_valences(mol)
    return mol


def _read_atom(smiles: str, i: int) -> tuple[Atom, int]:
    char = smiles[i]
    if char == "[":
        end = smiles.find("]", i)
        if end == -1:
            raise SmilesError("Unclosed '[' in SMILES")
        return _read_bracket_atom(smiles[i + 1:end]), end + 1
    if char == "*":
        return Atom("*"), i + 1
    two = smiles[i:i + 2]
    if two in ("Cl", "Br"):
        return Atom(two), i + 2
    if char in ORGANIC_SUBSET:
        return Atom(char), i + 1
    if char in AROMATIC_SYMBOLS:
        return Atom(AROMATIC_SYMBOLS[char], aromatic=True), i + 1
    raise SmilesError(f"Unsupported character {char!r} at position {i} in {smiles!r}")


def _read_bracket_atom(body: str) -> Atom:
    if not body:
        raise SmilesError("Empty bracket atom '[]'")
    j = 0
    while j < len(body) and body[j].isdigit():  # isotope: parsed, then dropped
        j += 1
    if j >= len(body):
        raise SmilesError(f"Bracket atom [{body}] has no element")
    aromatic = False
    if body[j] == "*":
        element, j = "*", j + 1
    elif body[j:j + 2] in ATOMIC_MASS and body[j:j + 2].istitle():
        element, j = body[j:j + 2], j + 2
    elif body[j].isupper():
        element, j = body[j], j + 1
    elif body[j].islower():
        if body[j] not in AROMATIC_SYMBOLS:
            raise SmilesError(f"Unsupported aromatic element in [{body}]")
        element, aromatic, j = AROMATIC_SYMBOLS[body[j]], True, j + 1
    else:
        raise SmilesError(f"Cannot read element from [{body}]")

    explicit_h = 0
    charge = 0
    while j < len(body):
        c = body[j]
        if c in "@":  # stereochemistry: consumed and deliberately discarded
            j += 1
            while j < len(body) and body[j] == "@":
                j += 1
        elif c == "H":
            j += 1
            digits = ""
            while j < len(body) and body[j].isdigit():
                digits += body[j]
                j += 1
            explicit_h = int(digits) if digits else 1
        elif c in "+-":
            sign = 1 if c == "+" else -1
            j += 1
            digits = ""
            while j < len(body) and body[j].isdigit():
                digits += body[j]
                j += 1
            if digits:
                charge = sign * int(digits)
            else:
                repeats = 1
                while j < len(body) and body[j] == c:
                    repeats += 1
                    j += 1
                charge = sign * repeats
        elif c == ":":  # atom map number, irrelevant here
            j += 1
            while j < len(body) and body[j].isdigit():
                j += 1
        else:
            raise SmilesError(f"Unsupported token {c!r} in bracket atom [{body}]")
    return Atom(element, aromatic=aromatic, charge=charge, explicit_h=explicit_h)


def _check_valences(mol: Molecule) -> None:
    """Reject structures whose drawn bonds exceed any allowed valence."""
    for atom in mol.atoms:
        options = VALENCES.get(atom.element)
        if options is None:
            continue
        if atom.element in _LONE_PAIR_DONORS:
            limit = max(options) + atom.charge
        else:
            limit = max(options) - abs(atom.charge)
        used = mol.bond_order_sum(atom.index) + (atom.explicit_h or 0)
        if used > limit:
            raise SmilesError(
                f"Atom {atom.index} ({atom.element}) carries {used} bonds, above its "
                f"maximum of {limit}"
            )


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def write_smiles(mol: Molecule) -> str:
    """Serialise a molecule back to SMILES via a depth-first traversal.

    The traversal runs twice over the same deterministic ordering: once to
    split the edges into a spanning forest plus back-edges, and once to emit.
    Back-edges become ring-bond numbers, written with their bond symbol at the
    opening atom only.  The output is valid and round-trips, but it is *not* a
    canonical SMILES: use :func:`canonical_smiles` when identity comparison is
    the point.
    """
    if not mol.atoms:
        return ""

    visited: set[int] = set()
    seen_bonds: set[int] = set()
    children: dict[int, list[tuple[int, Bond]]] = {a.index: [] for a in mol.atoms}
    ring_labels: dict[int, int] = {}
    ring_at_atom: dict[int, list[tuple[int, Bond, bool]]] = {a.index: [] for a in mol.atoms}
    counter = 0

    def ordered_bonds(node: int) -> list[Bond]:
        return sorted(mol.bonds_of(node), key=lambda b: b.other(node))

    def classify(node: int) -> None:
        nonlocal counter
        visited.add(node)
        for bond in ordered_bonds(node):
            if id(bond) in seen_bonds:
                continue
            seen_bonds.add(id(bond))
            neighbour = bond.other(node)
            if neighbour in visited:
                counter += 1
                ring_labels[id(bond)] = counter
                ring_at_atom[node].append((counter, bond, True))
                ring_at_atom[neighbour].append((counter, bond, False))
            else:
                children[node].append((neighbour, bond))
                classify(neighbour)

    for root in range(len(mol.atoms)):
        if root not in visited:
            classify(root)

    def bond_prefix(bond: Bond, a: int, b: int) -> str:
        if bond.aromatic:
            return ""
        if bond.order == 1 and mol.atoms[a].aromatic and mol.atoms[b].aromatic:
            return "-"  # a biaryl single bond must be written explicitly
        return BOND_SYMBOL[bond.order]

    parts: list[str] = []

    def emit(node: int) -> None:
        parts.append(_atom_token(mol, node))
        for label, bond, is_opening in sorted(ring_at_atom[node]):
            prefix = bond_prefix(bond, bond.a, bond.b) if is_opening else ""
            parts.append(prefix + _ring_token(label))
        kids = children[node]
        for position, (child, bond) in enumerate(kids):
            last = position == len(kids) - 1
            if not last:
                parts.append("(")
            parts.append(bond_prefix(bond, node, child))
            emit(child)
            if not last:
                parts.append(")")

    emitted: set[int] = set()
    for root in range(len(mol.atoms)):
        if root in emitted:
            continue
        component = _component(children, root)
        if component & emitted:
            continue
        if parts:
            parts.append(".")
        emit(root)
        emitted |= component
    return "".join(parts)


def _component(children: dict[int, list[tuple[int, "Bond"]]], root: int) -> set[int]:
    seen = {root}
    stack = [root]
    while stack:
        node = stack.pop()
        for child, _ in children[node]:
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return seen


def _ring_token(number: int) -> str:
    return str(number) if number < 10 else f"%{number:02d}"


def _atom_token(mol: Molecule, index: int) -> str:
    atom = mol.atoms[index]
    if atom.is_attachment:
        return "*"
    symbol = atom.element.lower() if atom.aromatic else atom.element
    needs_bracket = (
        atom.charge != 0
        or atom.element not in ORGANIC_SUBSET
        or (atom.explicit_h is not None and atom.explicit_h != _inferred_h(mol, index))
    )
    if not needs_bracket:
        return symbol
    hydrogens = atom.explicit_h if atom.explicit_h is not None else _inferred_h(mol, index)
    h_part = "" if hydrogens == 0 else ("H" if hydrogens == 1 else f"H{hydrogens}")
    if atom.charge > 0:
        charge_part = "+" if atom.charge == 1 else f"+{atom.charge}"
    elif atom.charge < 0:
        charge_part = "-" if atom.charge == -1 else f"-{abs(atom.charge)}"
    else:
        charge_part = ""
    return f"[{symbol}{h_part}{charge_part}]"


def _inferred_h(mol: Molecule, index: int) -> int:
    """Hydrogen count the parser would infer if the atom carried no bracket."""
    atom = mol.atoms[index]
    options = VALENCES.get(atom.element)
    if options is None:
        return 0
    if atom.element in _LONE_PAIR_DONORS:
        options = tuple(v + atom.charge for v in options)
    else:
        options = tuple(v - abs(atom.charge) for v in options)
    used = mol.bond_order_sum(index)
    for value in options:
        if value >= used:
            return max(0, value - used)
    return 0


# ---------------------------------------------------------------------------
# Fragment assembly
# ---------------------------------------------------------------------------

# A single bond between two of these, formed at a junction, gives a peroxide,
# hydrazine, hydroxylamine, disulfide, or sulfenamide. Each is either unstable,
# genotoxic, or both, and none is what the fragments intended.
_UNSTABLE_JUNCTION_ELEMENTS = {"N", "O", "S"}


def attach(core: Molecule, fragment: Molecule, core_point: int | None = None,
           fragment_point: int | None = None,
           allow_heteroatom_junction: bool = False) -> Molecule:
    """Fuse ``fragment`` onto ``core`` at one ``*`` attachment point on each.

    Both attachment atoms are consumed; a bond is formed between the heavy
    atoms they were marked on.  Remaining attachment points survive, which is
    what lets a multi-armed scaffold be decorated one arm at a time.

    Combinatorial assembly happily produces bonds no chemist would draw: join
    an ether linker to an arm that begins with oxygen and the result is a
    peroxide.  By default a new single bond between two heteroatoms is
    refused, which rules out peroxides, hydrazines, hydroxylamines,
    disulfides, and sulfenamides at the junction.  Pass
    ``allow_heteroatom_junction=True`` to build one deliberately.
    """
    core_points = core.attachment_points
    fragment_points = fragment.attachment_points
    if not core_points:
        raise SmilesError("Core has no '*' attachment point")
    if not fragment_points:
        raise SmilesError("Fragment has no '*' attachment point")
    core_point = core_points[0] if core_point is None else core_point
    fragment_point = fragment_points[0] if fragment_point is None else fragment_point

    core_anchor = _sole_neighbour(core, core_point)
    fragment_anchor = _sole_neighbour(fragment, fragment_point)

    merged = Molecule(name=core.name)
    core_map: dict[int, int] = {}
    for atom in core.atoms:
        if atom.index == core_point:
            continue
        core_map[atom.index] = merged.add_atom(
            Atom(atom.element, atom.aromatic, atom.charge, atom.explicit_h)
        )
    fragment_map: dict[int, int] = {}
    for atom in fragment.atoms:
        if atom.index == fragment_point:
            continue
        fragment_map[atom.index] = merged.add_atom(
            Atom(atom.element, atom.aromatic, atom.charge, atom.explicit_h)
        )
    for bond in core.bonds:
        if core_point in (bond.a, bond.b):
            continue
        merged.add_bond(core_map[bond.a], core_map[bond.b], bond.order, bond.aromatic)
    for bond in fragment.bonds:
        if fragment_point in (bond.a, bond.b):
            continue
        merged.add_bond(fragment_map[bond.a], fragment_map[bond.b], bond.order, bond.aromatic)
    if not allow_heteroatom_junction:
        pair = {core.atoms[core_anchor].element, fragment.atoms[fragment_anchor].element}
        if pair <= _UNSTABLE_JUNCTION_ELEMENTS:
            raise SmilesError(
                f"Refusing to bond {core.atoms[core_anchor].element} to "
                f"{fragment.atoms[fragment_anchor].element} at an assembly junction: "
                "a heteroatom-heteroatom single bond here is a peroxide, hydrazine, "
                "hydroxylamine, disulfide, or sulfenamide"
            )
    merged.add_bond(core_map[core_anchor], fragment_map[fragment_anchor], 1, False)
    _check_valences(merged)
    return merged


def _sole_neighbour(mol: Molecule, index: int) -> int:
    neighbours = mol.neighbors(index)
    if len(neighbours) != 1:
        raise SmilesError(
            f"Attachment point {index} must have exactly one neighbour, found {len(neighbours)}"
        )
    return neighbours[0]


def cap_attachments(mol: Molecule, element: str = "H") -> Molecule:
    """Remove leftover ``*`` points, capping each with hydrogen (or a group).

    A designed molecule must not ship with dangling valences; scoring and
    formula reporting run on the capped form.
    """
    capped = mol.copy()
    points = capped.attachment_points
    if not points:
        return capped
    keep = [a for a in capped.atoms if not a.is_attachment]
    remap = {}
    rebuilt = Molecule(name=capped.name)
    for atom in keep:
        remap[atom.index] = rebuilt.add_atom(
            Atom(atom.element, atom.aromatic, atom.charge, atom.explicit_h)
        )
    for bond in capped.bonds:
        if bond.a in points or bond.b in points:
            continue
        rebuilt.add_bond(remap[bond.a], remap[bond.b], bond.order, bond.aromatic)
    if element != "H":
        for point in points:
            anchor = remap[_sole_neighbour(capped, point)]
            new_atom = rebuilt.add_atom(Atom(element))
            rebuilt.add_bond(anchor, new_atom, 1, False)
    return rebuilt


# ---------------------------------------------------------------------------
# Optional RDKit delegation
# ---------------------------------------------------------------------------

def canonical_smiles(smiles: str) -> str:
    """Canonical form when RDKit is present; a round-tripped form otherwise.

    The fallback is deterministic but only canonical with respect to this
    module's traversal, so it is safe for caching and round-trip tests and
    unsafe as a cross-toolkit identity key.
    """
    if HAS_RDKIT:  # pragma: no cover - environment dependent
        mol = _rdkit_chem.MolFromSmiles(smiles)
        if mol is None:
            raise SmilesError(f"RDKit rejected SMILES: {smiles}")
        return _rdkit_chem.MolToSmiles(mol)
    return write_smiles(parse_smiles(smiles))


def is_valid(smiles: str) -> bool:
    try:
        canonical_smiles(smiles)
        return True
    except Exception:
        return False

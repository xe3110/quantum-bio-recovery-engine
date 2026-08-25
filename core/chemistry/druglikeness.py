"""Developability models: oral rules, CNS exposure, tractability, and alerts.

These are *filters and priors*, not predictions of behaviour in an animal.
Each one is a published heuristic or an explicitly-declared proxy:

``lipinski`` / ``veber``
    Lipinski (1997) and Veber (2002) oral-absorption rules, reported as the
    individual criteria plus a violation count rather than a pass/fail, so a
    deliberate rule-breaker is visible instead of silently discarded.

``cns_mpo``
    Wager et al. (2010), *ACS Chem. Neurosci.* 1:435 -- six desirability
    functions summing to 6.  Two of its six inputs, logD(7.4) and the most
    basic pKa, are measurements this engine has no way to make, so both are
    estimated from functional-group class (see :func:`basic_pka_estimate`).
    A CNS MPO computed from estimated pKa is a ranking aid within one
    campaign, not a number to compare against published MPO values.

``synthetic_tractability``
    A declared proxy, *not* the Ertl-Schuffenhauer SA score, which needs a
    fragment-frequency corpus this project does not ship.  It penalises size,
    stereocentres, and unusual ring topology -- the things that make a first
    synthesis slow -- and should be read as "cheap to try" rather than
    "makeable".

``structural_alerts``
    Reactive and assay-interfering motifs, matched by explicit graph
    predicates rather than SMARTS so every rule is readable and testable.
    Flagging a motif is a reason to look, not a verdict: several approved
    drugs carry alerts deliberately (a covalent warhead is the point).
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Callable

from core.chemistry.descriptors import (
    Descriptors, HALOGENS, compute_descriptors, ring_systems, aromatic_rings,
)
from core.chemistry.molecule import Molecule, parse_smiles


def _as_molecule(molecule: Molecule | str) -> Molecule:
    return parse_smiles(molecule) if isinstance(molecule, str) else molecule


# ---------------------------------------------------------------------------
# Oral-absorption rule sets
# ---------------------------------------------------------------------------

def lipinski(desc: Descriptors) -> dict[str, Any]:
    """Lipinski's rule of five, reported criterion by criterion."""
    criteria = {
        "mw_le_500": desc.molecular_weight <= 500,
        "clogp_le_5": desc.clogp <= 5,
        "hbd_le_5": desc.hbd <= 5,
        "hba_le_10": desc.hba <= 10,
    }
    violations = sum(1 for passed in criteria.values() if not passed)
    return {**criteria, "violations": violations, "compliant": violations <= 1}


def veber(desc: Descriptors) -> dict[str, Any]:
    """Veber's oral-bioavailability rules: flexibility and polar surface."""
    criteria = {
        "rotatable_le_10": desc.rotatable_bonds <= 10,
        "tpsa_le_140": desc.tpsa <= 140,
    }
    return {**criteria, "compliant": all(criteria.values())}


# ---------------------------------------------------------------------------
# Ionisation and CNS exposure
# ---------------------------------------------------------------------------

# Class-level basic pKa priors. Coarse by design: within a class the spread is
# 1-2 log units, which is smaller than the difference between classes.
_PKA_PRIORS = (
    ("amidine_guanidine", 11.5),
    ("primary_aliphatic_amine", 10.6),
    ("secondary_aliphatic_amine", 10.8),
    ("tertiary_aliphatic_amine", 9.8),
    ("piperazine", 9.0),
    ("morpholine", 8.4),
    ("imidazole", 6.9),
    ("pyridine", 5.2),
    ("aniline", 4.6),
)


def _amine_class(mol: Molecule, index: int) -> str | None:
    atom = mol.atoms[index]
    if atom.element != "N":
        return None
    bonds = mol.bonds_of(index)
    neighbours = [mol.atoms[b.other(index)] for b in bonds]
    hydrogens = mol.implicit_hydrogens(index)

    # An amide, sulfonamide, or nitrile nitrogen is not appreciably basic.
    for bond in bonds:
        other = bond.other(index)
        neighbour = mol.atoms[other]
        if neighbour.element == "C" and any(
            b.order == 2 and mol.atoms[b.other(other)].element in ("O", "S")
            for b in mol.bonds_of(other)
        ):
            if not any(mol.atoms[b.other(index)].element == "N" for b in bonds if b.other(index) != other):
                return None
            return "amidine_guanidine"
        if neighbour.element == "S" and sum(
            1 for b in mol.bonds_of(other) if b.order == 2 and mol.atoms[b.other(other)].element == "O"
        ) >= 1:
            return None
    if any(b.order >= 2 for b in bonds) and not atom.aromatic:
        return None

    if atom.aromatic:
        ring_elements = [
            mol.atoms[i].element
            for ring in aromatic_rings(mol) if index in ring for i in ring
        ]
        if len(bonds) >= 3:
            return None  # pyrrole-type nitrogen: not basic
        return "imidazole" if ring_elements.count("N") >= 2 else "pyridine"

    if any(n.aromatic for n in neighbours):
        return "aniline"

    ring_membership = [ring for ring in ring_systems(mol) if index in ring]
    if ring_membership:
        ring = ring_membership[0]
        elements = [mol.atoms[i].element for i in ring]
        if elements.count("N") >= 2:
            return "piperazine"
        if "O" in elements:
            return "morpholine"
    return {2: "primary_aliphatic_amine", 1: "secondary_aliphatic_amine"}.get(
        hydrogens, "tertiary_aliphatic_amine"
    )


def basic_pka_estimate(molecule: Molecule | str) -> tuple[float, str]:
    """Estimate the most basic pKa and name the group it came from.

    Returns ``(0.0, "none")`` for a molecule with no basic centre, which the
    CNS MPO treats as fully desirable on that axis.
    """
    mol = _as_molecule(molecule)
    priors = dict(_PKA_PRIORS)
    best, source = 0.0, "none"
    for atom in mol.atoms:
        group = _amine_class(mol, atom.index)
        if group is None:
            continue
        value = priors.get(group, 0.0)
        if value > best:
            best, source = value, group
    return round(best, 2), source


def logd_estimate(clogp: float, pka: float, ph: float = 7.4) -> float:
    """Distribution coefficient for a monoprotic base via Henderson-Hasselbalch.

    Neutral and acidic molecules are returned unchanged, so this is an
    estimate of the *basic* ionisation penalty only.
    """
    if pka <= 0:
        return round(clogp, 3)
    return round(clogp - _log10(1 + 10 ** (pka - ph)), 3)


def _log10(value: float) -> float:
    from math import log10

    return log10(value)


def _hump(value: float, lower: float, upper: float) -> float:
    """1.0 inside [lower, upper], falling linearly to 0 over the same width."""
    if lower <= value <= upper:
        return 1.0
    width = max(upper - lower, 1e-9)
    if value < lower:
        return max(0.0, 1.0 - (lower - value) / width)
    return max(0.0, 1.0 - (value - upper) / width)


def _ramp(value: float, best: float, worst: float) -> float:
    """1.0 at or beyond ``best``, 0.0 at or beyond ``worst``, linear between."""
    if best < worst:
        if value <= best:
            return 1.0
        if value >= worst:
            return 0.0
        return (worst - value) / (worst - best)
    if value >= best:
        return 1.0
    if value <= worst:
        return 0.0
    return (value - worst) / (best - worst)


def acidic_centres(molecule: Molecule | str) -> list[str]:
    """Groups that carry a negative charge at pH 7.4, and their names.

    This exists because CNS MPO cannot see them.  Wager's score takes the
    *most basic* pKa as its ionisation term, which is the right choice for the
    basic and neutral compounds it was derived on and leaves acids entirely
    unpenalised: a carboxylic acid with low molecular weight and modest logP
    can score above 5 out of 6 while being almost entirely anionic at
    physiological pH and effectively barred from the brain.  The delivery gate
    consults this separately.
    """
    mol = _as_molecule(molecule)
    found: list[str] = []
    for atom in mol.atoms:
        if atom.element == "O" and atom.charge < 0:
            found.append("carboxylate")
            continue
        if atom.element != "O" or mol.implicit_hydrogens(atom.index) != 1:
            continue
        for neighbour in mol.neighbors(atom.index):
            element = mol.atoms[neighbour].element
            oxygens = sum(
                1 for b in mol.bonds_of(neighbour)
                if b.order == 2 and mol.atoms[b.other(neighbour)].element == "O"
            )
            if element == "C" and oxygens >= 1:
                found.append("carboxylic_acid")
            elif element == "S" and oxygens >= 2:
                found.append("sulfonic_acid")
            elif element == "P" and oxygens >= 1:
                found.append("phosphonic_acid")
    for ring in ring_systems(mol):
        elements = [mol.atoms[i].element for i in ring]
        if len(ring) == 5 and elements.count("N") == 4:
            found.append("tetrazole")
    return sorted(set(found))


def cns_mpo(desc: Descriptors, pka: float | None = None) -> dict[str, Any]:
    """Wager CNS multi-parameter optimisation score, 0-6 (higher is better).

    A total at or above 4.0 is the threshold the original work associates with
    a higher chance of adequate brain exposure and a lower alignment with
    safety attrition.
    """
    if pka is None:
        pka = 0.0
    logd = logd_estimate(desc.clogp, pka)
    components = {
        "clogp": _ramp(desc.clogp, 3.0, 5.0),
        "clogd": _ramp(logd, 2.0, 4.0),
        "mw": _ramp(desc.molecular_weight, 360.0, 500.0),
        "tpsa": _hump(desc.tpsa, 40.0, 90.0),
        "hbd": _ramp(desc.hbd, 0.5, 3.5),
        "pka": _ramp(pka, 8.0, 10.0),
    }
    total = sum(components.values())
    return {
        "components": {k: round(v, 3) for k, v in components.items()},
        "estimated_logd74": logd,
        "estimated_basic_pka": pka,
        "total": round(total, 3),
        "normalised": round(total / 6.0, 4),
        "above_threshold": total >= 4.0,
        "caveat": (
            "Computed from an estimated basic pKa and the logD derived from it, not from "
            "measurements. The score is blind to acidic ionisation by construction; see "
            "acidic_centres()."
        ),
    }


# ---------------------------------------------------------------------------
# Structural alerts
# ---------------------------------------------------------------------------

def _has_carbonyl(mol: Molecule, carbon: int) -> bool:
    return mol.atoms[carbon].element == "C" and any(
        b.order == 2 and mol.atoms[b.other(carbon)].element == "O" for b in mol.bonds_of(carbon)
    )


def _alert_acyl_halide(mol: Molecule) -> bool:
    return any(
        _has_carbonyl(mol, a.index)
        and any(mol.atoms[n].element in HALOGENS for n in mol.neighbors(a.index))
        for a in mol.atoms if a.element == "C"
    )


def _alert_michael_acceptor(mol: Molecule) -> bool:
    """A C=C conjugated to a carbonyl, nitrile, or sulfonyl: a covalent warhead."""
    for bond in mol.bonds:
        if bond.order != 2:
            continue
        a, b = mol.atoms[bond.a], mol.atoms[bond.b]
        if a.aromatic or b.aromatic or {a.element, b.element} != {"C"}:
            continue
        for end, far in ((bond.a, bond.b), (bond.b, bond.a)):
            for neighbour in mol.neighbors(end):
                if neighbour == far:
                    continue
                atom = mol.atoms[neighbour]
                if atom.element == "C" and _has_carbonyl(mol, neighbour):
                    return True
                if atom.element == "S" and sum(
                    1 for x in mol.bonds_of(neighbour)
                    if x.order == 2 and mol.atoms[x.other(neighbour)].element == "O"
                ) >= 2:
                    return True
                if atom.element == "C" and any(
                    x.order == 3 and mol.atoms[x.other(neighbour)].element == "N"
                    for x in mol.bonds_of(neighbour)
                ):
                    return True
    return False


def _alert_small_strained_ring(mol: Molecule) -> bool:
    """Epoxide or aziridine: a hard electrophile and a genotoxicity concern."""
    return any(
        len(ring) == 3 and any(mol.atoms[i].element in ("O", "N") for i in ring)
        for ring in ring_systems(mol)
    )


def _alert_aldehyde(mol: Molecule) -> bool:
    return any(
        a.element == "C" and _has_carbonyl(mol, a.index)
        and mol.implicit_hydrogens(a.index) >= 1 and not a.aromatic
        for a in mol.atoms
    )


def _alert_hydrazine(mol: Molecule) -> bool:
    return any(
        bond.order == 1
        and mol.atoms[bond.a].element == "N" and mol.atoms[bond.b].element == "N"
        and not mol.atoms[bond.a].aromatic and not mol.atoms[bond.b].aromatic
        and not any(_has_carbonyl(mol, n) for n in (*mol.neighbors(bond.a), *mol.neighbors(bond.b)))
        for bond in mol.bonds
    )


def _alert_azide(mol: Molecule) -> bool:
    for atom in mol.atoms:
        if atom.element != "N":
            continue
        nitrogen_neighbours = [n for n in mol.neighbors(atom.index) if mol.atoms[n].element == "N"]
        if len(nitrogen_neighbours) >= 2 and any(
            b.order >= 2 for b in mol.bonds_of(atom.index)
        ):
            return True
    return False


def _alert_nitro(mol: Molecule) -> bool:
    for atom in mol.atoms:
        if atom.element != "N":
            continue
        oxygens = [n for n in mol.neighbors(atom.index) if mol.atoms[n].element == "O"]
        if len(oxygens) >= 2 and all(mol.degree(o) == 1 for o in oxygens):
            return True
    return False


def _alert_thiol(mol: Molecule) -> bool:
    return any(
        a.element == "S" and mol.implicit_hydrogens(a.index) >= 1 and mol.degree(a.index) == 1
        for a in mol.atoms
    )


def _alert_quinone(mol: Molecule) -> bool:
    for ring in ring_systems(mol):
        if len(ring) != 6:
            continue
        exocyclic = sum(
            1 for i in ring
            if mol.atoms[i].element == "C" and any(
                b.order == 2 and mol.atoms[b.other(i)].element == "O" and b.other(i) not in ring
                for b in mol.bonds_of(i)
            )
        )
        if exocyclic >= 2:
            return True
    return False


def _alert_anhydride(mol: Molecule) -> bool:
    for atom in mol.atoms:
        if atom.element != "O" or mol.degree(atom.index) != 2:
            continue
        if all(_has_carbonyl(mol, n) for n in mol.neighbors(atom.index)):
            return True
    return False


def _alert_isocyanate(mol: Molecule) -> bool:
    for atom in mol.atoms:
        if atom.element != "C":
            continue
        doubles = [b for b in mol.bonds_of(atom.index) if b.order == 2]
        elements = sorted(mol.atoms[b.other(atom.index)].element for b in doubles)
        if elements in (["N", "O"], ["N", "S"]):
            return True
    return False


ALERT_RULES: tuple[tuple[str, str, Callable[[Molecule], bool]], ...] = (
    ("acyl_halide", "Hydrolytically unstable and indiscriminately acylating", _alert_acyl_halide),
    ("michael_acceptor", "Covalent electrophile; intentional in some inhibitors", _alert_michael_acceptor),
    ("strained_heterocycle", "Epoxide/aziridine; alkylating and genotoxicity risk", _alert_small_strained_ring),
    ("aldehyde", "Reactive carbonyl; protein adducts and rapid clearance", _alert_aldehyde),
    ("hydrazine", "Hepatotoxicity and mutagenicity liability", _alert_hydrazine),
    ("azide", "Explosive and cytotoxic hazard", _alert_azide),
    ("nitro_group", "Nitroreduction to mutagenic intermediates", _alert_nitro),
    ("free_thiol", "Oxidation, disulfide scrambling, and assay interference", _alert_thiol),
    ("quinone", "Redox cycling and frequent-hitter behaviour", _alert_quinone),
    ("anhydride", "Hydrolytically unstable acylating agent", _alert_anhydride),
    ("isocyanate", "Highly electrophilic; sensitisation risk", _alert_isocyanate),
)


def structural_alerts(molecule: Molecule | str) -> list[dict[str, str]]:
    """Return every reactive or interfering motif matched, with its rationale."""
    mol = _as_molecule(molecule)
    return [
        {"alert": name, "concern": concern}
        for name, concern, rule in ALERT_RULES
        if rule(mol)
    ]


# ---------------------------------------------------------------------------
# Tractability and an overall desirability
# ---------------------------------------------------------------------------

def synthetic_tractability(molecule: Molecule | str, desc: Descriptors | None = None) -> dict[str, Any]:
    """A declared 0-1 proxy for how cheaply a first batch could be made.

    Not the Ertl-Schuffenhauer SA score: there is no fragment-frequency
    corpus behind it.  It is a size-, stereo-, and topology-based penalty,
    which is enough to keep a generator from drifting into 40-heavy-atom
    macrocycles with six stereocentres.
    """
    mol = _as_molecule(molecule)
    desc = desc or compute_descriptors(mol)
    rings = ring_systems(mol)

    size_penalty = max(0.0, (desc.heavy_atoms - 24) / 24.0)
    stereo_penalty = 0.18 * desc.stereocentres
    macrocycle_penalty = 0.35 * sum(1 for r in rings if len(r) > 7)
    unusual_ring_penalty = 0.12 * sum(1 for r in rings if len(r) in (3, 4))

    ring_atom_lists = [set(r) for r in rings]
    fusions = sum(
        1 for i, a in enumerate(ring_atom_lists)
        for b in ring_atom_lists[i + 1:] if len(a & b) >= 2
    )
    spiro = sum(
        1 for i, a in enumerate(ring_atom_lists)
        for b in ring_atom_lists[i + 1:] if len(a & b) == 1
    )
    topology_penalty = 0.08 * fusions + 0.20 * spiro
    exotic_penalty = 0.25 * sum(
        1 for a in mol.heavy_atoms if a.element not in ("C", "N", "O", "S", "F", "Cl", "Br", "I")
    )

    penalty = (size_penalty + stereo_penalty + macrocycle_penalty
               + unusual_ring_penalty + topology_penalty + exotic_penalty)
    score = max(0.0, min(1.0, 1.0 - penalty))
    return {
        "score": round(score, 4),
        "penalties": {
            "size": round(size_penalty, 3),
            "stereocentres": round(stereo_penalty, 3),
            "macrocycle": round(macrocycle_penalty, 3),
            "strained_ring": round(unusual_ring_penalty, 3),
            "ring_topology": round(topology_penalty, 3),
            "exotic_elements": round(exotic_penalty, 3),
        },
    }


# Desirability windows for the overall drug-likeness score. Declared here so
# they can be audited and overridden per campaign rather than buried in code.
DRUGLIKE_WINDOWS: dict[str, tuple[float, float]] = {
    "molecular_weight": (200.0, 450.0),
    "clogp": (0.5, 4.0),
    "tpsa": (40.0, 110.0),
    "hbd": (0.0, 3.0),
    "hba": (1.0, 8.0),
    "rotatable_bonds": (0.0, 8.0),
    "aromatic_rings": (1.0, 3.0),
    "fraction_sp3": (0.25, 1.0),
}


def drug_likeness(desc: Descriptors, windows: dict[str, tuple[float, float]] | None = None) -> dict[str, Any]:
    """A 0-1 desirability score: the geometric mean over declared windows.

    Deliberately *not* called QED -- it does not use Bickerton's fitted
    desirability functions.  The geometric mean is the point: one badly
    out-of-range property drags the whole score down instead of being
    averaged away by seven good ones.
    """
    windows = windows or DRUGLIKE_WINDOWS
    values = desc.as_dict()
    per_property = {}
    for name, (lower, upper) in windows.items():
        per_property[name] = round(_hump(float(values[name]), lower, upper), 4)
    product = 1.0
    for value in per_property.values():
        product *= max(value, 1e-6)
    score = product ** (1.0 / len(per_property))
    return {"score": round(score, 4), "per_property": per_property}


@dataclass(frozen=True)
class DevelopabilityProfile:
    """Everything the design engine knows about a molecule's developability."""

    descriptors: Descriptors
    lipinski: dict[str, Any]
    veber: dict[str, Any]
    cns_mpo: dict[str, Any]
    drug_likeness: dict[str, Any]
    tractability: dict[str, Any]
    alerts: list[dict[str, str]]
    acidic_centres: list[str] = field(default_factory=list)

    @property
    def alert_names(self) -> list[str]:
        return [a["alert"] for a in self.alerts]

    @property
    def anionic_at_physiological_ph(self) -> bool:
        """Whether the molecule carries a group that is ionised at pH 7.4."""
        return bool(self.acidic_centres)

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["descriptors"] = self.descriptors.as_dict()
        return payload


def profile_molecule(molecule: Molecule | str) -> DevelopabilityProfile:
    """Compute the full developability profile for one structure."""
    mol = _as_molecule(molecule)
    desc = compute_descriptors(mol)
    pka, _source = basic_pka_estimate(mol)
    return DevelopabilityProfile(
        descriptors=desc,
        lipinski=lipinski(desc),
        veber=veber(desc),
        cns_mpo=cns_mpo(desc, pka),
        drug_likeness=drug_likeness(desc),
        tractability=synthetic_tractability(mol, desc),
        alerts=structural_alerts(mol),
        acidic_centres=acidic_centres(mol),
    )

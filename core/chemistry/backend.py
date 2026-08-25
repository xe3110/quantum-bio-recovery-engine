"""Chemistry backend selection: RDKit where available, local where not.

The local implementation in this package exists because RDKit has no wheel for
the interpreter this project pins. That is a deployment constraint, not a
judgement about which is better. RDKit is a mature, widely-validated toolkit
and the local stack is a few hundred lines with a measured error bar; where
both are available, RDKit wins.

So the local stack is **fallback infrastructure**, and this module makes that
explicit rather than implicit:

* ``active_backend()`` reports which is in use, and every descriptor vector
  carries that string, so no number is ever read without knowing what computed
  it.
* ``QBRE_CHEMISTRY_BACKEND=local`` forces the fallback, which the test suite
  uses so that assertions do not depend on whether RDKit happens to be
  installed.
* The measured divergence between the two is not a matter of opinion. Run
  ``tools/validate_chemistry.py`` in an RDKit environment; the current numbers
  live in ``docs/chemistry_validation.json``.

The consequence worth stating plainly: **a campaign run with RDKit present and
one run without it can produce different numbers.** They are recorded in the
provenance block of every result file for exactly that reason.
"""

from __future__ import annotations

import os
from typing import Any

ENV_VAR = "QBRE_CHEMISTRY_BACKEND"

try:  # pragma: no cover - environment dependent
    from rdkit import Chem as _Chem
    from rdkit import RDLogger as _RDLogger
    from rdkit.Chem import Crippen as _Crippen
    from rdkit.Chem import Descriptors as _RDDescriptors
    from rdkit.Chem import Lipinski as _Lipinski
    from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors

    _RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except Exception:  # pragma: no cover
    _Chem = None
    HAS_RDKIT = False


def available_backends() -> list[str]:
    return (["rdkit"] if HAS_RDKIT else []) + ["local"]


def active_backend(preference: str | None = None) -> str:
    """Which backend descriptor calls will use.

    Resolution order: explicit argument, then ``QBRE_CHEMISTRY_BACKEND``, then
    RDKit if importable, then local. An explicit request for a backend that is
    not installed raises rather than silently falling back -- a run that
    believes it used RDKit and did not would be worse than a failed run.
    """
    requested = preference or os.environ.get(ENV_VAR)
    if requested:
        requested = requested.strip().lower()
        if requested not in ("rdkit", "local"):
            raise ValueError(
                f"Unknown chemistry backend {requested!r}; expected 'rdkit' or 'local'"
            )
        if requested == "rdkit" and not HAS_RDKIT:
            raise RuntimeError(
                "The rdkit backend was requested but RDKit is not importable in this "
                "interpreter. Install it in an environment that has a wheel (Python 3.11 "
                "here), or unset the preference to use the local fallback."
            )
        return requested
    return "rdkit" if HAS_RDKIT else "local"


def backend_versions() -> dict[str, Any]:
    """Provenance for the result file: what computed these numbers."""
    payload: dict[str, Any] = {
        "active": active_backend(),
        "available": available_backends(),
        "rdkit_available": HAS_RDKIT,
    }
    if HAS_RDKIT:
        payload["rdkit_version"] = _Chem.rdBase.rdkitVersion
    return payload


def rdkit_descriptor_values(smiles: str) -> dict[str, Any]:
    """Compute the shared descriptor set with RDKit.

    HBD and HBA use ``NHOHCount`` and ``NOCount`` -- the Lipinski counting
    convention the local implementation follows -- rather than
    ``NumHDonors``/``NumHAcceptors``, which apply a narrower definition. Mixing
    the two would make the backends disagree for reasons that have nothing to
    do with implementation quality.
    """
    if not HAS_RDKIT:  # pragma: no cover
        raise RuntimeError("RDKit is not available in this interpreter")
    molecule = _Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    return {
        "smiles": _Chem.MolToSmiles(molecule),
        "formula": _rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": round(_RDDescriptors.MolWt(molecule), 3),
        "heavy_atoms": molecule.GetNumHeavyAtoms(),
        "clogp": round(_Crippen.MolLogP(molecule), 3),
        "tpsa": round(_RDDescriptors.TPSA(molecule), 2),
        "hbd": int(_Lipinski.NHOHCount(molecule)),
        "hba": int(_Lipinski.NOCount(molecule)),
        "rotatable_bonds": int(_RDDescriptors.NumRotatableBonds(molecule)),
        "aromatic_rings": int(_rdMolDescriptors.CalcNumAromaticRings(molecule)),
        "rings": int(_rdMolDescriptors.CalcNumRings(molecule)),
        "fraction_sp3": round(float(_rdMolDescriptors.CalcFractionCSP3(molecule)), 4),
        "formal_charge": int(_Chem.GetFormalCharge(molecule)),
        "stereocentres": len(
            _Chem.FindMolChiralCenters(molecule, includeUnassigned=True, useLegacyImplementation=False)
        ),
    }


def rdkit_canonical(smiles: str) -> str:
    if not HAS_RDKIT:  # pragma: no cover
        raise RuntimeError("RDKit is not available in this interpreter")
    molecule = _Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    return _Chem.MolToSmiles(molecule)

"""Backend selection, and agreement between the two implementations.

The local chemistry stack is fallback infrastructure for interpreters that
cannot install RDKit. These tests pin the contract between the two: that the
choice is explicit, recorded, and never silent; and, where RDKit is present,
that the fallback's deviation stays inside the bounds documented in
docs/chemistry_validation.json.
"""

from __future__ import annotations

import pytest

from core.chemistry import backend as chemistry_backend
from core.chemistry.backend import HAS_RDKIT, active_backend, available_backends, backend_versions
from core.chemistry.descriptors import compute_descriptors

CORPUS = [
    "CC(=O)Oc1ccccc1C(=O)O", "CCCCCCCCc1ccc(CCC(N)(CO)CO)cc1",
    "COC(=O)C=CC(=O)OC", "CC12CC3CC(C)(C1)CC(N)(C3)C2",
    "N1(CCN(CC1)c2ccc(Cl)cc2)C(=O)Nc3cccnc3", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
]

requires_rdkit = pytest.mark.skipif(not HAS_RDKIT, reason="RDKit not installed in this interpreter")


def test_local_backend_is_always_available():
    assert "local" in available_backends()


def test_backend_choice_is_recorded_on_every_descriptor_vector():
    """A descriptor without its provenance cannot be interpreted."""
    assert compute_descriptors("CCO").backend == "local"


def test_requesting_an_unknown_backend_raises():
    with pytest.raises(ValueError, match="Unknown chemistry backend"):
        active_backend("openbabel")


def test_requesting_rdkit_when_absent_fails_loudly():
    """A run that believes it used RDKit and silently did not is worse than a failure."""
    if HAS_RDKIT:
        assert active_backend("rdkit") == "rdkit"
    else:
        with pytest.raises(RuntimeError, match="not importable"):
            active_backend("rdkit")


def test_backend_versions_reports_enough_to_reproduce():
    versions = backend_versions()
    assert versions["active"] in ("rdkit", "local")
    assert versions["rdkit_available"] == HAS_RDKIT
    if HAS_RDKIT:
        assert versions["rdkit_version"]


@requires_rdkit
@pytest.mark.parametrize("smiles", CORPUS)
def test_backends_agree_on_composition(smiles):
    """Formula and mass must be identical; these are not estimates."""
    local = compute_descriptors(smiles, backend="local")
    reference = compute_descriptors(smiles, backend="rdkit")
    assert local.formula == reference.formula
    assert local.molecular_weight == pytest.approx(reference.molecular_weight, abs=0.02)


@requires_rdkit
@pytest.mark.parametrize("smiles", CORPUS)
def test_backends_agree_exactly_on_hydrogen_bonding(smiles):
    """HBD and HBA are counting rules, not fitted models -- they must match."""
    local = compute_descriptors(smiles, backend="local")
    reference = compute_descriptors(smiles, backend="rdkit")
    assert local.hbd == reference.hbd
    assert local.hba == reference.hba


@requires_rdkit
@pytest.mark.parametrize("smiles", CORPUS)
def test_local_tpsa_stays_inside_its_documented_error_bar(smiles):
    """TPSA is a lookup over perceived environments, so it inherits perception.

    Exact on almost everything; the exception is a fused heteroaromatic whose
    aromaticity the local model takes as written rather than re-perceiving,
    where the table lookup lands on a different row. Measured worst case is
    3.53 A^2 (docs/chemistry_validation.json).
    """
    local = compute_descriptors(smiles, backend="local")
    reference = compute_descriptors(smiles, backend="rdkit")
    assert abs(local.tpsa - reference.tpsa) <= 4.0


@requires_rdkit
@pytest.mark.parametrize("smiles", CORPUS)
def test_local_clogp_stays_inside_its_documented_error_bar(smiles):
    """The reduced Crippen typing is an approximation with a measured bound.

    docs/chemistry_validation.json records a mean absolute error of 0.65-0.86
    depending on corpus, and a worst case of 2.23, against RDKit's full
    implementation. This test fails if the fallback drifts materially past
    that, which would invalidate the error bar the documentation quotes.

    That bound is wide enough to matter. It is the main reason RDKit is the
    preferred backend and the local logP should not be used to make a
    lipophilicity call on its own.
    """
    local = compute_descriptors(smiles, backend="local")
    reference = compute_descriptors(smiles, backend="rdkit")
    assert abs(local.clogp - reference.clogp) <= 2.5

"""Test-suite fixtures.

The chemistry backend is pinned to the local implementation for the whole
suite. Without this, assertions on descriptor values would pass or fail
depending on whether RDKit happens to be importable in the interpreter running
the tests, which is exactly the kind of environment-dependent result the
project reports provenance to avoid.

``tests/test_chemistry_backend.py`` covers the RDKit path explicitly and skips
itself where RDKit is absent.
"""

from __future__ import annotations

import pytest

from core.chemistry.backend import ENV_VAR


@pytest.fixture(autouse=True, scope="session")
def _pin_local_chemistry_backend() -> None:
    import os

    os.environ[ENV_VAR] = "local"

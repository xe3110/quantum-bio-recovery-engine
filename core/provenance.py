"""Run provenance: what produced a result file, and from exactly which inputs.

A result file is a **snapshot of one run**, not a curated dataset and not a
durable scientific claim. Telling the two apart matters as soon as anything is
regenerated: results under ``experiments/*/results/`` are disposable and
rebuilt on demand, while everything under ``data/`` is curated, reviewed, and
the thing a reviewer is expected to argue with.

What makes a snapshot interpretable later is knowing precisely what went into
it. This module records:

* a **content digest of every input file**, so a result can be matched to the
  exact signature, panel, interactome, and library that produced it -- a
  version string in a metadata block does not survive someone editing a CSV;
* the **environment**, including the chemistry backend, because RDKit being
  importable changes descriptor values and therefore rankings;
* the **git revision and whether the tree was dirty**, since a clean revision
  hash on a modified tree is worse than no hash at all.
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parents[1]

# Packages whose version can change a result. Recorded by import rather than by
# parsing requirements.txt, so the record reflects what actually ran.
TRACKED_PACKAGES = (
    "numpy", "scipy", "pandas", "networkx", "rdkit",
    "qiskit", "qiskit_aer", "qiskit_algorithms", "qiskit_optimization",
)


def file_digest(path: Path | str, algorithm: str = "sha256") -> dict[str, Any]:
    """Content digest and size for one input file."""
    path = Path(path)
    digest = hashlib.new(algorithm)
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
            size += len(chunk)
    try:
        name = str(path.resolve().relative_to(ROOT))
    except ValueError:
        name = str(path.resolve())
    return {"path": name, algorithm: digest.hexdigest(), "bytes": size}


def input_manifest(paths: Iterable[Path | str]) -> list[dict[str, Any]]:
    """Digest every input, skipping duplicates and reporting missing files.

    A missing input is recorded rather than raising: the run has already
    happened by the time provenance is written, and losing the whole record
    over one unreadable path would be the wrong trade.
    """
    seen: set[str] = set()
    manifest: list[dict[str, Any]] = []
    for path in paths:
        if path is None:
            continue
        resolved = Path(path).resolve()
        if str(resolved) in seen:
            continue
        seen.add(str(resolved))
        if not resolved.exists():
            manifest.append({"path": str(path), "error": "missing at write time"})
            continue
        manifest.append(file_digest(resolved))
    return sorted(manifest, key=lambda entry: entry["path"])


def package_versions() -> dict[str, str]:
    """Installed versions of the packages that can change a result."""
    from importlib.metadata import PackageNotFoundError, version

    found: dict[str, str] = {}
    for name in TRACKED_PACKAGES:
        try:
            found[name] = version(name.replace("_", "-"))
        except PackageNotFoundError:
            found[name] = "not installed"
        except Exception:
            found[name] = "unknown"
    return found


def git_state() -> dict[str, Any]:
    """Revision, branch, and whether the working tree had uncommitted changes.

    ``dirty`` is the field that matters. A revision hash recorded from a
    modified tree describes code that was never committed, which is worse than
    recording nothing, so the flag travels with the hash.
    """
    def run(*args: str) -> str | None:
        try:
            return subprocess.check_output(
                args, cwd=ROOT, stderr=subprocess.DEVNULL
            ).decode().strip()
        except Exception:
            return None

    status = run("git", "status", "--porcelain")
    return {
        "revision": run("git", "rev-parse", "HEAD") or "unknown",
        "branch": run("git", "rev-parse", "--abbrev-ref", "HEAD") or "unknown",
        "dirty": bool(status) if status is not None else None,
        "uncommitted_files": len(status.splitlines()) if status else 0,
    }


def environment() -> dict[str, Any]:
    """Interpreter, platform, packages, and the active chemistry backend."""
    from core.chemistry.backend import backend_versions

    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "executable": sys.executable,
        "packages": package_versions(),
        "chemistry_backend": backend_versions(),
    }


def run_provenance(
    inputs: Iterable[Path | str],
    command: Iterable[str] | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """The full provenance block for one run artifact."""
    payload = {
        "artifact_class": "run_snapshot",
        "artifact_note": (
            "A snapshot of one run, regenerated by re-running the command below. It is not "
            "a curated dataset: everything under data/ is curated and reviewed, everything "
            "under experiments/*/results/ is disposable output. Do not cite a number from "
            "this file without the input digests that accompany it."
        ),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(command) if command else " ".join(sys.argv),
        "git": git_state(),
        "environment": environment(),
        "inputs": input_manifest(inputs),
    }
    if extra:
        payload.update(extra)
    return payload

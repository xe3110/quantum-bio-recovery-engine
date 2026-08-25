"""Cross-validate the local chemistry model against RDKit.

Why this exists
---------------
`core/chemistry/` is a hand-written cheminformatics stack, present because
RDKit has no wheel for the interpreter this project pins. Formula and mass
checks confirm *composition*; they say nothing about whether the model
perceives rings, aromaticity, donors, or lipophilicity the way a reference
toolkit does. Those are the quantities decisions are actually made on.

This harness runs both implementations over the same corpus and reports the
deviation per descriptor, including the molecules the pipeline *generates*
rather than only the ones it was curated from. It is the evidence for -- or
against -- treating the local stack as usable.

Structural perception is checked the strictest way available: the local
writer's output is re-read by RDKit and canonicalised, then compared to
RDKit's canonicalisation of the original input. Agreement means the local
parser and writer together preserve the molecule as RDKit understands it.

Run in an environment that has RDKit (Python 3.11 here):

    PYTHONPATH=. qbio-quantum-env/bin/python -m tools.validate_chemistry
    PYTHONPATH=. qbio-quantum-env/bin/python -m tools.validate_chemistry --designs 400
"""

from __future__ import annotations

import argparse
import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs/chemistry_validation.json"

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
except ImportError as error:  # pragma: no cover
    raise SystemExit(
        "RDKit is required for this harness and is not installed in this interpreter.\n"
        "It has no wheel for the project's pinned Python; use the 3.11 environment:\n"
        "    qbio-quantum-env/bin/pip install rdkit\n"
        "    PYTHONPATH=. qbio-quantum-env/bin/python -m tools.validate_chemistry"
    ) from error

from core.chemistry import descriptors as local
from core.chemistry.molecule import parse_smiles, write_smiles


# Each entry: name, local callable, RDKit callable, tolerance, integer?
#
# HBD/HBA are compared against RDKit's NHOHCount and NOCount rather than
# NumHDonors/NumHAcceptors, because the local implementation uses the Lipinski
# counting convention (hydrogens on N/O; count of N and O) and RDKit's
# NumHDonors/NumHAcceptors apply a different, narrower definition.
DESCRIPTORS: tuple[tuple[str, Callable, Callable, float, bool], ...] = (
    ("molecular_weight", lambda m: m.molecular_weight(), Descriptors.MolWt, 0.05, False),
    ("clogp", local.clogp, Crippen.MolLogP, 0.0, False),
    ("tpsa", local.tpsa, Descriptors.TPSA, 0.0, False),
    ("hbd", local.hbd_count, Lipinski.NHOHCount, 0.0, True),
    ("hba", local.hba_count, Lipinski.NOCount, 0.0, True),
    ("rotatable_bonds", local.rotatable_bonds, Descriptors.NumRotatableBonds, 0.0, True),
    ("rings", lambda m: len(local.ring_systems(m)), rdMolDescriptors.CalcNumRings, 0.0, True),
    ("aromatic_rings", lambda m: len(local.aromatic_rings(m)),
     rdMolDescriptors.CalcNumAromaticRings, 0.0, True),
    ("fraction_sp3", local.fraction_sp3, rdMolDescriptors.CalcFractionCSP3, 0.0, False),
)


@dataclass
class Deviation:
    """Per-descriptor agreement between the two implementations."""

    name: str
    n: int = 0
    exact: int = 0
    within_tolerance: int = 0
    errors: list[float] = field(default_factory=list)
    worst: tuple[float, str] = (0.0, "")

    def add(self, local_value: float, reference: float, smiles: str, tolerance: float) -> None:
        error = float(local_value) - float(reference)
        self.n += 1
        self.errors.append(error)
        if abs(error) <= 1e-9:
            self.exact += 1
        if abs(error) <= tolerance + 1e-9:
            self.within_tolerance += 1
        if abs(error) > abs(self.worst[0]):
            self.worst = (round(error, 4), smiles)

    def summary(self) -> dict[str, Any]:
        absolute = [abs(e) for e in self.errors]
        return {
            "n": self.n,
            "exact_agreement": round(self.exact / self.n, 4) if self.n else None,
            "mean_absolute_error": round(statistics.fmean(absolute), 4) if absolute else None,
            "median_absolute_error": round(statistics.median(absolute), 4) if absolute else None,
            "p95_absolute_error": round(
                sorted(absolute)[int(0.95 * (len(absolute) - 1))], 4
            ) if absolute else None,
            "max_absolute_error": round(max(absolute), 4) if absolute else None,
            "bias_mean_signed_error": round(statistics.fmean(self.errors), 4) if self.errors else None,
            "worst_case": {"error": self.worst[0], "smiles": self.worst[1]},
        }


def canonical(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(molecule) if molecule is not None else None


@dataclass
class PerceptionResult:
    """Whether the local parser+writer preserve the molecule RDKit sees."""

    n: int = 0
    rdkit_accepts_local_output: int = 0
    canonical_match: int = 0
    mismatches: list[dict[str, str]] = field(default_factory=list)
    rejections: list[str] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "rdkit_parses_local_output": self.rdkit_accepts_local_output,
            "rdkit_rejection_rate": round(1 - self.rdkit_accepts_local_output / self.n, 4) if self.n else None,
            "canonical_smiles_match": self.canonical_match,
            "canonical_match_rate": round(self.canonical_match / self.n, 4) if self.n else None,
            "example_mismatches": self.mismatches[:15],
            "example_rejections": self.rejections[:15],
        }


def build_corpus(n_designs: int) -> dict[str, list[tuple[str, str]]]:
    """Assemble the validation corpus, grouped by provenance."""
    from core.chemistry.molecule import cap_attachments
    from core.design.pharmacophores import load_library

    corpus: dict[str, list[tuple[str, str]]] = {}

    known = json.loads((ROOT / "data/chemistry/ms_known_structures.json").read_text())
    corpus["curated_drugs"] = [
        (e["name"], e["smiles"]) for e in known["structures"] if e["smiles"]
    ]

    library = load_library()
    corpus["library_fragments"] = []
    for fragment in library.fragments:
        try:
            capped = cap_attachments(fragment.molecule())
        except Exception:
            continue
        # A one-atom linker or cap collapses to a single heavy atom when its
        # attachment points are removed -- an ether linker becomes water. Those
        # are artefacts of capping in isolation, never structures the pipeline
        # emits, and comparing descriptors on them measures nothing.
        if len(capped.heavy_atoms) < 2:
            continue
        corpus["library_fragments"].append((fragment.identifier, write_smiles(capped)))

    corpus["generated_designs"] = _sample_designs(n_designs)
    return corpus


def _sample_designs(limit: int) -> list[tuple[str, str]]:
    """Assemble real designs -- the molecules the pipeline actually emits.

    Curated inputs are hand-checked; generated structures are not, and they
    are the ones exercising unusual ring fusions and substitution patterns.
    Validating only the curated set would test the easy half of the corpus.
    """
    from itertools import product

    from core.chemistry.molecule import SmilesError
    from core.design.denovo import Recipe, assemble
    from core.design.pharmacophores import load_library

    library = load_library()
    arms = [f.identifier for f in library.pharmacophores]
    scaffolds = [s.identifier for s in library.scaffolds]
    linkers = [l.identifier for l in library.linkers]
    caps = [c.identifier for c in library.caps]

    seen: set[str] = set()
    out: list[tuple[str, str]] = []
    for index, (scaffold, arm_a, arm_b, linker_a, linker_b) in enumerate(
        product(scaffolds, arms, arms, linkers, linkers)
    ):
        if arm_a == arm_b:
            continue
        spare = library.get(scaffold).attachment_count - 2
        cap_choice = (caps[index % len(caps)],) * spare
        try:
            molecule = assemble(
                library, Recipe(scaffold, (arm_a, arm_b), (linker_a, linker_b), cap_choice)
            )
            smiles = molecule.to_smiles()
        except (SmilesError, KeyError, ValueError):
            continue
        if smiles in seen:
            continue
        seen.add(smiles)
        out.append((f"design_{len(out):04d}", smiles))
        if len(out) >= limit:
            break
    return out


def validate(corpus: dict[str, list[tuple[str, str]]]) -> dict[str, Any]:
    results: dict[str, Any] = {}
    for group, entries in corpus.items():
        perception = PerceptionResult()
        deviations = {name: Deviation(name) for name, *_ in DESCRIPTORS}
        local_failures: list[dict[str, str]] = []

        for name, smiles in entries:
            reference_molecule = Chem.MolFromSmiles(smiles)
            if reference_molecule is None:
                continue  # RDKit cannot read the input; not a local failure
            reference_canonical = Chem.MolToSmiles(reference_molecule)

            try:
                molecule = parse_smiles(smiles)
                written = write_smiles(molecule)
            except Exception as error:
                local_failures.append({"name": name, "smiles": smiles, "error": str(error)[:120]})
                continue

            perception.n += 1
            round_tripped = canonical(written)
            if round_tripped is not None:
                perception.rdkit_accepts_local_output += 1
                if round_tripped == reference_canonical:
                    perception.canonical_match += 1
                else:
                    perception.mismatches.append({
                        "input": smiles, "local_output": written,
                        "rdkit_of_input": reference_canonical,
                        "rdkit_of_local_output": round_tripped,
                    })
            else:
                perception.rejections.append(written)

            for descriptor, local_fn, reference_fn, tolerance, _integer in DESCRIPTORS:
                try:
                    deviations[descriptor].add(
                        local_fn(molecule), reference_fn(reference_molecule), smiles, tolerance
                    )
                except Exception:
                    continue

        results[group] = {
            "n_entries": len(entries),
            "local_parse_failures": local_failures,
            "structural_perception": perception.summary(),
            "descriptors": {k: v.summary() for k, v in deviations.items() if v.n},
        }
    return results


def validate_fingerprints(corpus: dict[str, list[tuple[str, str]]], sample: int = 120) -> dict[str, Any]:
    """Compare local circular fingerprints against RDKit Morgan fingerprints.

    Absolute Tanimoto values are not expected to agree: the local
    implementation keeps an unfolded identifier set while RDKit folds into a
    bit vector, and the invariants differ. What matters for novelty is whether
    the two *rank* similarity the same way, so this reports rank correlation
    alongside the raw agreement.
    """
    from rdkit.Chem import rdFingerprintGenerator

    from core.chemistry.fingerprint import circular_fingerprint, tanimoto

    entries = [e for group in corpus.values() for e in group][:sample]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    local_scores: list[float] = []
    reference_scores: list[float] = []
    for i in range(len(entries) - 1):
        for j in range(i + 1, min(i + 6, len(entries))):
            a, b = entries[i][1], entries[j][1]
            mol_a, mol_b = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
            if mol_a is None or mol_b is None:
                continue
            try:
                local_scores.append(tanimoto(circular_fingerprint(a), circular_fingerprint(b)))
            except Exception:
                continue
            from rdkit import DataStructs

            reference_scores.append(DataStructs.TanimotoSimilarity(
                generator.GetFingerprint(mol_a), generator.GetFingerprint(mol_b)
            ))

    if len(local_scores) < 3:
        return {"n_pairs": len(local_scores), "note": "insufficient pairs"}

    from scipy.stats import pearsonr, spearmanr

    spearman = spearmanr(local_scores, reference_scores)
    pearson = pearsonr(local_scores, reference_scores)
    differences = [abs(a - b) for a, b in zip(local_scores, reference_scores)]
    # The 0.6 line is what assess_novelty() uses to call a design "novel".
    disagreements = sum(
        1 for a, b in zip(local_scores, reference_scores) if (a >= 0.6) != (b >= 0.6)
    )
    return {
        "n_pairs": len(local_scores),
        "spearman_rho": round(float(spearman.statistic), 4),
        "pearson_r": round(float(pearson.statistic), 4),
        "mean_absolute_difference": round(statistics.fmean(differences), 4),
        "max_absolute_difference": round(max(differences), 4),
        "novelty_threshold_disagreements": disagreements,
        "novelty_threshold_disagreement_rate": round(disagreements / len(local_scores), 4),
        "interpretation": (
            "Absolute Tanimoto values are not expected to match: the local fingerprint keeps "
            "an unfolded identifier set, RDKit folds to 2048 bits, and the atom invariants "
            "differ. Rank correlation and the rate of disagreement at the 0.6 novelty "
            "threshold are the quantities that matter for how the fingerprint is used."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--designs", type=int, default=300,
                        help="generated structures to include in the corpus")
    parser.add_argument("--out", type=Path, default=REPORT)
    args = parser.parse_args()

    corpus = build_corpus(args.designs)
    print("Corpus:")
    for group, entries in corpus.items():
        print(f"  {group:20} {len(entries)}")

    results = validate(corpus)
    fingerprints = validate_fingerprints(corpus)

    for group, outcome in results.items():
        perception = outcome["structural_perception"]
        print(f"\n{group}  (n={perception['n']})")
        print(f"  RDKit parses local output : {perception['rdkit_parses_local_output']}/{perception['n']}")
        print(f"  canonical SMILES match    : {perception['canonical_smiles_match']}/{perception['n']}"
              f"  ({perception['canonical_match_rate']})")
        if outcome["local_parse_failures"]:
            print(f"  local parse failures      : {len(outcome['local_parse_failures'])}")
        print(f"  {'descriptor':18}{'exact':>8}{'MAE':>9}{'p95':>9}{'max':>9}{'bias':>9}")
        for name, stats in outcome["descriptors"].items():
            print(f"  {name:18}{stats['exact_agreement']:>8.3f}{stats['mean_absolute_error']:>9.3f}"
                  f"{stats['p95_absolute_error']:>9.3f}{stats['max_absolute_error']:>9.3f}"
                  f"{stats['bias_mean_signed_error']:>9.3f}")

    print(f"\nfingerprints: spearman={fingerprints.get('spearman_rho')} "
          f"MAD={fingerprints.get('mean_absolute_difference')} "
          f"threshold disagreements={fingerprints.get('novelty_threshold_disagreements')}"
          f"/{fingerprints.get('n_pairs')}")

    report = {
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "descriptor_comparison": results,
        "fingerprint_comparison": fingerprints,
        "method": {
            "structural_perception": (
                "The local writer's output is re-read by RDKit and canonicalised, then compared "
                "against RDKit's canonicalisation of the original input. Agreement means the "
                "local parser and writer together preserve the molecule as RDKit perceives it."
            ),
            "hbd_hba_convention": (
                "Compared against RDKit NHOHCount and NOCount, matching the Lipinski counting "
                "convention the local implementation uses, not NumHDonors/NumHAcceptors."
            ),
        },
    }
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nWrote {args.out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

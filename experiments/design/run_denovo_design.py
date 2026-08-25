"""Design a new molecule for a disease, end to end.

Pipeline
--------
1. Load the disease context: signature, interactome, panel, druggability.
2. Derive the Target Product Profile -- which proteins are worth engaging,
   in which direction, inside which physicochemical envelope.
3. Formulate pharmacophore selection as a QUBO and solve it, optionally
   benchmarking QAOA and an exact eigensolver against exhaustive enumeration.
4. Assemble, refine, and score real molecular structures for the winning arm
   sets, then filter for novelty and structural diversity.

Everything produced is a discovery-stage hypothesis.  No structure here has
been docked, simulated, synthesised, or assayed, and no number in the output
is a measured affinity.  See docs/denovo_design_protocol.md.

Usage:
    python -m experiments.design.run_denovo_design --disease multiple_sclerosis
    python -m experiments.design.run_denovo_design --quantum-benchmark --arms 3
    python -m experiments.design.run_denovo_design --quick
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from core.design.denovo import DesignWeights, run_design
from core.design.pharmacophores import load_library, unreachable_requirements
from core.design.quantum_assembly import (
    HAS_QISKIT, SelectionWeights, benchmark_backends, build_selection_problem,
    solve_enumeration,
)
from core.design.target_profile import build_target_profile
from core.models.disease import available_diseases, load_disease
from core.provenance import run_provenance

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "experiments/design/results"

CSV_COLUMNS = [
    "rank", "smiles", "molecular_formula", "fitness", "profile_coverage",
    "reversal_efficiency", "counter_therapeutic", "molecular_weight", "clogp", "tpsa",
    "hbd", "hba", "rotatable_bonds", "aromatic_rings", "fraction_sp3", "heavy_atoms",
    "cns_mpo", "meets_delivery_gate", "drug_likeness", "synthetic_tractability",
    "ligand_efficiency", "lipinski_violations", "structural_alerts",
    "nearest_known_compound", "max_tanimoto", "is_novel",
    "scaffold", "arms", "linkers", "caps", "window_violations",
]


def _display(path: Path) -> str:
    """Repo-relative where possible, absolute otherwise."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


def _csv_row(rank: int, candidate: dict) -> dict:
    descriptors = candidate["descriptors"]
    novelty = candidate["novelty"] or {}
    recipe = candidate["recipe"]
    return {
        "rank": rank,
        "smiles": candidate["smiles"],
        "molecular_formula": candidate["molecular_formula"],
        "fitness": candidate["fitness"],
        "profile_coverage": candidate["profile_coverage"],
        "reversal_efficiency": candidate["signature_metrics"].get("reversal_efficiency"),
        "counter_therapeutic": candidate["signature_metrics"].get("counter_therapeutic"),
        "molecular_weight": descriptors["molecular_weight"],
        "clogp": descriptors["clogp"],
        "tpsa": descriptors["tpsa"],
        "hbd": descriptors["hbd"],
        "hba": descriptors["hba"],
        "rotatable_bonds": descriptors["rotatable_bonds"],
        "aromatic_rings": descriptors["aromatic_rings"],
        "fraction_sp3": descriptors["fraction_sp3"],
        "heavy_atoms": descriptors["heavy_atoms"],
        "cns_mpo": candidate["cns_mpo"]["total"],
        "meets_delivery_gate": candidate["meets_delivery_gate"],
        "drug_likeness": candidate["drug_likeness"]["score"],
        "synthetic_tractability": candidate["synthetic_tractability"]["score"],
        "ligand_efficiency": candidate["ligand_efficiency"],
        "lipinski_violations": candidate["lipinski"]["violations"],
        "structural_alerts": ";".join(a["alert"] for a in candidate["structural_alerts"]),
        "nearest_known_compound": novelty.get("nearest_known_compound", ""),
        "max_tanimoto": novelty.get("max_tanimoto", ""),
        "is_novel": novelty.get("is_novel", ""),
        "scaffold": recipe["scaffold"],
        "arms": ";".join(recipe["arms"]),
        "linkers": ";".join(recipe["linkers"]),
        "caps": ";".join(recipe["caps"]),
        "window_violations": ";".join(candidate["window_violations"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--disease", default="multiple_sclerosis",
                        help=f"registry entry to design against (available: {', '.join(available_diseases())})")
    parser.add_argument("--targets", type=int, default=14, help="size of the target product profile")
    parser.add_argument("--k", type=int, default=3, help="pharmacophore arms per molecule")
    parser.add_argument("--arms", type=int, default=3, help="how many arm sets to carry into assembly")
    parser.add_argument("--qubits", type=int, default=10,
                        help="fragments in the Hamiltonian; one binary variable each")
    parser.add_argument("--budget", type=int, default=None,
                        help="heavy-atom budget for the arms (default: derived from the "
                             "profile's molecular-weight ceiling, less the scaffold and linkers)")
    parser.add_argument("--top", type=int, default=10, help="designs to report")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--iterations", type=int, default=400, help="refinement steps per restart")
    parser.add_argument("--quantum-benchmark", action="store_true",
                        help="also solve the Hamiltonian with QAOA and an exact eigensolver")
    parser.add_argument("--qaoa-reps", type=int, nargs="+", default=[1, 2, 3, 4, 5],
                        help="QAOA circuit depths to sweep; the best feasible result is kept")
    parser.add_argument("--quick", action="store_true", help="fast smoke run")
    parser.add_argument("--outdir", type=Path, default=RESULTS)
    args = parser.parse_args()

    if args.quick:
        args.targets, args.k, args.arms, args.qubits, args.top, args.iterations = 8, 2, 1, 8, 3, 40

    disease = load_disease(args.disease)
    library = load_library()

    print(f"Disease      : {disease.name}")
    print(f"Signature    : {len(disease.signature().genes)} directional genes")
    print(f"Library      : {len(library)} fragments "
          f"({len(library.pharmacophores)} pharmacophores)")

    profile = build_target_profile(disease, top_n=args.targets)
    unreachable = unreachable_requirements(library, profile.genes)
    print(f"\nTarget profile ({len(profile.requirements)} requirements):")
    for requirement in profile.requirements:
        arrow = "up  " if requirement.desired_direction > 0 else "down"
        print(f"  {requirement.gene:9} {arrow} priority={requirement.priority:.4f} "
              f"[{requirement.target_class}]")
    if unreachable:
        print(f"  no chemical matter in the library for: {', '.join(unreachable)}")

    problem = build_selection_problem(
        profile, library, k=args.k, heavy_atom_budget=args.budget, max_variables=args.qubits
    )
    summary = problem.hamiltonian_summary()
    print(f"\nHamiltonian  : {summary['n_variables']} binary variables, "
          f"{summary['n_quadratic_terms']} couplings, choose {args.k} "
          f"({summary['feasible_states']} feasible of {summary['search_space']} states)")

    quantum_report = None
    if args.quantum_benchmark:
        if not HAS_QISKIT:
            print("  qiskit is not installed in this interpreter; skipping the quantum backends.")
        else:
            quantum_report = benchmark_backends(
                problem, reps_range=tuple(args.qaoa_reps), seed=args.seed
            )
            for backend, outcome in quantum_report["backends"].items():
                if "error" in outcome or "skipped" in outcome:
                    print(f"  {backend:12} {outcome.get('error') or outcome.get('skipped')}")
                    continue
                flag = "match" if outcome["matches_enumeration_optimum"] else "DIFFERS"
                print(f"  {backend:12} {flag:8} objective={outcome['qubo_objective']:+.4f} "
                      f"feasible={outcome['feasible']} {outcome['seconds']}s")
                if outcome.get("note"):
                    print(f"               {outcome['note']}")
                sweep = outcome.get("detail", {}).get("depth_sweep")
                if sweep:
                    summary = "  ".join(
                        f"d{a['reps']}={'-' if a.get('error') else format(a['qubo_objective'], '+.4f')}"
                        f"{'' if a.get('feasible', True) else '!'}"
                        for a in sweep
                    )
                    print(f"               depth sweep: {summary}   (! = constraint violated)")

    selections = solve_enumeration(problem, top=args.arms)
    print(f"\nArm sets carried into assembly:")
    for index, selection in enumerate(selections, 1):
        print(f"  {index}. objective={selection.qubo_objective:+.4f}  "
              f"{' + '.join(selection.identifiers)}")

    print("\nAssembling...")
    design = run_design(
        disease=disease,
        profile=profile,
        library=library,
        arm_sets=[s.identifiers for s in selections],
        weights=DesignWeights(),
        seed=args.seed,
        top=args.top,
        iterations=args.iterations,
    )
    search = design["search"]
    print(f"  {search['recipes_attempted']} recipes attempted, "
          f"{search['assembly_failures']} rejected as chemically invalid, "
          f"{search['unique_structures']} distinct structures")

    print(f"\nTop {len(design['candidates'])} designs:")
    for index, candidate in enumerate(design["candidates"], 1):
        novelty = candidate["novelty"] or {}
        gate = "" if candidate["meets_delivery_gate"] else "  [below CNS gate]"
        print(f"\n  {index}. {candidate['molecular_formula']}  fitness={candidate['fitness']:.4f}{gate}")
        print(f"     {candidate['smiles']}")
        print(f"     MW={candidate['descriptors']['molecular_weight']:.1f} "
              f"cLogP={candidate['descriptors']['clogp']:.2f} "
              f"TPSA={candidate['descriptors']['tpsa']:.1f} "
              f"CNS-MPO={candidate['cns_mpo']['total']:.2f} "
              f"coverage={candidate['profile_coverage']:.3f}")
        print(f"     arms: {' + '.join(candidate['recipe']['arms'])}")
        print(f"     nearest known: {novelty.get('nearest_known_compound','?')} "
              f"(Tanimoto {novelty.get('max_tanimoto','?')}, "
              f"novel={novelty.get('is_novel','?')})")
        if candidate["structural_alerts"]:
            print(f"     alerts: {', '.join(a['alert'] for a in candidate['structural_alerts'])}")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    report = {
        "disease": disease.as_dict(),
        "target_product_profile": profile.as_dict(),
        "unreachable_requirements": unreachable,
        "hamiltonian": summary,
        "selection_weights": problem.weights.as_dict(),
        "arm_sets": [s.as_dict() for s in selections],
        "quantum_benchmark": quantum_report,
        "design": design,
        "library_metadata": library.metadata,
        "provenance": run_provenance(
            inputs=[
                disease.signature_path, disease.network_path, disease.panel_path,
                disease.druggability_path, disease.structures_path,
                ROOT / "data/chemistry/pharmacophore_library.json",
                ROOT / f"data/diseases/{disease.identifier}.json",
            ],
            extra={"qiskit_available": HAS_QISKIT, "seed": args.seed},
        ),
        "interpretation": (
            "Discovery-stage structure hypotheses. Target engagement is inherited from curated "
            "fragment annotations, not predicted from structure; no docking, free-energy, or "
            "ADMET model was run; no compound was synthesised or assayed. Read the target "
            "product profile and the mechanism combinations, not the ranking of individual "
            "molecules."
        ),
    }
    json_path = outdir / "denovo_designs.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n")

    csv_path = outdir / "denovo_candidates.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for index, candidate in enumerate(design["all_ranked"], 1):
            writer.writerow(_csv_row(index, candidate))

    print(f"\nWrote {_display(json_path)}")
    print(f"Wrote {_display(csv_path)}")


if __name__ == "__main__":
    main()

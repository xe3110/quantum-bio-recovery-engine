"""Run a deterministic MS combination-candidate prioritisation analysis.

Usage: python -m experiments.ms.run_publishable_prioritization --top 25
"""
import argparse
import csv
import json
from pathlib import Path

from core.biology.ms_prioritization import PrioritizationConfig, load_signature, rank_pairs, sensitivity_interval

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser(description="Research-use-only MS combination candidate prioritisation")
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--draws", type=int, default=500)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--include-preclinical", action="store_true")
    args = parser.parse_args()

    with (ROOT / "data/drugs/ms_panel_v2.json").open() as handle:
        panel = json.load(handle)
    with (ROOT / "data/ms_expression_v2.csv").open(newline="") as handle:
        signature = load_signature(csv.DictReader(handle))
    config = PrioritizationConfig(minimum_evidence="preclinical" if args.include_preclinical else "phase_2")
    rows = rank_pairs(panel["drugs"], signature, config)
    by_name = {d["name"]: d for d in panel["drugs"]}
    primary = [row for row in rows if not row["excluded_from_primary_ranking"]][: args.top]
    for row in primary:
        row["sensitivity_95ci"] = sensitivity_interval(by_name[row["drug_a"]], by_name[row["drug_b"]], signature, args.draws, args.seed)

    output = ROOT / "experiments/ms/results"
    output.mkdir(exist_ok=True)
    with (output / "ms_candidate_pairs.json").open("w") as handle:
        json.dump({"metadata": panel["metadata"], "config": config.__dict__, "n_candidates": len(rows), "primary_results": primary}, handle, indent=2)
    print(f"Screened {len(panel['drugs'])} candidates and {len(rows)} pairs; wrote {len(primary)} eligible pairs to {output / 'ms_candidate_pairs.json'}.")


if __name__ == "__main__":
    main()

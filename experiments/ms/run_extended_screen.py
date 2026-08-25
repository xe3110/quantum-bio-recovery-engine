"""Run the extended, multi-parameter MS combination-candidate screen.

Discovery-stage hypothesis generation only.  Nothing produced here is a
clinical combination recommendation; see docs/ms_publication_protocol.md.

Usage:
    python -m experiments.ms.run_extended_screen --top 30 --seed 7
    python -m experiments.ms.run_extended_screen --quick    # fast smoke run
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np

from core.biology.ms_scoring import (
    RISK_DOMAINS, ScoringConfig, eligible_drugs, load_signature, rank_pairs,
)
from core.biology.network_proximity import ProximityScorer, classify_exposure, load_network
from core.biology.screen_statistics import (
    axis_stratum_enrichment, benjamini_hochberg, bootstrap_stability, pareto_front,
    permutation_pvalues, weight_sensitivity,
)

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "experiments/ms/results"

CSV_COLUMNS = [
    "rank", "drug_a", "drug_b", "priority_score", "p_empirical", "q_value", "permutation_z",
    "signed_reversal", "reversal_efficiency", "reversal_gain_over_best_single",
    "counter_therapeutic", "gene_coverage", "pathway_coverage", "axis_coverage",
    "target_complementarity", "compartment_complementarity", "cns_reach", "regimen_compatibility",
    "network_separation", "proximity_z_a", "proximity_z_b", "complementary_exposure",
    "evidence_score", "uncertainty", "safety_union", "safety_overlap", "worst_risk_domain",
    "same_mechanism", "shared_safety_class", "redundant_targets",
    "excluded_from_primary_ranking", "pareto_optimal",
]


def _git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


def build_exposures(drugs, signature_genes, seed, permutations):
    """Pre-compute the interactome exposure classification for every pair."""
    scorer = ProximityScorer(load_network(), seed=seed, permutations=permutations)
    exposures = {}
    for a, b in combinations(drugs, 2):
        key = (min(a["name"], b["name"]), max(a["name"], b["name"]))
        exposures[key] = classify_exposure(
            scorer, list(a["target_effects"]), list(b["target_effects"]), signature_genes
        )
    return exposures


def control_report(rows, controls, top_k):
    """Check the screen against pre-declared positive and negative controls.

    Declared before ranking, this is the screen's own falsification test: a
    redundant pair ranking highly, or a withdrawn agent surfacing in the top
    list, indicates the scoring is wrong rather than the biology surprising.
    """
    primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
    position = {(r["drug_a"], r["drug_b"]): i + 1 for i, r in enumerate(primary)}
    top_names = {n for r in primary[:top_k] for n in (r["drug_a"], r["drug_b"])}

    def best_rank(name):
        hits = [p for (a, b), p in position.items() if name in (a, b)]
        return min(hits) if hits else None

    redundant = controls.get("positive_redundancy", [])
    redundant_excluded = None
    if len(redundant) >= 2:
        key = (min(redundant[0], redundant[1]), max(redundant[0], redundant[1]))
        match = next((r for r in rows if (min(r["drug_a"], r["drug_b"]), max(r["drug_a"], r["drug_b"])) == key), None)
        redundant_excluded = bool(match["excluded_from_primary_ranking"]) if match else None

    return {
        "n_primary_pairs": len(primary),
        "redundant_pair": {
            "pair": redundant,
            "excluded_as_expected": redundant_excluded,
            "interpretation": "A mechanistically duplicate pair must be excluded from the primary ranking.",
        },
        "safety_controls": {
            name: {"best_primary_rank": best_rank(name), "in_top_k": name in top_names}
            for name in controls.get("safety_penalty", [])
        },
        "negative_efficacy_controls": {
            name: {"best_primary_rank": best_rank(name), "in_top_k": name in top_names}
            for name in controls.get("negative_efficacy", [])
        },
        "note": (
            "Negative-efficacy controls are agents whose MS trials failed. Because the "
            "screen scores mechanism and transcriptional direction rather than trial "
            "outcome, some may still rank highly; that is a stated limitation of the "
            "method, and any such appearance is reported rather than filtered away."
        ),
    }


def diversity_filter(rows, top_k, max_per_drug):
    """Greedy top-K selection capping how often any single agent may appear.

    One broadly-acting agent can otherwise occupy most of the leaderboard,
    which hides the mechanistic diversity the screen exists to surface.
    """
    counts, chosen = {}, []
    for row in rows:
        if row["excluded_from_primary_ranking"]:
            continue
        a, b = row["drug_a"], row["drug_b"]
        if counts.get(a, 0) >= max_per_drug or counts.get(b, 0) >= max_per_drug:
            continue
        counts[a] = counts.get(a, 0) + 1
        counts[b] = counts.get(b, 0) + 1
        chosen.append(row)
        if len(chosen) >= top_k:
            break
    return chosen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--permutations", type=int, default=400, help="draws for the significance null")
    parser.add_argument("--proximity-permutations", type=int, default=200)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--weight-trials", type=int, default=40)
    parser.add_argument("--max-per-drug", type=int, default=3, help="cap on appearances per agent in the diverse top list")
    parser.add_argument("--min-evidence", default="phase_2",
                        choices=["approved", "phase_3", "phase_2", "preclinical"])
    parser.add_argument("--quick", action="store_true", help="reduced draws for a fast smoke run")
    args = parser.parse_args()

    if args.quick:
        args.permutations, args.bootstrap, args.weight_trials, args.proximity_permutations = 50, 25, 5, 50

    panel = json.loads((ROOT / "data/drugs/ms_panel_v3.json").read_text())
    signature = load_signature(ROOT / "data/ms_expression_v3.csv")
    config = ScoringConfig(minimum_evidence=args.min_evidence)
    drugs = eligible_drugs(panel["drugs"], config)
    by_name = {d["name"]: d for d in panel["drugs"]}
    print(f"Panel: {len(panel['drugs'])} candidates -> {len(drugs)} eligible at '{args.min_evidence}'")
    print(f"Signature: {len(signature.genes)} genes")

    print("Computing interactome exposure ...")
    exposures = build_exposures(drugs, signature.genes, args.seed, args.proximity_permutations)

    print(f"Scoring {len(drugs) * (len(drugs) - 1) // 2} pairs ...")
    rows = rank_pairs(drugs, signature, config, exposures)

    print(f"Permutation null ({args.permutations} draws/pair) ...")
    rows = permutation_pvalues(rows, by_name, signature, config, args.permutations, args.seed)
    rows = benjamini_hochberg(rows)
    rows = pareto_front(rows)
    rows.sort(key=lambda r: (r["excluded_from_primary_ranking"], -r["priority_score"], r["drug_a"], r["drug_b"]))

    print("Bootstrap stability ...")
    stability = bootstrap_stability(drugs, signature, config, exposures, args.bootstrap, args.top, args.seed)
    print("Weight sensitivity ...")
    sensitivity = weight_sensitivity(drugs, signature, config, exposures, args.weight_trials, 0.5, args.top, args.seed)
    print("Mechanism-stratum enrichment ...")
    strata = axis_stratum_enrichment(rows, by_name, drugs, signature, config, exposures, args.bootstrap, args.seed)

    primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
    diverse = diversity_filter(rows, args.top, args.max_per_drug)
    controls = control_report(rows, panel["metadata"].get("controls", {}), args.top)

    RESULTS.mkdir(parents=True, exist_ok=True)
    for index, row in enumerate(rows, start=1):
        row["rank"] = index
    with (RESULTS / "ms_pairs_full.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    significant = [r for r in primary if r["q_value"] < 0.05]
    report = {
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "git_revision": _git_revision(),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "seed": args.seed,
            "command": " ".join(sys.argv),
        },
        "inputs": {
            "panel": "data/drugs/ms_panel_v3.json",
            "panel_version": panel["metadata"]["panel_version"],
            "signature": "data/ms_expression_v3.csv",
            "network": "data/networks/string_ms_network.tsv",
            "network_meta": json.loads((ROOT / "data/networks/string_ms_network.meta.json").read_text()),
        },
        "config": config.as_dict(),
        "summary": {
            "n_candidates_total": len(panel["drugs"]),
            "n_candidates_eligible": len(drugs),
            "n_pairs_scored": len(rows),
            "n_pairs_primary": len(primary),
            "n_excluded_redundant": len(rows) - len(primary),
            "n_pareto_optimal": sum(r["pareto_optimal"] for r in rows),
            "n_complementary_exposure": sum(bool(r["complementary_exposure"]) for r in rows),
            "n_q_below_0.05": len(significant),
            "risk_domains_scored": list(RISK_DOMAINS),
        },
        "controls": controls,
        "mechanism_strata": strata,
        "stability": stability,
        "weight_sensitivity": sensitivity,
        "top_pairs_unfiltered": primary[: args.top],
        "top_pairs_diverse": diverse,
        "caveats": [
            "Discovery-stage hypothesis generation. Not a clinical combination recommendation.",
            "target_effects are curated directional hypotheses, not affinities, exposures, or doses.",
            "The expression signature is illustrative and must be replaced with a versioned, "
            "phenotype-stratified human dataset before publication.",
            "Drug targets were curated from the signature gene set, so disease-module proximity "
            "is partly circular; the separation term carries the independent topological signal.",
            "Concurrent immunomodulation carries additive infection and malignancy risk that no "
            "in-silico score can quantify.",
        ],
    }
    (RESULTS / "ms_extended_screen.json").write_text(json.dumps(report, indent=2, default=str) + "\n")

    print(f"\n{'=' * 78}\nTop {min(args.top, len(diverse))} pairs (diversity-capped at {args.max_per_drug}/agent)\n{'=' * 78}")
    print(f"{'#':>3} {'combination':<50} {'score':>7} {'q':>8} {'sep':>6}")
    for i, row in enumerate(diverse, start=1):
        combo = f"{row['drug_a']} + {row['drug_b']}"
        sep = row["network_separation"]
        sep_s = f"{sep:6.2f}" if not np.isnan(sep) else "   nan"
        print(f"{i:>3} {combo:<50} {row['priority_score']:7.3f} {row['q_value']:8.4f} {sep_s}")

    jaccard = stability[f"top{args.top}_jaccard_mean"]
    print(f"\n{'=' * 78}\nMechanism strata, ranked by median pair score (the stable unit of inference)\n{'=' * 78}")
    print(f"{'axis pair':<44} {'n':>5} {'median':>7} {'95% CI':>18} {'q':>8}")
    for stratum in strata["strata"]:
        low, high = stratum["median_bootstrap_95ci"]
        print(f"{stratum['axis_pair']:<44} {stratum['n_pairs']:>5} {stratum['median_score']:7.3f} "
              f"{f'[{low:.3f}, {high:.3f}]':>18} {stratum['q_greater']:8.2g}")

    print(f"\nPairs with q < 0.05: {len(significant)} / {len(primary)}")
    print(f"Weight-sensitivity Spearman (mean/min): {sensitivity['spearman_mean']} / {sensitivity['spearman_min']}")
    print(f"Top-{args.top} bootstrap Jaccard: {jaccard}")
    if jaccard < 0.5:
        print(
            "  WARNING: individual pair ranks are NOT stable under the curated target-effect\n"
            "  uncertainty. Report the mechanism strata above as the primary result and treat\n"
            "  named pairs as illustrative of their stratum, not as ranked recommendations."
        )
    print(f"\nWrote {RESULTS / 'ms_extended_screen.json'}")
    print(f"Wrote {RESULTS / 'ms_pairs_full.csv'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()

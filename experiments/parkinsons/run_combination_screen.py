"""Parkinson's disease: monotherapy, pair, and triple combination screen.

Discovery-stage hypothesis generation only. Nothing produced here is a clinical
combination recommendation; see docs/parkinsons_screen.md, and
docs/disease_campaign_protocol.md for the procedure this run instantiates.

The screen runs the same scoring contract at every order, so the
monotherapy-versus-combination comparison is a comparison of biology rather
than of two implementations. It answers three questions:

1. Which single agents move the Parkinson's signature furthest?
2. Which pairs beat their own best member, and by how much?
3. Does a third agent earn its place in an elderly patient's daily regimen?

Usage:
    python -m experiments.parkinsons.run_combination_screen --top 25 --seed 7
    python -m experiments.parkinsons.run_combination_screen --quick
    python -m experiments.parkinsons.run_combination_screen --max-order 2
"""

from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from core.biology.combination_scoring import (
    CombinationConfig, combination_key, eligible_drugs, rank_combinations,
    attach_subset_gain,
)
from core.biology.combination_statistics import (
    axis_stratum_enrichment, benjamini_hochberg, bootstrap_stability,
    combination_pareto_front, monotherapy_comparison, permutation_pvalues,
    weight_sensitivity,
)
from core.biology.network_proximity import ComplementaryExposure, ProximityScorer, classify_exposure
from core.models.disease import load_disease
from core.provenance import run_provenance

DISEASE = "parkinsons"
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "experiments/parkinsons/results"

CSV_COLUMNS = [
    "rank", "order", "combination", "priority_score", "p_empirical", "q_value", "permutation_z",
    "signed_reversal", "reversal_efficiency", "gene_coverage", "pathway_coverage", "axis_coverage",
    "counter_therapeutic", "best_single_reversal", "reversal_gain_over_best_single",
    "additivity_ratio", "best_subset_score", "score_gain_over_best_subset",
    "target_complementarity", "compartment_complementarity", "cns_reach", "regimen_burden",
    "network_separation", "proximity_z_min", "proximity_z_max", "complementary_exposure",
    "evidence_score", "uncertainty", "safety_union", "safety_overlap", "worst_risk_domain",
    "same_mechanism", "shared_safety_class", "redundant_targets",
    "excluded_from_primary_ranking", "pareto_optimal",
]


# ---------------------------------------------------------------------------
# Interactome exposure
# ---------------------------------------------------------------------------

def build_exposures(drugs, signature_genes, disease, seed, permutations):
    """Classify complementary exposure for every pair, then assemble triples.

    Pairwise separation is the only network quantity that has to be measured.
    A triple's binding separation is the smallest of its three pairwise values
    and its member proximities are already known, so the higher orders are
    assembled from the pairwise cache rather than recomputed -- which is what
    keeps the triple enumeration cheap enough to run exhaustively.
    """
    scorer = ProximityScorer(disease.network(), seed=seed, permutations=permutations)
    pair_exposure: dict[tuple[str, ...], ComplementaryExposure] = {}
    for a, b in combinations(drugs, 2):
        key = combination_key((a["name"], b["name"]))
        pair_exposure[key] = classify_exposure(
            scorer, list(a["target_effects"]), list(b["target_effects"]), signature_genes
        )
    return scorer, pair_exposure


def triple_exposures(pair_exposure, triples):
    """Derive triple exposure from the cached pairwise classifications."""
    out = {}
    for members in triples:
        key = combination_key(members)
        faces = [pair_exposure.get(combination_key(f)) for f in combinations(key, 2)]
        if any(f is None for f in faces):
            continue
        seps = [f.separation for f in faces]
        zs = [z for f in faces for z in (f.z_a, f.z_b)]
        finite_seps = [s for s in seps if not np.isnan(s)]
        finite_zs = [z for z in zs if not np.isnan(z)]
        complete = len(finite_seps) == len(seps) and len(finite_zs) == len(zs)
        out[key] = ComplementaryExposure(
            separation=round(min(finite_seps), 4) if finite_seps else float("nan"),
            z_a=round(min(finite_zs), 4) if finite_zs else float("nan"),
            z_b=round(max(finite_zs), 4) if finite_zs else float("nan"),
            # Every pair inside the triple must itself be complementary. The
            # conservative reading: a triple containing one overlapping pair is
            # not exploiting three distinct neighbourhoods.
            is_complementary=bool(complete and all(f.is_complementary for f in faces)),
        )
    return out


# ---------------------------------------------------------------------------
# Controls and presentation
# ---------------------------------------------------------------------------

def control_report(rows_by_order, controls, top_k):
    """Check the screen against controls declared before ranking.

    Parkinson's has an unusually well-documented record of neuroprotection
    failures, which is why the panel carries five negative-efficacy controls
    where the MS panel carries four. They test a known limitation rather than a
    pass/fail criterion: the screen scores mechanism and transcriptional
    direction, not trial outcome, so a mechanistically coherent agent that
    failed in phase 3 can still rank well. Such appearances are reported, never
    filtered.
    """
    primary = [r for order in sorted(rows_by_order) if order > 1
               for r in rows_by_order[order] if not r["excluded_from_primary_ranking"]]
    primary.sort(key=lambda r: -r["priority_score"])
    position = {}
    for i, row in enumerate(primary, start=1):
        for name in row["members"]:
            position.setdefault(name, i)
    top_names = {n for r in primary[:top_k] for n in r["members"]}

    redundant = controls.get("positive_redundancy", [])
    redundant_excluded = None
    if len(redundant) >= 2:
        key = combination_key(redundant[:2])
        match = next((r for r in rows_by_order.get(2, [])
                      if combination_key(r["members"]) == key), None)
        redundant_excluded = bool(match["excluded_from_primary_ranking"]) if match else None

    return {
        "redundant_pair": {
            "pair": redundant,
            "excluded_as_expected": redundant_excluded,
            "interpretation": "Two non-ergot D2/D3 agonists are pharmacodynamically duplicate "
                              "and must be excluded from the primary ranking.",
        },
        "safety_controls": {
            name: {"best_primary_rank": position.get(name), "in_top_k": name in top_names}
            for name in controls.get("safety_penalty", [])
        },
        "negative_efficacy_controls": {
            name: {"best_primary_rank": position.get(name), "in_top_k": name in top_names}
            for name in controls.get("negative_efficacy", [])
        },
        "note": (
            "Negative-efficacy controls are agents whose Parkinson's trials failed "
            "(Creatine, Coenzyme Q10, Isradipine, Inosine, Cinpanemab). The screen scores "
            "mechanism and transcriptional direction rather than trial outcome, so some may "
            "still rank highly; that is a stated limitation of the method and any such "
            "appearance is reported rather than filtered away."
        ),
    }


def diversity_filter(rows, top_k, max_per_drug):
    """Greedy top-K selection capping how often any one agent may appear.

    Without it, levodopa and the LRRK2 inhibitor occupy most of the leaderboard
    and hide the mechanistic diversity the screen exists to surface.
    """
    counts, chosen = {}, []
    for row in rows:
        if row["excluded_from_primary_ranking"]:
            continue
        if any(counts.get(n, 0) >= max_per_drug for n in row["members"]):
            continue
        for n in row["members"]:
            counts[n] = counts.get(n, 0) + 1
        chosen.append(row)
        if len(chosen) >= top_k:
            break
    return chosen


def _fmt(value, width=7, places=3):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return " " * (width - 3) + "nan"
    return f"{value:{width}.{places}f}"


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--max-order", type=int, default=3, choices=[1, 2, 3],
                        help="highest combination order to enumerate")
    parser.add_argument("--permutations", type=int, default=400,
                        help="draws for the significance null at order 1 and 2")
    parser.add_argument("--triple-permutations", type=int, default=200,
                        help="draws at order 3, where there are an order of magnitude more rows")
    parser.add_argument("--proximity-permutations", type=int, default=200)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--weight-trials", type=int, default=40)
    parser.add_argument("--max-per-drug", type=int, default=3)
    parser.add_argument("--min-evidence", default="phase_2",
                        choices=["approved", "phase_3", "phase_2", "preclinical"])
    parser.add_argument("--quick", action="store_true", help="reduced draws for a fast smoke run")
    args = parser.parse_args()

    if args.quick:
        args.permutations, args.triple_permutations = 50, 25
        args.bootstrap, args.weight_trials, args.proximity_permutations = 25, 5, 50

    disease = load_disease(DISEASE)
    signature = disease.signature()
    panel = disease.panel()
    panel_metadata = disease.panel_metadata()
    config = CombinationConfig(minimum_evidence=args.min_evidence)
    drugs = eligible_drugs(panel, config)
    by_name = {d["name"]: d for d in panel}

    print(f"Disease: {disease.name}")
    print(f"Panel: {len(panel)} candidates -> {len(drugs)} eligible at '{args.min_evidence}'")
    print(f"Signature: {len(signature.genes)} genes across {len(disease.pathways)} pathways")
    print(f"Axes: {', '.join(disease.therapeutic_axes)}")
    print(f"Risk domains: {', '.join(disease.risk_domains)}\n")

    print("Computing interactome exposure ...")
    _, pair_exposure = build_exposures(drugs, signature.genes, disease, args.seed,
                                       args.proximity_permutations)

    rows_by_order: dict[int, list[dict]] = {}
    stability: dict[int, dict] = {}
    sensitivity: dict[int, dict] = {}
    strata: dict[int, dict] = {}
    prefilter: set[tuple[str, ...]] | None = None
    exposures: dict[tuple[str, ...], ComplementaryExposure] = {}

    for order in range(1, args.max_order + 1):
        if order == 2:
            exposures = pair_exposure
        elif order == 3:
            pool = [d["name"] for d in drugs]
            exposures = triple_exposures(pair_exposure, combinations(sorted(pool), 3))

        # The prefilter that built *this* order, held separately: `prefilter` is
        # reassigned below to seed the next order, and the robustness analyses
        # must re-enumerate exactly the space that was scored here.
        order_prefilter = prefilter
        print(f"\n--- order {order} ---")
        rows = rank_combinations(drugs, signature, disease, config, order, exposures, order_prefilter)
        print(f"Scored {len(rows)} combinations of order {order}")
        if order > 1:
            rows = attach_subset_gain(rows, rows_by_order[order - 1])

        draws = args.triple_permutations if order >= 3 else args.permutations
        print(f"Permutation null ({draws} draws/combination) ...")
        rows = permutation_pvalues(rows, by_name, signature, disease, config, draws, args.seed)
        rows = benjamini_hochberg(rows)
        rows = combination_pareto_front(rows)
        rows.sort(key=lambda r: (r["excluded_from_primary_ranking"], -r["priority_score"],
                                 r["combination"]))
        rows_by_order[order] = rows

        primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
        # Only combinations that survived redundancy exclusion may seed the next
        # order, so a triple can never be built around a duplicate pair.
        prefilter = {combination_key(r["members"]) for r in primary}

        if order > 1:
            print("Bootstrap stability ...")
            stability[order] = bootstrap_stability(
                drugs, signature, disease, config, order, exposures,
                prefilter=order_prefilter,
                draws=args.bootstrap, top_k=args.top, seed=args.seed)
            print("Weight sensitivity ...")
            sensitivity[order] = weight_sensitivity(
                drugs, signature, disease, config, order, exposures,
                prefilter=order_prefilter,
                trials=args.weight_trials, top_k=args.top, seed=args.seed)
            print("Mechanism-stratum enrichment ...")
            strata[order] = axis_stratum_enrichment(
                rows, by_name, drugs, signature, disease, config, order, exposures,
                prefilter=order_prefilter,
                bootstrap_draws=args.bootstrap, seed=args.seed)

    comparison = monotherapy_comparison(rows_by_order, top_k=10)
    controls = control_report(rows_by_order, panel_metadata.get("controls", {}), args.top)

    # ----------------------------------------------------------------- output
    RESULTS.mkdir(parents=True, exist_ok=True)
    all_rows = []
    for order in sorted(rows_by_order):
        for index, row in enumerate(rows_by_order[order], start=1):
            row["rank"] = index
        all_rows.extend(rows_by_order[order])
    csv_path = RESULTS / "pd_combinations_full.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    report = {
        "provenance": run_provenance(
            [disease.signature_path, disease.panel_path, disease.network_path],
            command=["python", "-m", "experiments.parkinsons.run_combination_screen",
                     *[a for a in __import__("sys").argv[1:]]],
            extra={"seed": args.seed, "disease": DISEASE},
        ),
        "disease": disease.as_dict(),
        "inputs": {
            "signature": str(disease.signature_path.relative_to(ROOT)),
            "panel": str(disease.panel_path.relative_to(ROOT)),
            "panel_version": panel_metadata.get("panel_version"),
            "network": str(disease.network_path.relative_to(ROOT)),
        },
        "config": config.as_dict(),
        "summary": {
            "n_candidates_total": len(panel),
            "n_candidates_eligible": len(drugs),
            "orders": {
                str(order): {
                    "n_scored": len(rows),
                    "n_primary": sum(not r["excluded_from_primary_ranking"] for r in rows),
                    "n_excluded_redundant": sum(r["excluded_from_primary_ranking"] for r in rows),
                    "n_pareto_optimal": sum(r["pareto_optimal"] for r in rows),
                    "n_complementary_exposure": sum(bool(r["complementary_exposure"]) for r in rows),
                    "n_q_below_0.05": sum(r["q_value"] < 0.05 for r in rows
                                          if not r["excluded_from_primary_ranking"]),
                }
                for order, rows in rows_by_order.items()
            },
        },
        "monotherapy_vs_combination": comparison,
        "controls": controls,
        "mechanism_strata": {str(k): v for k, v in strata.items()},
        "stability": {str(k): v for k, v in stability.items()},
        "weight_sensitivity": {str(k): v for k, v in sensitivity.items()},
        "top_by_order": {
            str(order): [r for r in rows if not r["excluded_from_primary_ranking"]][: args.top]
            for order, rows in rows_by_order.items()
        },
        "top_diverse_by_order": {
            str(order): diversity_filter(rows, args.top, args.max_per_drug)
            for order, rows in rows_by_order.items() if order > 1
        },
        "caveats": [
            "Discovery-stage hypothesis generation. Not a clinical combination recommendation.",
            "target_effects are curated directional hypotheses, not affinities, exposures, or doses.",
            "The expression signature is curated from literature rather than derived from a "
            "cohort, and is weaker evidence than the MS entry's. Replace it with a versioned, "
            "stage-stratified human dataset before publication use.",
            "No pharmacokinetic or drug-drug-interaction model. Co-administration cost is "
            "represented only by route burden and half-life spread in regimen_burden. This is "
            "the largest gap for a disease treated with chronic oral polypharmacy in an "
            "elderly population, and is the top-priority next model.",
            "priority_score is comparable only within a fixed order. Across orders, compare "
            "the efficacy block; see monotherapy_vs_combination.",
            "Panel targets were curated from the signature gene set, so disease-module "
            "proximity is partly circular; the separation term carries the independent signal.",
            "Adding a third chronic medication to an elderly patient carries adherence, "
            "falls, and anticholinergic-burden costs that no in-silico score can quantify.",
        ],
    }
    json_path = RESULTS / "pd_combination_screen.json"
    json_path.write_text(json.dumps(report, indent=2, default=str) + "\n")

    # ---------------------------------------------------------------- printout
    line = "=" * 88
    print(f"\n{line}\nMonotherapy vs combination -- efficacy block only (signed reversal)\n{line}")
    print(f"best monotherapy signed_reversal = {comparison['best_monotherapy_signed_reversal']}")
    print(f"{'order':>5} {'n':>6} {'median rev':>11} {'max rev':>9} {'>best mono':>11} "
          f"{'max gain':>9} {'sub-add':>9} {'p':>9}")
    for order, block in comparison["by_order"].items():
        dist = block["efficacy_distribution"]["signed_reversal"]
        p = block["p_greater_than_monotherapy"]
        print(f"{order:>5} {block['n']:>6} {dist['median']:>11.4f} {dist['max']:>9.4f} "
              f"{block['fraction_beating_best_monotherapy']:>10.1%} "
              f"{block['max_reversal_gain_over_best_monotherapy']:>9.4f} "
              f"{block['fraction_subadditive']:>8.1%} "
              f"{('%.2g' % p) if p is not None else '-':>9}")

    for order in sorted(rows_by_order):
        rows = rows_by_order[order]
        primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
        shown = primary[: args.top] if order == 1 else diversity_filter(rows, args.top, args.max_per_drug)
        label = ("monotherapy" if order == 1 else
                 f"order-{order} combinations (diversity-capped at {args.max_per_drug}/agent)")
        print(f"\n{line}\nTop {len(shown)} {label}\n{line}")
        header = f"{'#':>3} {'combination':<52} {'score':>7} {'rev':>7} {'q':>8} {'sep':>7}"
        if order > 1:
            header += f" {'gain/best':>10}"
        if order > 2:
            header += f" {'gain/pair':>10}"
        print(header)
        for i, row in enumerate(shown, start=1):
            out = (f"{i:>3} {row['combination'][:52]:<52} {row['priority_score']:7.3f} "
                   f"{row['signed_reversal']:7.4f} {row['q_value']:8.4f} "
                   f"{_fmt(row['network_separation'])}")
            if order > 1:
                out += f" {row['reversal_gain_over_best_single']:10.4f}"
            if order > 2:
                out += f" {_fmt(row.get('score_gain_over_best_subset'), 10, 4)}"
            print(out)

    for order, block in strata.items():
        print(f"\n{line}\nOrder-{order} mechanism strata, ranked by median score "
              f"(the stable unit of inference)\n{line}")
        print(f"{'axis stratum':<58} {'n':>6} {'median':>7} {'95% CI':>18} {'q':>8}")
        for stratum in block["strata"]:
            low, high = stratum["median_bootstrap_95ci"]
            print(f"{stratum['axis_stratum'][:58]:<58} {stratum['n_combinations']:>6} "
                  f"{stratum['median_score']:7.3f} {f'[{low:.3f}, {high:.3f}]':>18} "
                  f"{stratum['q_greater']:8.2g}")

    print(f"\n{line}\nRobustness\n{line}")
    for order in sorted(stability):
        jac = stability[order][f"top{args.top}_jaccard_mean"]
        spear = sensitivity[order]["spearman_mean"]
        print(f"order {order}: top-{args.top} bootstrap Jaccard = {jac}; "
              f"weight-sensitivity Spearman (mean/min) = {spear} / "
              f"{sensitivity[order]['spearman_min']}")
        if jac < 0.5:
            print("  WARNING: individual ranks at this order are NOT stable under the curated\n"
                  "  target-effect uncertainty. Report the mechanism strata above as the primary\n"
                  "  result and treat named combinations as illustrative of their stratum.")

    print(f"\nWrote {json_path}")
    print(f"Wrote {csv_path} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()

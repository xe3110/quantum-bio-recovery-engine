"""Significance, stability, and stratum inference for combinations of any order.

The pairwise screen's statistics live in ``core.biology.screen_statistics`` and
are written against ``(drug_a, drug_b)`` rows. This module does the same work
for the k-ary rows produced by ``core.biology.combination_scoring``, and reuses
that module's two genuinely order-independent pieces -- Benjamini-Hochberg and
the Pareto front -- rather than restating them.

A ranked list is not a result on its own. What a reviewer will ask for is:

* **Is the score better than chance?** A per-combination permutation null that
  preserves each agent's target-set size and its multiset of effect
  magnitudes, so the null isolates *which* genes an agent hits from *how hard*
  and *how many*. Converted to an empirical p-value, then to a BH q-value.
* **Is the ranking stable?** Bootstrap resampling of the curated target effects
  under each record's own stated ``target_uncertainty``, reported as score and
  rank intervals plus a top-K Jaccard index.
* **Does anything survive a different weighting?** A weight-free Pareto front
  and a Spearman sweep with every weight jittered.
* **What is stable, if the individual ranks are not?** Aggregation by
  therapeutic-axis stratum, which is the level at which this class of screen
  has historically had something defensible to say.
"""

from __future__ import annotations

from itertools import product
from typing import Any, Sequence

import numpy as np

from core.biology.combination_scoring import (
    CombinationConfig, PARETO_OBJECTIVES, combination_key, combination_metrics,
    rank_combinations,
)
from core.biology.screen_statistics import _permute_drug, _stable_seed, benjamini_hochberg, pareto_front  # noqa: F401
from core.biology.signature import Signature
from core.models.disease import DiseaseContext

__all__ = [
    "benjamini_hochberg", "pareto_front", "combination_pareto_front",
    "permutation_pvalues", "bootstrap_stability", "axis_stratum_enrichment",
    "weight_sensitivity", "monotherapy_comparison",
]


def combination_pareto_front(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Pareto front over the combination objectives (adds regimen burden)."""
    return pareto_front(rows, objectives=PARETO_OBJECTIVES)


def _perturb(drugs: Sequence[dict], rng: np.random.Generator) -> list[dict]:
    """Resample every curated target effect under its record's own uncertainty."""
    perturbed = []
    for drug in drugs:
        sd = float(drug["target_uncertainty"])
        clone = dict(drug)
        clone["target_effects"] = {
            g: float(np.clip(v + rng.normal(0.0, sd), -1.0, 1.0))
            for g, v in drug["target_effects"].items()
        }
        perturbed.append(clone)
    return perturbed


# ---------------------------------------------------------------------------
# Permutation significance
# ---------------------------------------------------------------------------

def permutation_pvalues(
    rows: list[dict[str, Any]],
    by_name: dict[str, dict],
    signature: Signature,
    disease: DiseaseContext,
    config: CombinationConfig,
    draws: int = 400,
    seed: int = 7,
) -> list[dict[str, Any]]:
    """Attach an empirical p-value to each row (add-one corrected, never 0).

    Topology is deliberately excluded from the null: ``exposure=None`` means
    the observed and null scores are compared on the terms the permutation
    actually randomises. Including a fixed separation term on both sides would
    inflate both equally and tell a reader nothing.
    """
    universe = signature.genes
    for row in rows:
        rng = np.random.default_rng(_stable_seed(seed, *row["members"]))
        members = [by_name[n] for n in row["members"]]
        null = np.empty(draws)
        for i in range(draws):
            shuffled = [_permute_drug(d, universe, rng) for d in members]
            null[i] = combination_metrics(
                shuffled, signature, disease, config, exposure=None
            )["priority_score"]
        observed = row["priority_score"]
        sd = float(null.std())
        row["null_mean"] = round(float(null.mean()), 4)
        row["null_sd"] = round(sd, 4)
        row["permutation_z"] = round(float((observed - null.mean()) / sd), 4) if sd > 0 else 0.0
        row["p_empirical"] = round(float((np.sum(null >= observed) + 1) / (draws + 1)), 5)
    return rows


# ---------------------------------------------------------------------------
# Bootstrap stability
# ---------------------------------------------------------------------------

def bootstrap_stability(
    drugs: Sequence[dict],
    signature: Signature,
    disease: DiseaseContext,
    config: CombinationConfig,
    order: int,
    exposures: dict,
    prefilter: set[tuple[str, ...]] | None = None,
    draws: int = 200,
    top_k: int = 25,
    seed: int = 7,
) -> dict[str, Any]:
    """Resample curated target effects and measure how much the ranking moves."""
    rng = np.random.default_rng(seed)
    base = [r for r in rank_combinations(drugs, signature, disease, config, order, exposures, prefilter)
            if not r["excluded_from_primary_ranking"]]
    ordered_top = [combination_key(r["members"]) for r in base][:top_k]
    top_set = set(ordered_top)

    scores: dict[tuple, list[float]] = {}
    ranks: dict[tuple, list[int]] = {}
    jaccards = []
    for _ in range(draws):
        trial = [r for r in rank_combinations(
            _perturb(drugs, rng), signature, disease, config, order, exposures, prefilter)
            if not r["excluded_from_primary_ranking"]]
        for position, row in enumerate(trial, start=1):
            key = combination_key(row["members"])
            scores.setdefault(key, []).append(row["priority_score"])
            ranks.setdefault(key, []).append(position)
        trial_top = {combination_key(r["members"]) for r in trial[:top_k]}
        union = top_set | trial_top
        jaccards.append(len(top_set & trial_top) / len(union) if union else 0.0)

    intervals = {}
    for key in ordered_top:
        s = np.array(scores.get(key, [np.nan]))
        r = np.array(ranks.get(key, [np.nan]))
        intervals[" + ".join(key)] = {
            "score_95ci": [round(float(np.nanquantile(s, 0.025)), 4),
                           round(float(np.nanquantile(s, 0.975)), 4)],
            "rank_95ci": [int(np.nanquantile(r, 0.025)), int(np.nanquantile(r, 0.975))],
            "rank_median": int(np.nanmedian(r)),
        }
    return {
        "draws": draws,
        "order": order,
        "top_k": top_k,
        f"top{top_k}_jaccard_mean": round(float(np.mean(jaccards)), 4) if jaccards else 0.0,
        f"top{top_k}_jaccard_min": round(float(np.min(jaccards)), 4) if jaccards else 0.0,
        "intervals": intervals,
    }


# ---------------------------------------------------------------------------
# Weight sensitivity
# ---------------------------------------------------------------------------

def weight_sensitivity(
    drugs: Sequence[dict],
    signature: Signature,
    disease: DiseaseContext,
    base_config: CombinationConfig,
    order: int,
    exposures: dict,
    prefilter: set[tuple[str, ...]] | None = None,
    trials: int = 40,
    jitter: float = 0.5,
    top_k: int = 25,
    seed: int = 7,
) -> dict[str, Any]:
    """Jitter every weight multiplicatively and measure ranking agreement.

    The weights were chosen by hand. This reports how much of the conclusion
    depends on that choice, as a Spearman correlation of the full ranking and a
    top-K overlap against the declared weights.
    """
    from scipy.stats import spearmanr

    rng = np.random.default_rng(seed)
    base = rank_combinations(drugs, signature, disease, base_config, order, exposures, prefilter)
    base_keys = [combination_key(r["members"]) for r in base]
    base_scores = {k: r["priority_score"] for k, r in zip(base_keys, base)}
    base_top = set(base_keys[:top_k])

    tunable = [f for f in vars(base_config) if f.endswith(("_weight", "_penalty"))]
    correlations, overlaps = [], []
    for _ in range(trials):
        kwargs = {f: max(0.0, getattr(base_config, f) * float(rng.uniform(1 - jitter, 1 + jitter)))
                  for f in tunable}
        cfg = CombinationConfig(
            minimum_evidence=base_config.minimum_evidence,
            max_target_overlap=base_config.max_target_overlap,
            broad_coverage_reference=base_config.broad_coverage_reference,
            **kwargs,
        )
        trial = rank_combinations(drugs, signature, disease, cfg, order, exposures, prefilter)
        trial_scores = {combination_key(r["members"]): r["priority_score"] for r in trial}
        correlations.append(spearmanr(
            [base_scores[k] for k in base_keys], [trial_scores[k] for k in base_keys]
        ).statistic)
        trial_top = {combination_key(r["members"]) for r in trial[:top_k]}
        overlaps.append(len(base_top & trial_top) / max(1, top_k))
    return {
        "trials": trials,
        "weight_jitter": jitter,
        "spearman_mean": round(float(np.mean(correlations)), 4),
        "spearman_min": round(float(np.min(correlations)), 4),
        f"top{top_k}_overlap_mean": round(float(np.mean(overlaps)), 4),
        f"top{top_k}_overlap_min": round(float(np.min(overlaps)), 4),
    }


# ---------------------------------------------------------------------------
# Mechanism-stratum inference
# ---------------------------------------------------------------------------

def _strata_for(row: dict, by_name: dict[str, dict]) -> set[tuple[str, ...]]:
    """Every axis tuple a combination realises, as an unordered label.

    A combination of agents each carrying several therapeutic axes belongs to
    several strata at once; all of them are credited, which is why stratum
    counts sum to more than the number of combinations.
    """
    axis_lists = [by_name[n].get("therapeutic_axes", []) or ["unspecified"] for n in row["members"]]
    return {tuple(sorted(combo)) for combo in product(*axis_lists)}


def axis_stratum_enrichment(
    rows: list[dict[str, Any]],
    by_name: dict[str, dict],
    drugs: Sequence[dict],
    signature: Signature,
    disease: DiseaseContext,
    config: CombinationConfig,
    order: int,
    exposures: dict,
    prefilter: set[tuple[str, ...]] | None = None,
    bootstrap_draws: int = 200,
    min_members: int = 5,
    seed: int = 7,
) -> dict[str, Any]:
    """Aggregate scores by therapeutic-axis stratum and test each stratum.

    Individual combination ranks are unstable under curated target-effect
    uncertainty, so the defensible unit of inference is the *mechanism
    stratum*. For every axis tuple this reports the score distribution, a
    Mann-Whitney test against all other strata, and a bootstrap interval on the
    stratum median.
    """
    from scipy.stats import mannwhitneyu

    primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
    buckets: dict[tuple[str, ...], list[float]] = {}
    for row in primary:
        for stratum in _strata_for(row, by_name):
            buckets.setdefault(stratum, []).append(row["priority_score"])

    rng = np.random.default_rng(seed)
    medians: dict[tuple[str, ...], list[float]] = {k: [] for k in buckets}
    for _ in range(bootstrap_draws):
        trial = [r for r in rank_combinations(
            _perturb(drugs, rng), signature, disease, config, order, exposures, prefilter)
            if not r["excluded_from_primary_ranking"]]
        trial_buckets: dict[tuple[str, ...], list[float]] = {}
        for row in trial:
            for stratum in _strata_for(row, by_name):
                trial_buckets.setdefault(stratum, []).append(row["priority_score"])
        for stratum, values in trial_buckets.items():
            if stratum in medians:
                medians[stratum].append(float(np.median(values)))

    results = []
    for stratum, values in buckets.items():
        others = [v for k, vals in buckets.items() if k != stratum for v in vals]
        if len(values) >= min_members and others:
            stat = mannwhitneyu(values, others, alternative="greater")
            p = float(stat.pvalue)
            # Rank-biserial effect: probability a stratum member outscores a
            # non-member, interpretable without reference to the sample size.
            effect = float(stat.statistic) / (len(values) * len(others))
        else:
            p, effect = float("nan"), float("nan")
        boot = np.array(medians.get(stratum, []))
        results.append({
            "axis_stratum": " + ".join(stratum),
            "n_combinations": len(values),
            "median_score": round(float(np.median(values)), 4),
            "mean_score": round(float(np.mean(values)), 4),
            "median_bootstrap_95ci": (
                [round(float(np.quantile(boot, 0.025)), 4), round(float(np.quantile(boot, 0.975)), 4)]
                if boot.size > 1 else [float("nan"), float("nan")]
            ),
            "p_greater": round(p, 6) if not np.isnan(p) else None,
            "prob_outscores_other_strata": round(effect, 4) if not np.isnan(effect) else None,
        })

    results = [r for r in results if r["p_greater"] is not None]
    benjamini_hochberg(results, key="p_greater", out="q_greater")
    results.sort(key=lambda r: -r["median_score"])
    return {
        "bootstrap_draws": bootstrap_draws,
        "order": order,
        "interpretation": (
            "Stratum medians carry bootstrap intervals that are narrow relative to the "
            "spread between strata, whereas individual combination ranks are not stable. "
            "Treat the mechanism stratum, not the named combination, as the unit of "
            "inference."
        ),
        "strata": results,
    }


# ---------------------------------------------------------------------------
# Monotherapy versus combination
# ---------------------------------------------------------------------------

# Defined identically at every order, so these are the only fields on which a
# single agent and a combination may be compared. priority_score is not among
# them, by construction: a single agent scores zero on the complementarity and
# separation terms and would be penalised for being one drug rather than for
# being a worse one.
EFFICACY_FIELDS = (
    "signed_reversal", "reversal_efficiency", "gene_coverage",
    "pathway_coverage", "axis_coverage", "counter_therapeutic",
)


def monotherapy_comparison(
    by_order: dict[int, list[dict[str, Any]]],
    top_k: int = 10,
) -> dict[str, Any]:
    """Compare orders on the efficacy block, and only on the efficacy block.

    Reports the full distribution per order, the best combination found at each
    order, and -- the question the whole screen is built to answer -- what
    fraction of combinations at order k beat *every* monotherapy on signed
    reversal, and by how much.
    """
    from scipy.stats import mannwhitneyu

    singles = [r for r in by_order.get(1, [])]
    best_mono = max((r["signed_reversal"] for r in singles), default=0.0)
    mono_reversals = [r["signed_reversal"] for r in singles]

    per_order = {}
    for order in sorted(by_order):
        rows = [r for r in by_order[order] if not r["excluded_from_primary_ranking"]] or by_order[order]
        reversals = np.array([r["signed_reversal"] for r in rows], dtype=float)
        beats = [r for r in rows if r["signed_reversal"] > best_mono]
        if order > 1 and mono_reversals and reversals.size:
            p = float(mannwhitneyu(reversals, mono_reversals, alternative="greater").pvalue)
        else:
            p = None
        per_order[order] = {
            "n": len(rows),
            "efficacy_distribution": {
                field: {
                    "median": round(float(np.median([r[field] for r in rows])), 4),
                    "max": round(float(np.max([r[field] for r in rows])), 4),
                }
                for field in EFFICACY_FIELDS
            },
            "n_beating_best_monotherapy": len(beats),
            "fraction_beating_best_monotherapy": round(len(beats) / len(rows), 4) if rows else 0.0,
            "max_reversal_gain_over_best_monotherapy": round(
                float(max((r["signed_reversal"] - best_mono for r in rows), default=0.0)), 4),
            # The median is 1.0 wherever most members have disjoint target sets,
            # which makes it the uninformative half of this statistic. What
            # matters is the tail: the fraction of combinations whose members
            # are covering the same signal, and how badly.
            "median_additivity_ratio": round(
                float(np.median([r["additivity_ratio"] for r in rows])), 4),
            "min_additivity_ratio": round(
                float(np.min([r["additivity_ratio"] for r in rows])), 4),
            "p05_additivity_ratio": round(
                float(np.quantile([r["additivity_ratio"] for r in rows], 0.05)), 4),
            "fraction_subadditive": round(
                float(np.mean([r["additivity_ratio"] < 0.999 for r in rows])), 4),
            "p_greater_than_monotherapy": round(p, 6) if p is not None else None,
            "best": [
                {k: r[k] for k in ("combination", "priority_score", *EFFICACY_FIELDS,
                                   "additivity_ratio", "reversal_gain_over_best_single",
                                   "safety_union", "regimen_burden")}
                for r in sorted(rows, key=lambda r: -r["signed_reversal"])[:top_k]
            ],
        }
    return {
        "best_monotherapy_signed_reversal": round(float(best_mono), 4),
        "comparable_fields": list(EFFICACY_FIELDS),
        "note": (
            "Orders are compared on the efficacy block only. priority_score includes "
            "terms a single agent cannot earn (target and compartment complementarity, "
            "network separation) and is meaningful only within a fixed order."
        ),
        "by_order": per_order,
    }

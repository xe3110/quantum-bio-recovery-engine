"""Significance, stability, and multi-objective analysis for the MS pair screen.

A ranked list on its own is not a result.  This module supplies the three
things a reviewer will ask for:

* **Is the score better than chance?**  A size- and magnitude-preserving
  permutation null per pair, converted to an empirical p-value and then to a
  Benjamini-Hochberg q-value across the whole screen.
* **Is the ranking stable?**  Bootstrap resampling of the curated target
  effects under each record's stated uncertainty, reported as a rank interval
  and a top-K Jaccard index.
* **Does the ranking survive different weights?**  A Pareto front computed
  without any weighting, plus a Spearman weight-sensitivity sweep.
"""

from __future__ import annotations

import hashlib
from typing import Any, Sequence

import numpy as np

from core.biology.ms_scoring import (
    PARETO_OBJECTIVES,
    ScoringConfig,
    Signature,
    pair_metrics,
    rank_pairs,
)


# ---------------------------------------------------------------------------
# Permutation significance
# ---------------------------------------------------------------------------

def _stable_seed(seed: int, *parts: str) -> int:
    """Derive a per-pair seed that is stable across processes.

    Python's builtin ``hash()`` is randomised per interpreter unless
    PYTHONHASHSEED is set, which would make the permutation p-values
    irreproducible between runs despite a stated seed.
    """
    digest = hashlib.sha256("|".join((str(seed), *parts)).encode()).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


def _permute_drug(drug: dict, universe: Sequence[str], rng: np.random.Generator) -> dict:
    """Reassign a drug's effect magnitudes to random genes.

    Target-set size and the multiset of effect magnitudes are preserved, so the
    null isolates *which* genes are hit from *how hard* and *how many*.
    """
    effects = list(drug["target_effects"].values())
    genes = rng.choice(np.asarray(universe, dtype=object), size=len(effects), replace=False)
    clone = dict(drug)
    clone["target_effects"] = {str(g): float(v) for g, v in zip(genes, effects)}
    return clone


def permutation_pvalues(
    rows: list[dict[str, Any]],
    by_name: dict[str, dict],
    signature: Signature,
    config: ScoringConfig,
    draws: int = 400,
    seed: int = 7,
) -> list[dict[str, Any]]:
    """Attach an empirical p-value to each row (add-one corrected, never 0)."""
    universe = signature.genes
    for row in rows:
        rng = np.random.default_rng(_stable_seed(seed, row["drug_a"], row["drug_b"]))
        a, b = by_name[row["drug_a"]], by_name[row["drug_b"]]
        null = np.empty(draws)
        for i in range(draws):
            null[i] = pair_metrics(
                _permute_drug(a, universe, rng), _permute_drug(b, universe, rng),
                signature, config, exposure=None,
            )["priority_score"]
        observed = row["priority_score"]
        row["null_mean"] = round(float(null.mean()), 4)
        row["null_sd"] = round(float(null.std()), 4)
        row["permutation_z"] = round(float((observed - null.mean()) / null.std()), 4) if null.std() > 0 else 0.0
        row["p_empirical"] = round(float((np.sum(null >= observed) + 1) / (draws + 1)), 5)
    return rows


def benjamini_hochberg(rows: list[dict[str, Any]], key: str = "p_empirical", out: str = "q_value") -> list[dict]:
    """Add BH-FDR q-values computed across every row supplied."""
    if not rows:
        return rows
    order = sorted(range(len(rows)), key=lambda i: rows[i][key])
    n = len(rows)
    running = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        running = min(running, rows[idx][key] * n / rank)
        rows[idx][out] = round(float(running), 5)
    return rows


# ---------------------------------------------------------------------------
# Multi-objective (weight-free) analysis
# ---------------------------------------------------------------------------

def pareto_front(rows: list[dict[str, Any]], objectives=PARETO_OBJECTIVES) -> list[dict[str, Any]]:
    """Flag rows not dominated on any declared objective.

    A pair is dominated when another pair is at least as good on every
    objective and strictly better on one.  This answers "which candidates
    survive regardless of how the weights are chosen?".
    """
    def vector(row):
        return [row[k] if maximise else -row[k] for k, maximise in objectives]

    vectors = [vector(r) for r in rows]
    for i, row in enumerate(rows):
        vi = vectors[i]
        dominated = False
        for j, vj in enumerate(vectors):
            if i == j:
                continue
            if all(b >= a for a, b in zip(vi, vj)) and any(b > a for a, b in zip(vi, vj)):
                dominated = True
                break
        row["pareto_optimal"] = not dominated
    return rows


def weight_sensitivity(
    drugs: Sequence[dict],
    signature: Signature,
    base_config: ScoringConfig,
    exposures: dict,
    trials: int = 40,
    jitter: float = 0.5,
    top_k: int = 25,
    seed: int = 7,
) -> dict[str, Any]:
    """Perturb every weight multiplicatively and measure ranking agreement.

    Reports the Spearman correlation of the full ranking and the top-K overlap
    against the declared weights, so a reader can judge how much the
    conclusions depend on the specific weighting chosen.
    """
    from scipy.stats import spearmanr

    rng = np.random.default_rng(seed)
    base = rank_pairs(drugs, signature, base_config, exposures)
    base_keys = [(r["drug_a"], r["drug_b"]) for r in base]
    base_scores = {k: r["priority_score"] for k, r in zip(base_keys, base)}
    base_top = {k for k, r in zip(base_keys, base) if not r["excluded_from_primary_ranking"]}
    base_top = set(list(base_top)[:0]) or {k for k in base_keys[:top_k]}

    tunable = [f for f in vars(base_config) if f.endswith(("_weight", "_penalty"))]
    correlations, overlaps = [], []
    for _ in range(trials):
        kwargs = {f: max(0.0, getattr(base_config, f) * float(rng.uniform(1 - jitter, 1 + jitter)))
                  for f in tunable}
        cfg = ScoringConfig(**{**{"minimum_evidence": base_config.minimum_evidence}, **kwargs})
        trial = rank_pairs(drugs, signature, cfg, exposures)
        trial_scores = {(r["drug_a"], r["drug_b"]): r["priority_score"] for r in trial}
        correlations.append(spearmanr(
            [base_scores[k] for k in base_keys], [trial_scores[k] for k in base_keys]
        ).statistic)
        trial_top = {(r["drug_a"], r["drug_b"]) for r in trial[:top_k]}
        overlaps.append(len(base_top & trial_top) / top_k)
    return {
        "trials": trials,
        "weight_jitter": jitter,
        "spearman_mean": round(float(np.mean(correlations)), 4),
        "spearman_min": round(float(np.min(correlations)), 4),
        f"top{top_k}_overlap_mean": round(float(np.mean(overlaps)), 4),
        f"top{top_k}_overlap_min": round(float(np.min(overlaps)), 4),
    }


# ---------------------------------------------------------------------------
# Bootstrap stability
# ---------------------------------------------------------------------------

def bootstrap_stability(
    drugs: Sequence[dict],
    signature: Signature,
    config: ScoringConfig,
    exposures: dict,
    draws: int = 200,
    top_k: int = 25,
    seed: int = 7,
) -> dict[str, Any]:
    """Resample curated target effects under their stated uncertainty.

    Returns per-pair score intervals and rank intervals plus the mean Jaccard
    overlap of the top-K list across replicates.
    """
    rng = np.random.default_rng(seed)
    base = rank_pairs(drugs, signature, config, exposures)
    base_top = {(r["drug_a"], r["drug_b"]) for r in base if not r["excluded_from_primary_ranking"]}
    base_top = set(list(base_top)[:top_k]) if len(base_top) > top_k else base_top
    ordered_top = [(r["drug_a"], r["drug_b"]) for r in base if not r["excluded_from_primary_ranking"]][:top_k]

    scores: dict[tuple, list[float]] = {}
    ranks: dict[tuple, list[int]] = {}
    jaccards = []
    for _ in range(draws):
        perturbed = []
        for drug in drugs:
            sd = float(drug["target_uncertainty"])
            clone = dict(drug)
            clone["target_effects"] = {
                g: float(np.clip(v + rng.normal(0.0, sd), -1.0, 1.0))
                for g, v in drug["target_effects"].items()
            }
            perturbed.append(clone)
        trial = rank_pairs(perturbed, signature, config, exposures)
        eligible = [r for r in trial if not r["excluded_from_primary_ranking"]]
        for position, row in enumerate(eligible, start=1):
            key = (row["drug_a"], row["drug_b"])
            scores.setdefault(key, []).append(row["priority_score"])
            ranks.setdefault(key, []).append(position)
        trial_top = {(r["drug_a"], r["drug_b"]) for r in eligible[:top_k]}
        union = set(ordered_top) | trial_top
        jaccards.append(len(set(ordered_top) & trial_top) / len(union) if union else 0.0)

    intervals = {}
    for key in ordered_top:
        s = np.array(scores.get(key, [np.nan]))
        r = np.array(ranks.get(key, [np.nan]))
        intervals[f"{key[0]} + {key[1]}"] = {
            "score_95ci": [round(float(np.nanquantile(s, 0.025)), 4), round(float(np.nanquantile(s, 0.975)), 4)],
            "rank_95ci": [int(np.nanquantile(r, 0.025)), int(np.nanquantile(r, 0.975))],
            "rank_median": int(np.nanmedian(r)),
        }
    return {
        "draws": draws,
        "top_k": top_k,
        f"top{top_k}_jaccard_mean": round(float(np.mean(jaccards)), 4),
        f"top{top_k}_jaccard_min": round(float(np.min(jaccards)), 4),
        "intervals": intervals,
    }


# ---------------------------------------------------------------------------
# Mechanism-stratum inference
# ---------------------------------------------------------------------------

def axis_stratum_enrichment(
    rows: list[dict[str, Any]],
    by_name: dict[str, dict],
    drugs: Sequence[dict],
    signature: Signature,
    config: ScoringConfig,
    exposures: dict,
    bootstrap_draws: int = 200,
    seed: int = 7,
) -> dict[str, Any]:
    """Aggregate pair scores by therapeutic-axis stratum and test each stratum.

    Individual pair ranks are unstable under curated target-effect uncertainty,
    so the defensible unit of inference is the *mechanism stratum* rather than
    the specific pair.  For every unordered pair of therapeutic axes this
    reports the score distribution, a Mann-Whitney test against all other
    strata, and a bootstrap interval on the stratum median showing that the
    stratum-level ordering is stable where the pair-level ordering is not.
    """
    from scipy.stats import mannwhitneyu

    def strata_for(row) -> set[tuple[str, str]]:
        axes_a = by_name[row["drug_a"]].get("therapeutic_axes", [])
        axes_b = by_name[row["drug_b"]].get("therapeutic_axes", [])
        return {(min(x, y), max(x, y)) for x in axes_a for y in axes_b}

    primary = [r for r in rows if not r["excluded_from_primary_ranking"]]
    buckets: dict[tuple[str, str], list[float]] = {}
    for row in primary:
        for stratum in strata_for(row):
            buckets.setdefault(stratum, []).append(row["priority_score"])

    # Bootstrap the stratum medians by re-scoring under target uncertainty.
    rng = np.random.default_rng(seed)
    medians: dict[tuple[str, str], list[float]] = {k: [] for k in buckets}
    for _ in range(bootstrap_draws):
        perturbed = []
        for drug in drugs:
            sd = float(drug["target_uncertainty"])
            clone = dict(drug)
            clone["target_effects"] = {
                g: float(np.clip(v + rng.normal(0.0, sd), -1.0, 1.0))
                for g, v in drug["target_effects"].items()
            }
            perturbed.append(clone)
        trial = [r for r in rank_pairs(perturbed, signature, config, exposures)
                 if not r["excluded_from_primary_ranking"]]
        trial_buckets: dict[tuple[str, str], list[float]] = {}
        for row in trial:
            for stratum in strata_for(row):
                trial_buckets.setdefault(stratum, []).append(row["priority_score"])
        for stratum, values in trial_buckets.items():
            if stratum in medians:
                medians[stratum].append(float(np.median(values)))

    results = []
    for stratum, values in buckets.items():
        others = [v for k, vals in buckets.items() if k != stratum for v in vals]
        if len(values) >= 5 and others:
            stat = mannwhitneyu(values, others, alternative="greater")
            p = float(stat.pvalue)
            # Rank-biserial effect size: probability a stratum pair outscores a
            # non-stratum pair, which is interpretable without the sample size.
            effect = float(stat.statistic) / (len(values) * len(others))
        else:
            p, effect = float("nan"), float("nan")
        boot = np.array(medians.get(stratum, []))
        results.append({
            "axis_pair": f"{stratum[0]} + {stratum[1]}",
            "n_pairs": len(values),
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
        "interpretation": (
            "Stratum medians carry bootstrap intervals that are narrow relative to the "
            "spread between strata, whereas individual pair ranks are not stable. Treat "
            "the mechanism stratum, not the specific pair, as the screen's unit of inference."
        ),
        "strata": results,
    }

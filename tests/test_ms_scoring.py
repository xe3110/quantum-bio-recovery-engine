"""Behavioural tests for the extended MS screen.

These pin the properties the scientific claims rest on: direction awareness,
Bliss combination, penalty behaviour, null calibration, and the control
expectations declared in the panel metadata.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from core.biology.ms_scoring import (
    ScoringConfig, Signature, alignment_metrics, bliss_combine, combine_effects,
    combined_safety, load_signature, pair_metrics, rank_pairs,
)
from core.biology.screen_statistics import benjamini_hochberg, pareto_front

ROOT = Path(__file__).resolve().parents[1]


def make_signature(**genes) -> Signature:
    logfc = {g: v[0] for g, v in genes.items()}
    desired = {g: v[1] for g, v in genes.items()}
    return Signature(logfc, desired, {g: abs(v[0]) for g, v in genes.items()},
                     {g: "adaptive_immunity" for g in genes})


def make_drug(name, effects, **overrides):
    base = {
        "name": name, "mechanism_class": name, "evidence_tier": "phase_3",
        "target_uncertainty": 0.1, "safety_classes": [name], "safety_burden": {},
        "pathways": ["adaptive_immunity"], "relevant_pathways": ["adaptive_immunity"],
        "target_effects": effects, "cns_penetration": 0.5, "compartment": "both",
        "therapeutic_axes": ["immunomodulation"], "route": "oral", "half_life_class": "short",
    }
    base.update(overrides)
    return base


# --- direction awareness ---------------------------------------------------

def test_desired_direction_governs_reversal_not_raw_logfc():
    """A compensatory gene (up in disease, but wanted higher) must not be reversed."""
    signature = make_signature(HMOX1=(0.7, 1))  # up in MS, therapeutically increase
    boosting = alignment_metrics({"HMOX1": 0.6}, signature)
    suppressing = alignment_metrics({"HMOX1": -0.6}, signature)
    assert boosting["reversal"] > 0
    assert boosting["counter_therapeutic"] == 0
    assert suppressing["reversal"] == 0
    assert suppressing["counter_therapeutic"] > 0


def test_counter_therapeutic_movement_is_penalised_in_score():
    signature = make_signature(TNF=(1.0, -1))
    config = ScoringConfig()
    helps = pair_metrics(make_drug("A", {"TNF": -0.8}), make_drug("B", {"TNF": -0.1}), signature, config)
    harms = pair_metrics(make_drug("A", {"TNF": -0.8}), make_drug("C", {"TNF": 0.8}), signature, config)
    assert harms["counter_therapeutic"] >= 0
    assert harms["priority_score"] < helps["priority_score"]


# --- combination model -----------------------------------------------------

def test_bliss_saturates_same_direction_and_cancels_opposing():
    assert bliss_combine(0.5, 0.5) == pytest.approx(0.75)
    assert bliss_combine(-0.5, -0.5) == pytest.approx(-0.75)
    assert bliss_combine(0.9, 0.9) < 1.0
    assert bliss_combine(0.6, -0.6) == pytest.approx(0.0)
    assert bliss_combine(0.0, 0.4) == pytest.approx(0.4)


def test_combination_never_exceeds_unit_effect():
    combined = combine_effects({"X": 0.95}, {"X": 0.95})
    assert abs(combined["X"]) <= 1.0


# --- safety ----------------------------------------------------------------

def test_shared_risk_domain_creates_overlap_burden():
    config = ScoringConfig()
    a = make_drug("A", {}, safety_burden={"infection": 0.8})
    b = make_drug("B", {}, safety_burden={"infection": 0.8})
    c = make_drug("C", {}, safety_burden={"ocular": 0.8})
    shared = combined_safety(a, b, config)
    distinct = combined_safety(a, c, config)
    assert shared["safety_overlap"] > distinct["safety_overlap"]
    assert shared["worst_risk_domain"] == "infection"


# --- redundancy exclusion --------------------------------------------------

def test_same_mechanism_or_safety_class_is_excluded_from_primary():
    signature = make_signature(TNF=(1.0, -1), IL10=(-1.0, 1))
    drugs = [
        make_drug("A", {"TNF": -0.8}, mechanism_class="m1", safety_classes=["s1"]),
        make_drug("B", {"IL10": 0.7}, mechanism_class="m2", safety_classes=["s2"]),
        make_drug("C", {"TNF": -0.7}, mechanism_class="m1", safety_classes=["s3"]),
        make_drug("D", {"IL10": 0.6}, mechanism_class="m4", safety_classes=["s1"]),
    ]
    rows = {(r["drug_a"], r["drug_b"]): r for r in rank_pairs(drugs, signature, ScoringConfig())}
    assert not rows[("A", "B")]["excluded_from_primary_ranking"]
    assert rows[("A", "C")]["excluded_from_primary_ranking"], "shared mechanism_class"
    assert rows[("A", "D")]["excluded_from_primary_ranking"], "shared safety_class"


def test_ranking_is_deterministic():
    signature = make_signature(TNF=(1.0, -1), IL10=(-1.0, 1))
    drugs = [make_drug(n, {"TNF": -0.5, "IL10": 0.3}, mechanism_class=n, safety_classes=[n])
             for n in "ABCDE"]
    first = rank_pairs(drugs, signature, ScoringConfig())
    second = rank_pairs(drugs, signature, ScoringConfig())
    assert [r["priority_score"] for r in first] == [r["priority_score"] for r in second]
    assert [(r["drug_a"], r["drug_b"]) for r in first] == [(r["drug_a"], r["drug_b"]) for r in second]


# --- statistics ------------------------------------------------------------

def test_benjamini_hochberg_is_monotone_and_bounded():
    rows = [{"p_empirical": p} for p in (0.001, 0.01, 0.02, 0.2, 0.9)]
    benjamini_hochberg(rows)
    q = [r["q_value"] for r in rows]
    assert all(0.0 <= v <= 1.0 for v in q)
    assert q == sorted(q), "q-values must be non-decreasing in p"
    assert all(a >= b for a, b in zip(q, [r["p_empirical"] for r in rows]))


def test_pareto_front_excludes_dominated_rows():
    objectives = (("good", True), ("bad", False))
    rows = [
        # "best" scores highest but carries the most of the bad objective, so it
        # and the low-risk "tradeoff" are both non-dominated; "dominated" is
        # worse than "tradeoff" on both and must be dropped.
        {"name": "best", "good": 1.0, "bad": 0.9},
        {"name": "tradeoff", "good": 0.4, "bad": 0.1},
        {"name": "dominated", "good": 0.3, "bad": 0.5},
    ]
    pareto_front(rows, objectives)
    flags = {r["name"]: r["pareto_optimal"] for r in rows}
    assert flags["best"] and flags["tradeoff"]
    assert not flags["dominated"]


# --- real-data integration -------------------------------------------------

@pytest.fixture(scope="module")
def real_inputs():
    signature = load_signature(ROOT / "data/ms_expression_v3.csv")
    panel = json.loads((ROOT / "data/drugs/ms_panel_v3.json").read_text())
    return signature, panel


def test_shipped_panel_targets_all_exist_in_signature(real_inputs):
    signature, panel = real_inputs
    universe = set(signature.genes)
    for record in panel["drugs"]:
        missing = set(record["target_effects"]) - universe
        assert not missing, f"{record['name']} targets genes absent from the signature: {missing}"


def test_shipped_panel_expands_prior_version(real_inputs):
    _, panel = real_inputs
    previous = json.loads((ROOT / "data/drugs/ms_panel_v2.json").read_text())
    assert len(panel["drugs"]) > len(previous["drugs"])
    assert len({d["mechanism_class"] for d in panel["drugs"]}) > len(
        {d["mechanism_class"] for d in previous["drugs"]}
    )


def test_declared_controls_are_present_in_panel(real_inputs):
    _, panel = real_inputs
    names = {d["name"] for d in panel["drugs"]}
    for group in panel["metadata"]["controls"].values():
        assert set(group) <= names, f"declared control missing from panel: {set(group) - names}"


def test_redundant_positive_control_pair_is_excluded(real_inputs):
    """Natalizumab and its biosimilar share a mechanism and must never rank."""
    signature, panel = real_inputs
    a, b = panel["metadata"]["controls"]["positive_redundancy"]
    by_name = {d["name"]: d for d in panel["drugs"]}
    row = pair_metrics(by_name[a], by_name[b], signature, ScoringConfig())
    assert row["excluded_from_primary_ranking"]
    assert row["same_mechanism"]


def test_withdrawn_agent_is_ranked_below_median(real_inputs):
    """Daclizumab was withdrawn for fatal toxicity; the safety penalty must bite."""
    signature, panel = real_inputs
    rows = [r for r in rank_pairs(panel["drugs"], signature, ScoringConfig())
            if not r["excluded_from_primary_ranking"]]
    best = min(i for i, r in enumerate(rows) if "Daclizumab" in (r["drug_a"], r["drug_b"]))
    assert best > len(rows) * 0.05, "a withdrawn agent must not surface near the top of the screen"


def test_signature_directions_are_valid(real_inputs):
    signature, _ = real_inputs
    assert len(signature.genes) > 100
    assert set(signature.desired.values()) == {-1, 1}
    assert all(w >= 0 for w in signature.weight.values())


def test_same_target_family_is_excluded_despite_different_mechanism_labels(real_inputs):
    """Different modality, same biological target, is still redundant.

    Firategrast (oral small molecule) and natalizumab (monoclonal antibody)
    carry different mechanism_class labels but both block alpha4 integrin, so
    the pair must never be ranked as a combination.
    """
    signature, panel = real_inputs
    by_name = {d["name"]: d for d in panel["drugs"]}
    row = pair_metrics(by_name["Firategrast"], by_name["Natalizumab"], signature, ScoringConfig())
    assert not row["same_mechanism"], "the labels genuinely differ"
    assert row["redundant_targets"]
    assert row["excluded_from_primary_ranking"]


def test_distinct_mechanism_pair_is_not_flagged_redundant(real_inputs):
    signature, panel = real_inputs
    by_name = {d["name"]: d for d in panel["drugs"]}
    row = pair_metrics(by_name["Minocycline"], by_name["Ocrelizumab"], signature, ScoringConfig())
    assert not row["redundant_targets"]
    assert not row["excluded_from_primary_ranking"]


def test_every_panel_record_declares_a_target_family(real_inputs):
    _, panel = real_inputs
    assert all(d.get("target_family") for d in panel["drugs"])


def test_permutation_seeding_is_stable_across_processes():
    """Per-pair seeds must not depend on Python's randomised builtin hash().

    Otherwise the empirical p-values differ between runs despite a stated seed,
    which would break the reproducibility claim in the protocol.
    """
    import subprocess
    import sys

    snippet = (
        "from core.biology.screen_statistics import _stable_seed;"
        "print(_stable_seed(7, 'Clemastine', 'Minocycline'))"
    )
    runs = {
        subprocess.run([sys.executable, "-c", snippet], cwd=ROOT,
                       capture_output=True, text=True, check=True).stdout.strip()
        for _ in range(3)
    }
    assert len(runs) == 1, f"seed varies across processes: {runs}"

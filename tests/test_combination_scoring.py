"""The k-ary combination screen, held to the claims its docstrings make.

Three of these tests are the ones that matter, because they check properties
the screen's conclusions depend on rather than the arithmetic of any one term:

* k = 2 through :func:`combination_metrics` must equal the pairwise screen's
  own definitions, or the monotherapy-versus-pair comparison is a comparison
  of two implementations;
* the efficacy block must be defined identically at every order, because it is
  the only block compared across orders;
* the vocabulary must come from the registry, because that is the whole reason
  this module exists rather than a second copy of ``ms_scoring``.
"""

from __future__ import annotations

import numpy as np
import pytest

from core.biology.combination_scoring import (
    CombinationConfig, attach_subset_gain, combination_key, combination_metrics,
    combine_many, eligible_drugs, rank_combinations, redundancy_flags,
)
from core.biology.combination_statistics import (
    EFFICACY_FIELDS, benjamini_hochberg, monotherapy_comparison,
)
from core.biology.signature import alignment_metrics, combine_effects
from core.models.disease import load_disease

DISEASES = ["parkinsons", "multiple_sclerosis"]


@pytest.fixture(scope="module")
def pd_disease():
    return load_disease("parkinsons")


@pytest.fixture(scope="module")
def pd_panel(pd_disease):
    return pd_disease.panel()


@pytest.fixture(scope="module")
def pd_signature(pd_disease):
    return pd_disease.signature()


@pytest.fixture(scope="module")
def config():
    return CombinationConfig()


def _named(panel, name):
    return next(d for d in panel if d["name"] == name)


# ---------------------------------------------------------------------------
# The scorer runs for every registered disease, on that disease's vocabulary
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("identifier", DISEASES)
@pytest.mark.parametrize("order", [1, 2, 3])
def test_every_order_scores_for_every_registered_disease(identifier, order, config):
    """The generality claim, at every order the screen offers."""
    disease = load_disease(identifier)
    drugs = eligible_drugs(disease.panel(), config)[:6]
    rows = rank_combinations(drugs, disease.signature(), disease, config, order=order)
    assert rows
    assert all(row["order"] == order for row in rows)
    assert all(len(row["members"]) == order for row in rows)


def test_axis_and_risk_vocabulary_come_from_the_registry(pd_panel, pd_signature, pd_disease, config):
    """Scoring PD against PD's vocabulary must not reach for MS's.

    ``axis_coverage`` is a fraction of the disease's own declared axes, and
    ``worst_risk_domain`` must name one of its own risk domains. If either had
    stayed a module constant, this would fail with an MS label.
    """
    row = combination_metrics(
        [_named(pd_panel, "Levodopa"), _named(pd_panel, "Minocycline")],
        pd_signature, pd_disease, config,
    )
    assert row["worst_risk_domain"] in pd_disease.risk_domains
    assert row["axis_coverage"] == pytest.approx(2 / len(pd_disease.therapeutic_axes))
    ms = load_disease("multiple_sclerosis")
    assert row["worst_risk_domain"] not in ms.risk_domains


# ---------------------------------------------------------------------------
# k = 2 must reduce to the pairwise definitions
# ---------------------------------------------------------------------------

def test_target_complementarity_at_k2_is_the_jaccard_distance(pd_panel, pd_signature, pd_disease, config):
    a, b = _named(pd_panel, "Levodopa"), _named(pd_panel, "Minocycline")
    ga, gb = set(a["target_effects"]), set(b["target_effects"])
    expected = 1.0 - len(ga & gb) / len(ga | gb)
    row = combination_metrics([a, b], pd_signature, pd_disease, config)
    assert row["target_complementarity"] == pytest.approx(expected, abs=1e-4)


def test_safety_union_at_k2_is_the_pairwise_union(pd_panel, pd_signature, pd_disease, config):
    """1 - prod(1 - a_i) must reduce to 1 - (1-a)(1-b) at two members."""
    a, b = _named(pd_panel, "Levodopa"), _named(pd_panel, "Pioglitazone")
    weights = pd_disease.risk_weights
    expected = sum(
        weights[d] * (1.0 - (1.0 - a["safety_burden"].get(d, 0.0)) * (1.0 - b["safety_burden"].get(d, 0.0)))
        for d in pd_disease.risk_domains
    ) / sum(weights.values())
    row = combination_metrics([a, b], pd_signature, pd_disease, config)
    assert row["safety_union"] == pytest.approx(expected, abs=1e-4)


def test_bliss_folding_is_order_independent(pd_panel):
    """A combination is a set; scoring it must not depend on how it was listed."""
    a, b, c = (_named(pd_panel, n) for n in ("Levodopa", "Minocycline", "Deferiprone"))
    effects = [d["target_effects"] for d in (a, b, c)]
    forward = combine_many(effects)
    backward = combine_many(list(reversed(effects)))
    assert set(forward) == set(backward)
    for gene in forward:
        assert forward[gene] == pytest.approx(backward[gene], abs=1e-9)


def test_combine_many_at_two_members_equals_combine_effects(pd_panel):
    a, b = _named(pd_panel, "Levodopa"), _named(pd_panel, "Exenatide")
    pair = combine_effects(a["target_effects"], b["target_effects"])
    assert combine_many([a["target_effects"], b["target_effects"]]) == pair


def test_scoring_a_combination_is_invariant_to_member_order(pd_panel, pd_signature, pd_disease, config):
    members = [_named(pd_panel, n) for n in ("Levodopa", "Minocycline", "Deferiprone")]
    a = combination_metrics(members, pd_signature, pd_disease, config)
    b = combination_metrics(list(reversed(members)), pd_signature, pd_disease, config)
    assert a["priority_score"] == b["priority_score"]
    for field in EFFICACY_FIELDS:
        assert a[field] == b[field]


# ---------------------------------------------------------------------------
# Cross-order comparison: the point of the screen
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("field", EFFICACY_FIELDS)
def test_efficacy_fields_exist_at_every_order(field, pd_panel, pd_signature, pd_disease, config):
    """The comparable block has to be present and finite at k = 1, 2, and 3."""
    members = [_named(pd_panel, n) for n in ("Levodopa", "Minocycline", "Deferiprone")]
    for order in (1, 2, 3):
        row = combination_metrics(members[:order], pd_signature, pd_disease, config)
        assert field in row and np.isfinite(row[field])


def test_a_single_agents_reversal_equals_its_own_alignment(pd_panel, pd_signature, pd_disease, config):
    """Monotherapy must run through the same path, not a special case."""
    drug = _named(pd_panel, "Levodopa")
    row = combination_metrics([drug], pd_signature, pd_disease, config)
    expected = alignment_metrics(drug["target_effects"], pd_signature)
    assert row["signed_reversal"] == pytest.approx(expected["reversal"], abs=1e-4)
    assert row["reversal_gain_over_best_single"] == pytest.approx(0.0, abs=1e-9)
    assert row["additivity_ratio"] == pytest.approx(1.0, abs=1e-4)


def test_a_single_agent_with_no_therapeutic_effect_is_not_called_redundant(pd_panel, pd_signature, pd_disease, config):
    """Carbidopa's only modelled action is counter-therapeutic, so its reversal is 0.

    A ratio of 0/0 reported as 0.0 would rank it alongside a genuinely redundant
    combination. A monotherapy is trivially additive with itself; what carbidopa
    has is no signal, which the reversal field already says.
    """
    row = combination_metrics([_named(pd_panel, "Carbidopa")], pd_signature, pd_disease, config)
    assert row["signed_reversal"] == 0.0
    assert row["counter_therapeutic"] > 0.0
    assert row["additivity_ratio"] == 1.0


def test_disjoint_targets_are_additive_and_shared_targets_are_not(pd_panel, pd_signature, pd_disease, config):
    """``additivity_ratio`` is the field that says whether members overlap.

    Two agents on disjoint genes sum exactly; two agents pushing the same gene
    the same way saturate under Bliss and recover less than their sum. Without
    this, a redundant pair would look as good as a complementary one.
    """
    a, b = _named(pd_panel, "Pramipexole"), _named(pd_panel, "Ropinirole")
    overlap_row = combination_metrics([a, b], pd_signature, pd_disease, config)
    assert set(a["target_effects"]) & set(b["target_effects"])
    assert overlap_row["additivity_ratio"] < 1.0

    disjoint = [_named(pd_panel, "Carbidopa"), _named(pd_panel, "Minocycline")]
    assert not set(disjoint[0]["target_effects"]) & set(disjoint[1]["target_effects"])
    assert combination_metrics(disjoint, pd_signature, pd_disease, config)["additivity_ratio"] == pytest.approx(1.0, abs=1e-4)


def test_monotherapy_comparison_refuses_to_compare_priority_scores(pd_panel, pd_signature, pd_disease, config):
    """The composite includes terms a single agent cannot earn.

    If ``priority_score`` ever appears in the comparable set, a reader will
    conclude combinations beat singles when what they actually did was collect
    complementarity bonuses that are zero by construction at k = 1.
    """
    by_order = {
        order: rank_combinations(pd_panel, pd_signature, pd_disease, config, order=order)
        for order in (1, 2)
    }
    report = monotherapy_comparison(by_order)
    assert "priority_score" not in report["comparable_fields"]
    assert set(report["comparable_fields"]) == set(EFFICACY_FIELDS)
    assert report["by_order"][2]["p_greater_than_monotherapy"] is not None
    assert report["by_order"][1]["p_greater_than_monotherapy"] is None


def test_a_single_agent_scores_zero_on_the_combination_only_terms(pd_panel, pd_signature, pd_disease, config):
    row = combination_metrics([_named(pd_panel, "Levodopa")], pd_signature, pd_disease, config)
    assert row["target_complementarity"] == 0.0
    assert row["safety_overlap"] == 0.0
    assert np.isnan(row["network_separation"])


# ---------------------------------------------------------------------------
# Redundancy, prefiltering, and subset gain
# ---------------------------------------------------------------------------

def test_two_dopamine_agonists_are_excluded_as_redundant(pd_panel, config):
    """The panel's declared positive-redundancy control.

    Pramipexole and ropinirole are both non-ergot D2/D3 agonists. A screen that
    ranks them as a combination is wrong, not surprising.
    """
    flags = redundancy_flags(
        [_named(pd_panel, "Pramipexole"), _named(pd_panel, "Ropinirole")], config
    )
    assert flags["same_mechanism"]
    assert flags["excluded_from_primary_ranking"]


def test_a_triple_containing_a_redundant_pair_is_itself_excluded(pd_panel, config):
    members = [_named(pd_panel, n) for n in ("Pramipexole", "Ropinirole", "Minocycline")]
    assert redundancy_flags(members, config)["excluded_from_primary_ranking"]


def test_the_prefilter_stops_triples_being_built_on_excluded_pairs(pd_panel, pd_signature, pd_disease, config):
    """This is what keeps the order-3 enumeration honest and affordable."""
    pairs = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=2)
    surviving = {combination_key(r["members"]) for r in pairs
                 if not r["excluded_from_primary_ranking"]}
    filtered = rank_combinations(pd_panel, pd_signature, pd_disease, config,
                                 order=3, prefilter=surviving)
    unfiltered = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=3)
    assert 0 < len(filtered) < len(unfiltered)
    for row in filtered:
        for face in (combination_key(f) for f in
                     [row["members"][:2], row["members"][1:], row["members"][::2]]):
            assert face in surviving


def test_subset_gain_measures_what_the_third_agent_adds(pd_panel, pd_signature, pd_disease, config):
    pairs = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=2)
    triples = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=3)[:50]
    attach_subset_gain(triples, pairs)
    lookup = {combination_key(r["members"]): r["priority_score"] for r in pairs}
    for row in triples:
        faces = [combination_key(f) for f in
                 (row["members"][:2], row["members"][1:], row["members"][::2])]
        assert row["best_subset_score"] == pytest.approx(max(lookup[f] for f in faces), abs=1e-4)
        assert row["score_gain_over_best_subset"] == pytest.approx(
            row["priority_score"] - row["best_subset_score"], abs=1e-4)


def test_some_triples_do_not_earn_their_third_agent(pd_panel, pd_signature, pd_disease, config):
    """A screen that finds every addition worthwhile is not measuring anything."""
    pairs = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=2)
    triples = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=3)
    attach_subset_gain(triples, pairs)
    gains = [r["score_gain_over_best_subset"] for r in triples
             if not np.isnan(r["score_gain_over_best_subset"])]
    assert gains
    assert min(gains) < 0.0


# ---------------------------------------------------------------------------
# Determinism and statistics
# ---------------------------------------------------------------------------

def test_ranking_is_deterministic(pd_panel, pd_signature, pd_disease, config):
    a = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=2)
    b = rank_combinations(pd_panel, pd_signature, pd_disease, config, order=2)
    assert [r["combination"] for r in a] == [r["combination"] for r in b]


def test_evidence_floor_excludes_lower_tiers(pd_panel):
    strict = eligible_drugs(pd_panel, CombinationConfig(minimum_evidence="approved"))
    loose = eligible_drugs(pd_panel, CombinationConfig(minimum_evidence="phase_2"))
    assert 0 < len(strict) < len(loose)
    assert all(d["evidence_tier"] == "approved" for d in strict)


def test_benjamini_hochberg_is_monotone_and_bounded():
    rows = [{"p_empirical": p} for p in (0.001, 0.01, 0.04, 0.2, 0.9)]
    benjamini_hochberg(rows)
    q = [r["q_value"] for r in rows]
    assert q == sorted(q)
    assert all(0.0 <= v <= 1.0 for v in q)
    assert all(r["q_value"] >= r["p_empirical"] for r in rows)


def test_empty_combination_is_rejected(pd_signature, pd_disease, config):
    with pytest.raises(ValueError):
        combination_metrics([], pd_signature, pd_disease, config)

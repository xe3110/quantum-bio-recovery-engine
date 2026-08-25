from core.biology.ms_prioritization import PrioritizationConfig, rank_pairs


def test_reversal_and_safety_exclusion_are_explicit():
    signature = {"TNF": 1.0, "IL10": -1.0}
    drugs = [
        {"name": "A", "target_effects": {"TNF": -0.8}, "pathways": ["innate"], "relevant_pathways": ["innate"], "evidence_tier": "phase_3", "target_uncertainty": 0.1, "safety_classes": ["b_cell"], "mechanism_class": "x"},
        {"name": "B", "target_effects": {"IL10": 0.7}, "pathways": ["repair"], "relevant_pathways": ["repair"], "evidence_tier": "phase_3", "target_uncertainty": 0.1, "safety_classes": ["s1p"], "mechanism_class": "y"},
        {"name": "C", "target_effects": {"TNF": -0.7}, "pathways": ["innate"], "relevant_pathways": ["innate"], "evidence_tier": "phase_3", "target_uncertainty": 0.1, "safety_classes": ["b_cell"], "mechanism_class": "z"},
    ]
    rows = rank_pairs(drugs, signature, PrioritizationConfig(minimum_evidence="phase_3"))
    ab = next(row for row in rows if {row["drug_a"], row["drug_b"]} == {"A", "B"})
    ac = next(row for row in rows if {row["drug_a"], row["drug_b"]} == {"A", "C"})
    assert ab["signed_reversal"] > 0
    assert not ab["excluded_from_primary_ranking"]
    assert ac["excluded_from_primary_ranking"]

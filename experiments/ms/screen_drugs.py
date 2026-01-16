from core.biology.ppi_network import PPINetworkBuilder
from core.biology.disease_state import DiseaseStateModel
from core.biology.system_distance import system_distance
from core.chemistry.drug_effects import DrugModel
from core.probability.bayesian_success import BayesianSuccessModel

# ---------------------------
# 1. MS Protein Set
# ---------------------------
MS_PROTEINS = [
    "IL2RA",
    "HLA-DRB1",
    "MBP",
    "TNF",
    "IFNG",
    "CXCL10",
    "CD40",
    "STAT3",
    "MOG",
    "VCAM1"
]

# ---------------------------
# 2. Drug Panel
# ---------------------------
DRUG_PANEL = [
    DrugModel(
        name="Fingolimod",
        targets={
            "IFNG": 0.6,
            "TNF": 0.7,
            "STAT3": 0.75,
            "CXCL10": 0.65,
            "VCAM1": 0.8
        }
    ),
    DrugModel(
        name="Dimethyl Fumarate",
        targets={
            "TNF": 0.75,
            "STAT3": 0.8,
            "IFNG": 0.7,
            "CD40": 0.85
        }
    ),
    DrugModel(
        name="Natalizumab",
        targets={
            "VCAM1": 0.5,
            "CXCL10": 0.7,
            "IL2RA": 0.8
        }
    ),
    DrugModel(
        name="Interferon-beta",
        targets={
            "IFNG": 0.85,
            "STAT3": 0.9,
            "CXCL10": 0.9
        }
    )
]

# ---------------------------
# 3. Build System
# ---------------------------
print("\nBuilding MS System...")
builder = PPINetworkBuilder(score_threshold=0.7)
graph = builder.build_graph(MS_PROTEINS)

model = DiseaseStateModel()
expr = model.load_expression_data("data/ms_expression.csv")

healthy, diseased = model.build_state_vector(graph, expr)
baseline_distance = system_distance(healthy, diseased, graph)

bayes = BayesianSuccessModel(prior_success=2, prior_failure=2)

# ---------------------------
# 4. Screen Drugs
# ---------------------------
results = []

for drug in DRUG_PANEL:
    post_state = drug.apply(diseased)
    post_distance = system_distance(healthy, post_state, graph)

    recovery_score = baseline_distance - post_distance

    alpha, beta_param = bayes.update(
        recovery_score=recovery_score,
        trials=50,
        noise=0.1
    )

    prob = bayes.probability(alpha, beta_param)

    results.append({
        "Drug": drug.name,
        "RecoveryScore": round(recovery_score, 3),
        "MeanProbability": round(prob["mean"], 3),
        "CI_Low": round(prob["ci_low"], 3),
        "CI_High": round(prob["ci_high"], 3)
    })

# ---------------------------
# 5. Rank & Display
# ---------------------------
results = sorted(results, key=lambda x: x["MeanProbability"], reverse=True)

print("\n=== MS Drug Screening Results ===")
for i, r in enumerate(results, 1):
    print(
        f"{i}. {r['Drug']:20} | "
        f"Recovery: {r['RecoveryScore']:>6} | "
        f"P(Success): {r['MeanProbability']} "
        f"[{r['CI_Low']}, {r['CI_High']}]"
    )

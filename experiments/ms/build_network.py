from core.biology.ppi_network import PPINetworkBuilder
from core.biology.disease_state import DiseaseStateModel
from core.biology.system_distance import system_distance
from core.chemistry.drug_effects import DrugModel
from core.probability.bayesian_success import BayesianSuccessModel

import pandas as pd

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

builder = PPINetworkBuilder()
graph = builder.build_graph(MS_PROTEINS)

print("Nodes:", graph.number_of_nodes())
print("Edges:", graph.number_of_edges())

print("\nSample interactions:")
for u, v, data in list(graph.edges(data=True))[:10]:
    print(f"{u} ↔ {v} | confidence={data['weight']:.3f}")

model = DiseaseStateModel()
expr = model.load_expression_data("data/ms_expression.csv")

healthy, diseased = model.build_state_vector(graph, expr)

distance = system_distance(healthy, diseased, graph)

print("\nSystem Disease Distance Score:", round(distance, 4))

fingolimod = DrugModel(
    name="Fingolimod",
    targets={
        "IFNG": 0.6,
        "TNF": 0.7,
        "STAT3": 0.75,
        "CXCL10": 0.65,
        "VCAM1": 0.8
    }
)

post_drug = fingolimod.apply(diseased)

recovery_score = (
    system_distance(healthy, diseased, graph)
    - system_distance(healthy, post_drug, graph)
)

print("\nDrug:", fingolimod.name)
print("Recovery Score:", round(recovery_score, 4))

bayes = BayesianSuccessModel(prior_success=2, prior_failure=2)

alpha, beta_param = bayes.update(recovery_score, trials=30, noise=0.1)
prob = bayes.probability(alpha, beta_param)

print("\nProbability of Therapeutic Success:")
print("Mean:", round(prob["mean"], 3))
print("95% CI:", (round(prob["ci_low"], 3), round(prob["ci_high"], 3)))

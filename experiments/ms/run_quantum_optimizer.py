from core.quantum.qaoa_optimizer import QuantumDrugOptimizer

# Drugs from your panel
drug_names = [
    "Fingolimod",
    "Dimethyl Fumarate",
    "Natalizumab",
    "Interferon-beta"
]

# Recovery scores (from Day 6 results)
recovery_scores = [
    0.629,  # Fingolimod
    0.380,  # Dimethyl Fumarate
    0.237,  # Natalizumab
    0.166   # Interferon-beta
]

# Synergy matrix (from Day 7 results)
synergy_matrix = [
    [0,     0.252, 0.156, 0.104],
    [0.252, 0,     0.237, 0.129],
    [0.156, 0.237, 0,     0.155],
    [0.104, 0.129, 0.155, 0]
]

optimizer = QuantumDrugOptimizer()

qp = optimizer.build_problem(
    drug_names=drug_names,
    recovery_scores=recovery_scores,
    synergy_matrix=synergy_matrix,
    k=3
)

result = optimizer.solve(qp)

print("\n=== Quantum Optimization Result ===")
print("Selected drugs:")
for i, val in enumerate(result.x):
    if round(val) == 1:
        print("-", drug_names[i])

print("Objective value:", round(result.fval, 4))

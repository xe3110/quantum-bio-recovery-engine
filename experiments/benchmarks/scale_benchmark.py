import time
import random
import numpy as np

from core.quantum.qaoa_optimizer import QuantumDrugOptimizer

# ---------------------------
# Configuration
# ---------------------------
SIZES = [8, 12, 16, 20, 24]   # Number of drugs
K = 4                     # Select k drugs (k >= 3 enters hard regime)
RANDOM_TRIALS = 500     # Random search attempts per N


# ---------------------------
# Problem Generator
# ---------------------------
def generate_problem(n):
    """
    Generate synthetic recovery + synergy landscape
    """
    recovery = np.random.uniform(0.1, 1.0, n)
    synergy = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            synergy[i][j] = synergy[j][i] = np.random.uniform(0, 0.5)

    return recovery, synergy


# ---------------------------
# Objective Function
# ---------------------------
def score_set(choice, recovery, synergy):
    """
    Compute:
    sum(recovery[i]) + sum(synergy[i][j]) for all i<j in choice
    """
    s = 0.0

    for i in choice:
        s += recovery[i]

    for i in range(len(choice)):
        for j in range(i + 1, len(choice)):
            s += synergy[choice[i]][choice[j]]

    return s


# ---------------------------
# Greedy Solver (k-selection)
# ---------------------------
def greedy_solver(recovery, synergy, k):
    n = len(recovery)
    selected = []

    # Pick best first drug
    first = max(range(n), key=lambda i: recovery[i])
    selected.append(first)

    while len(selected) < k:
        best = None
        best_score = -1e9

        for j in range(n):
            if j in selected:
                continue

            candidate = selected + [j]
            score = score_set(candidate, recovery, synergy)

            if score > best_score:
                best_score = score
                best = j

        selected.append(best)

    return selected


# ---------------------------
# Random Solver (k-selection)
# ---------------------------
def random_solver(recovery, synergy, k, trials):
    n = len(recovery)
    best_score = -1e9
    best_choice = None

    for _ in range(trials):
        choice = random.sample(range(n), k)
        score = score_set(choice, recovery, synergy)

        if score > best_score:
            best_score = score
            best_choice = choice

    return best_choice, best_score


# ---------------------------
# Benchmark Runner
# ---------------------------
def run_benchmark():
    print("\n=== Scaling Benchmark (k = {}) ===".format(K))
    results = []

    for n in SIZES:
        print(f"\n--- Benchmark N = {n} drugs ---")
        recovery, synergy = generate_problem(n)
        drug_names = [f"D{i}" for i in range(n)]

        # -------------------
        # Exact Solver (Hamiltonian Baseline)
        # -------------------
        optimizer = QuantumDrugOptimizer()
        qp = optimizer.build_problem(
            drug_names=drug_names,
            recovery_scores=recovery.tolist(),
            synergy_matrix=synergy.tolist(),
            k=K
        )

        t0 = time.time()
        result = optimizer.solve(qp)
        exact_time = time.time() - t0
        exact_value = result.fval

        # -------------------
        # Greedy Solver
        # -------------------
        t0 = time.time()
        gsel = greedy_solver(recovery, synergy, K)
        greedy_value = score_set(gsel, recovery, synergy)
        greedy_time = time.time() - t0

        # -------------------
        # Random Solver
        # -------------------
        t0 = time.time()
        rsel, random_value = random_solver(
            recovery, synergy, K, RANDOM_TRIALS
        )
        random_time = time.time() - t0

        print(f"Exact   | value={exact_value:.3f} | time={exact_time:.4f}s")
        print(f"Greedy | value={greedy_value:.3f} | time={greedy_time:.4f}s")
        print(f"Random | value={random_value:.3f} | time={random_time:.4f}s")

        results.append({
            "N": n,
            "ExactValue": exact_value,
            "ExactTime": exact_time,
            "GreedyValue": greedy_value,
            "GreedyTime": greedy_time,
            "RandomValue": random_value,
            "RandomTime": random_time
        })

    return results


# ---------------------------
# Entry Point
# ---------------------------
if __name__ == "__main__":
    run_benchmark()

import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("experiments/benchmarks/results_k4.csv")

N = df["N"]

# ---------------------------
# Runtime Plot
# ---------------------------
plt.figure()
plt.plot(N, df["ExactTime"], marker="o", label="Exact Solver")
plt.plot(N, [0]*len(N), marker="x", label="Greedy (≈0)")
plt.plot(N, [0]*len(N), marker="^", label="Random (≈0)")

plt.yscale("log")
plt.xlabel("Number of Drugs (N)")
plt.ylabel("Runtime (seconds, log scale)")
plt.title("Runtime Scaling for k = 4 Drug Selection")
plt.legend()
plt.grid(True)
plt.savefig("fig_runtime.png", dpi=300)

# ---------------------------
# Optimality Gap Plot
# ---------------------------
greedy_gap = df["ExactValue"] - df["GreedyValue"]
random_gap = df["ExactValue"] - df["RandomValue"]

plt.figure()
plt.plot(N, greedy_gap, marker="o", label="Greedy Gap")
plt.plot(N, random_gap, marker="^", label="Random Gap")

plt.xlabel("Number of Drugs (N)")
plt.ylabel("Optimality Gap (Exact − Heuristic)")
plt.title("Solution Quality Gap for k = 4 Drug Selection")
plt.legend()
plt.grid(True)
plt.savefig("fig_quality.png", dpi=300)

print("Saved: fig_runtime.png, fig_quality.png")

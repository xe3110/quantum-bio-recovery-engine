from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit_algorithms.minimum_eigensolvers import NumPyMinimumEigensolver


class QuantumDrugOptimizer:
    """
    Quantum-ready optimizer.
    Currently uses a classical exact eigensolver to solve
    the same Hamiltonian that would be sent to QAOA / VQE.

    This is the standard benchmark baseline used in
    quantum optimization research.
    """

    def __init__(self):
        self.solver = NumPyMinimumEigensolver()

    def build_problem(self, drug_names, recovery_scores, synergy_matrix, k=3):
        qp = QuadraticProgram()

        n = len(drug_names)
        for i in range(n):
            qp.binary_var(name=f"x{i}")

        # Objective: maximize recovery + synergy
        linear = {f"x{i}": float(recovery_scores[i]) for i in range(n)}

        quadratic = {}
        for i in range(n):
            for j in range(i + 1, n):
                quadratic[(f"x{i}", f"x{j}")] = float(synergy_matrix[i][j])

        qp.maximize(linear=linear, quadratic=quadratic)

        # Hard constraint: pick exactly k drugs
        qp.linear_constraint(
            linear={f"x{i}": 1 for i in range(n)},
            sense="==",
            rhs=k,
            name="pick_k"
        )

        return qp



    def solve(self, qp):
        optimizer = MinimumEigenOptimizer(self.solver)
        result = optimizer.solve(qp)
        return result

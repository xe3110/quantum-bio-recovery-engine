import copy
import math


class DrugModel:
    def __init__(self, name, targets):
        """
        targets:
          dict[str, float]
          protein -> effect strength (0–1)
        """
        self.name = name
        self.targets = targets or {}

    def _diffuse(self, state, graph, alpha=0.6, steps=3):
        """
        Network diffusion (random walk with restart)
        state: dict[node, value]
        graph: dict[node, dict[neighbor, weight]]
        alpha: retention factor (0–1)
        steps: propagation steps
        """
        new_state = copy.deepcopy(state)

        for _ in range(steps):
            propagated = {}

            for node, value in new_state.items():
                neighbors = graph.get(node, {})
                if not neighbors:
                    propagated[node] = value
                    continue

                influence = 0.0
                total_w = 0.0

                for nbr, w in neighbors.items():
                    influence += w * new_state.get(nbr, value)
                    total_w += w

                avg_influence = influence / total_w if total_w > 0 else value

                # Random walk with restart
                propagated[node] = alpha * value + (1 - alpha) * avg_influence

            new_state = propagated

        return new_state

    def apply(self, diseased_state, graph=None):
        """
        Apply drug effect + propagate through network
        """
        state = copy.deepcopy(diseased_state)

        # Direct target effects
        for protein, strength in self.targets.items():
            if protein in state:
                # Pull protein back toward healthy baseline (1.0)
                state[protein] = state[protein] + strength * (1.0 - state[protein])

        # Network propagation
        if graph:
            state = self._diffuse(state, graph)

        return state

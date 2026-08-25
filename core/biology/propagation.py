import copy
import math


def normalize_adjacency(ppi):
    """
    Normalize adjacency weights so each node's outgoing edges sum to 1.
    ppi: dict[str, dict[str, float]]
         {node: {neighbor: confidence}}
    """
    norm = {}
    for node, neighbors in ppi.items():
        total = sum(neighbors.values())
        if total == 0:
            continue
        norm[node] = {
            nbr: w / total
            for nbr, w in neighbors.items()
        }
    return norm


def propagate_state(
    initial_state,
    ppi_network,
    alpha=0.4,
    steps=8
):
    """
    Diffuse protein activity across the PPI network.

    initial_state: dict[str, float]
    ppi_network: dict[str, dict[str, float]]
    Returns: dict[str, float]
    """
    x0 = copy.deepcopy(initial_state)
    x = copy.deepcopy(initial_state)

    A = normalize_adjacency(ppi_network)

    for _ in range(steps):
        x_next = copy.deepcopy(x0)

        for node, neighbors in A.items():
            influence = 0.0
            for nbr, w in neighbors.items():
                if nbr in x:
                    influence += w * x[nbr]

            if node in x_next:
                x_next[node] = (
                    alpha * influence +
                    (1 - alpha) * x0.get(node, 0.0)
                )

        x = x_next

    return x

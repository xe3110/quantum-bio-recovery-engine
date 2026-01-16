import numpy as np

def system_distance(state_a, state_b, graph):
    """
    Network-aware distance:
    Weigh node differences by interaction strength
    """
    diffs = []

    for u, v, data in graph.edges(data=True):
        w = data["weight"]
        da = abs(state_a[u] - state_b[u])
        db = abs(state_a[v] - state_b[v])
        diffs.append(w * (da + db) / 2)

    return np.mean(diffs)

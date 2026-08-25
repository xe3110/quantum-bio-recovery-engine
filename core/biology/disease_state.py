import random
import math


class DiseaseStateModel:
    def __init__(self, ppi_network):
        """
        ppi_network:
          - networkx.Graph
          OR
          - dict[str, dict[str, float]]

        Internally normalized to:
          self.graph = {node: {neighbor: weight}}
        """
        self.network = ppi_network

        # ----------------------------
        # Normalize graph
        # ----------------------------
        if hasattr(ppi_network, "nodes"):
            # NetworkX graph
            self.graph = {
                n: {
                    nbr: ppi_network[n][nbr].get("weight", 1.0)
                    for nbr in ppi_network.neighbors(n)
                }
                for n in ppi_network.nodes
            }
        else:
            # Already adjacency dict
            self.graph = ppi_network

        # ----------------------------
        # Master regulator weights
        # ----------------------------
        # These nodes contribute more to disease recovery
        MASTER_NODES = {
            "SNCA",
            "LRRK2",
            "TH",
            "MAOB",
            "IFNG",
            "TNF",
            "STAT3"
        }

        self.node_weights = {
            node: 2.0 if node in MASTER_NODES else 1.0
            for node in self.graph
        }

    # ============================
    # STATE GENERATORS
    # ============================

    def healthy_state(self):
        """
        Baseline healthy protein activity state

        Returns:
          dict[str, float]  (node → activity level)
        """
        return {node: 1.0 for node in self.graph}

    def diseased_state(self, severity=0.7):
        """
        Structured Parkinson's disease signature
        """
        # Known directional biology
        UP = {"SNCA", "LRRK2", "MAOB", "COMT"}
        DOWN = {"TH", "SLC6A3", "PINK1", "PRKN"}

        state = {}

        for node in self.graph:
            base = 1.0

            if node in UP:
                val = base + severity
            elif node in DOWN:
                val = base - severity
            else:
                # mild random drift for non-core proteins
                val = base + random.uniform(-0.1, 0.1)

            state[node] = max(0.1, min(2.0, val))

        return state


    # ============================
    # SCORING
    # ============================

    def recovery_score(self, healthy, post_state):
        """
        Weighted comparison of post-drug state to healthy baseline

        Returns:
          float (0–1, higher = better recovery)
        """
        total = 0.0

        for node in healthy:
            h = healthy[node]
            p = post_state.get(node, h)

            w = self.node_weights.get(node, 1.0)
            total += w * abs(h - p)

        max_dist = sum(self.node_weights.values())

        if max_dist == 0:
            return 0.0

        return round(1.0 - (total / max_dist), 3)
    
    def disease_signature(self, healthy, diseased):
        """
        Vector of disease-induced changes per node
        """
        return {
            node: diseased[node] - healthy[node]
            for node in healthy
        }
    
    def directional_recovery(self, healthy, diseased, post_state):
        disease_vec = []
        drug_vec = []

        for node in healthy:
            d = diseased[node] - healthy[node]
            r = post_state[node] - diseased[node]

            w = self.node_weights.get(node, 1.0)

            disease_vec.append(w * d)
            drug_vec.append(w * r)

        dot = sum(d * r for d, r in zip(disease_vec, drug_vec))
        mag_d = math.sqrt(sum(d * d for d in disease_vec))
        mag_r = math.sqrt(sum(r * r for r in drug_vec))

        if mag_d == 0 or mag_r == 0:
            return 0.0

        cos_sim = dot / (mag_d * mag_r)

        # -------------------------
        # Pharmacological response curve
        # -------------------------
        # Raw reversal signal
        reversal = -cos_sim  # higher = better

        # Logistic calibration
        k = 4.0      # steepness (increase to spread more)
        x0 = 0.2     # midpoint (shift to control difficulty)

        score = 1.0 / (1.0 + math.exp(-k * (reversal - x0)))

        return round(score, 3)



    def success_probability(self, recovery, sigma=0.15):
        """
        Convert recovery score into probability + confidence interval

        Returns:
          (mean, (low, high))
        """
        mean = max(0.0, min(1.0, recovery))

        low = max(0.0, mean - sigma)
        high = min(1.0, mean + sigma)

        return round(mean, 3), (round(low, 3), round(high, 3))

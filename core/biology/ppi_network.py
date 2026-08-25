import requests
import networkx as nx

STRING_API_URL = "https://string-db.org/api/json/network"


class PPINetworkBuilder:
    def __init__(self, species=9606, score_threshold=0.4):
        """
        species: NCBI species ID (9606 = Homo sapiens)
        score_threshold: Minimum interaction confidence (0–1)
        """
        self.species = species
        self.score_threshold = score_threshold

        self.graph = None
        self.ppi = None

    def fetch_interactions(self, proteins):
        params = {
            "identifiers": "%0d".join(proteins),
            "species": self.species,
            "required_score": int(self.score_threshold * 1000)
        }

        response = requests.get(STRING_API_URL, params=params, timeout=15)
        response.raise_for_status()
        return response.json()

    def build_graph(self, proteins):
        print("🌐 Fetching PPI network from STRING...")

        data = self.fetch_interactions(proteins)
        G = nx.Graph()

        edge_count = 0

        for item in data:
            p1 = item["preferredName_A"]
            p2 = item["preferredName_B"]
            score = item["score"]

            G.add_node(p1)
            G.add_node(p2)

            if score >= self.score_threshold:
                G.add_edge(p1, p2, weight=score)
                edge_count += 1

        self.graph = G
        self.ppi = self._to_adjacency_dict(G)

        print(f"Nodes: {G.number_of_nodes()}")
        print(f"Edges: {edge_count}")

        self._print_sample_edges(G)

        return G

    def _to_adjacency_dict(self, G):
        """
        Convert NetworkX graph to:
        {node: {neighbor: weight}}
        """
        ppi = {}

        for node in G.nodes:
            ppi[node] = {}

        for u, v, data in G.edges(data=True):
            w = float(data.get("weight", 1.0))
            ppi[u][v] = w
            ppi[v][u] = w

        return ppi

    def _print_sample_edges(self, G, n=10):
        """
        Print a few sample PPI interactions for debugging/inspection
        """
        print("\nSample interactions:")
        for i, (u, v, data) in enumerate(G.edges(data=True)):
            if i >= n:
                break
            print(f"{u} ↔ {v} | confidence={data.get('weight', 0):.3f}")
        print()

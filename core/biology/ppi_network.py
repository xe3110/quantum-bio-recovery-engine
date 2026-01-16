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

    def fetch_interactions(self, proteins):
        params = {
            "identifiers": "%0d".join(proteins),
            "species": self.species,
            "required_score": int(self.score_threshold * 1000)
        }

        response = requests.get(STRING_API_URL, params=params)
        response.raise_for_status()
        return response.json()

    def build_graph(self, proteins):
        data = self.fetch_interactions(proteins)
        G = nx.Graph()

        for item in data:
            p1 = item["preferredName_A"]
            p2 = item["preferredName_B"]
            score = item["score"]

            G.add_node(p1)
            G.add_node(p2)
            G.add_edge(p1, p2, weight=score)

        return G

import pandas as pd
import numpy as np

class DiseaseStateModel:
    def __init__(self, fold_change_col="logFC"):
        self.fold_change_col = fold_change_col

    def load_expression_data(self, csv_path):
        """
        CSV should contain:
        - 'gene' column
        - fold change column (default: logFC)
        """
        return pd.read_csv(csv_path)

    def build_state_vector(self, graph, expression_df):
        """
        Returns:
        - healthy_state: dict {protein: 1.0}
        - disease_state: dict {protein: weighted by expression}
        """
        healthy_state = {}
        disease_state = {}

        for node in graph.nodes():
            healthy_state[node] = 1.0

            match = expression_df[expression_df["gene"] == node]
            if not match.empty:
                fc = match[self.fold_change_col].values[0]
                disease_state[node] = float(np.exp(fc))  # log fold-change → linear
            else:
                disease_state[node] = 1.0  # unknown = neutral

        return healthy_state, disease_state

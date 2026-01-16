class DrugModel:
    def __init__(self, name, targets):
        """
        targets: dict {protein: effect}
        effect:
          <1.0 = suppression
          >1.0 = activation
        """
        self.name = name
        self.targets = targets

    def apply(self, disease_state):
        post_state = disease_state.copy()

        for protein, effect in self.targets.items():
            if protein in post_state:
                post_state[protein] *= effect

        return post_state

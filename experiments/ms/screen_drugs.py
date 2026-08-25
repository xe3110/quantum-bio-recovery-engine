import argparse
from core.biology.ppi_network import PPINetworkBuilder
from core.biology.disease_state import DiseaseStateModel
from core.drugs.chembl_loader import fetch_drugs_for_targets
from core.biology.target_loader import load_targets

def parse_args():
    parser = argparse.ArgumentParser(
        description="Quantum Bio Recovery Engine — Drug Screening CLI"
    )

    parser.add_argument(
        "--disease",
        type=str,
        default="ms",
        help="Disease name (loads data/seeds/<disease>_drugs.json)"
    )

    parser.add_argument(
        "--mode",
        type=str,
        default="therapeutic",
        choices=["binding", "therapeutic"],
        help="Drug discovery mode"
    )

    parser.add_argument(
        "--max-drugs",
        type=int,
        default=50,
        help="Maximum number of drugs to screen"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    print("\nBuilding Disease System...")
    print(f"Disease: {args.disease}")
    print(f"Mode: {args.mode}")
    print(f"Max Drugs: {args.max_drugs}\n")

    TARGETS = load_targets(args.disease)

    if not TARGETS:
        print("❌ No disease targets found. Add data/targets/<disease>.json")
        return

    # ----------------------------
    # Build PPI Network
    # ----------------------------
    builder = PPINetworkBuilder()
    G = builder.build_graph(TARGETS)

    print(f"\nPPI Network Loaded: {len(G.nodes)} nodes, {len(G.edges)} edges\n")

    # ----------------------------
    # Disease State Model
    # ----------------------------
    disease_model = DiseaseStateModel(G)
    healthy = disease_model.healthy_state()
    diseased = disease_model.diseased_state()

    # ----------------------------
    # Load Drug Panel
    # ----------------------------
    DRUG_PANEL = fetch_drugs_for_targets(
        targets=TARGETS,
        disease=args.disease,
        max_drugs=args.max_drugs,
        mode=args.mode
    )

    # ----------------------------
    # Screen Drugs
    # ----------------------------
    print("\n=== Drug Screening Results ===\n")

    results = []
    for drug in DRUG_PANEL:
        post_state = drug.apply(diseased, graph=disease_model.graph)
        recovery = disease_model.directional_recovery(
            healthy,
            diseased,
            post_state
        )
        prob, ci = disease_model.success_probability(recovery)

        results.append((drug.name, recovery, prob, ci))

    results.sort(key=lambda x: x[1], reverse=True)

    for i, (name, recovery, prob, ci) in enumerate(results[:10], 1):
        print(
            f"{i:2d}. {name:25s} | "
            f"Recovery: {recovery:0.3f} | "
            f"P(Success): {prob:0.3f} {ci}"
        )

    print("\nDone.\n")


if __name__ == "__main__":
    main()

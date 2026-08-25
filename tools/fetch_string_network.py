"""Fetch and cache the STRING PPI subnetwork used for network-proximity scoring.

The cached edge list is committed so the screen is reproducible offline and so a
reviewer can inspect the exact interactome the results depend on.

Run: python -m tools.fetch_string_network [--score 400] [--add-nodes 150]
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
STRING_URL = "https://string-db.org/api/json/network"
OUT_EDGES = ROOT / "data/networks/string_ms_network.tsv"
OUT_META = ROOT / "data/networks/string_ms_network.meta.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score", type=int, default=400, help="STRING confidence x1000")
    parser.add_argument("--add-nodes", type=int, default=150,
                        help="STRING first-shell expansion, for shortest-path connectivity")
    parser.add_argument("--species", type=int, default=9606)
    args = parser.parse_args()

    with (ROOT / "data/ms_expression_v3.csv").open(newline="") as handle:
        genes = [row["gene"] for row in csv.DictReader(handle)]

    response = requests.get(
        STRING_URL,
        params={
            "identifiers": "%0d".join(genes),
            "species": args.species,
            "required_score": args.score,
            "add_nodes": args.add_nodes,
            "caller_identity": "quantum-bio-recovery-engine",
        },
        timeout=120,
    )
    response.raise_for_status()
    records = response.json()

    edges = set()
    for item in records:
        a, b = item["preferredName_A"], item["preferredName_B"]
        if a == b:
            continue
        edges.add((min(a, b), max(a, b), round(float(item["score"]), 3)))

    OUT_EDGES.parent.mkdir(parents=True, exist_ok=True)
    with OUT_EDGES.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_a", "protein_b", "combined_score"])
        writer.writerows(sorted(edges))

    nodes = {n for a, b, _ in edges for n in (a, b)}
    missing = sorted(set(genes) - nodes)
    meta = {
        "source": "STRING v12 REST API",
        "url": STRING_URL,
        "species": args.species,
        "required_score": args.score,
        "add_nodes": args.add_nodes,
        "query_genes": len(genes),
        "nodes": len(nodes),
        "edges": len(edges),
        "signature_genes_in_network": len(set(genes) & nodes),
        "signature_genes_absent": missing,
        "note": (
            "Absent genes are excluded from network-proximity metrics but still "
            "contribute to signature-reversal metrics."
        ),
    }
    OUT_META.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()

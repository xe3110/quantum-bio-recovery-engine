"""Fetch and cache a STRING PPI subnetwork for network-proximity scoring.

The cached edge list is committed so the screen is reproducible offline and so a
reviewer can inspect the exact interactome the results depend on.

Disease-agnostic: point it at any signature CSV and output path. It defaults to
the MS signature so existing invocations keep working.

Run: python -m tools.fetch_string_network [--score 400] [--add-nodes 150]
     python -m tools.fetch_string_network --signature data/pd_expression.csv \
         --out data/networks/string_pd_network.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
STRING_URL = "https://string-db.org/api/json/network"
DEFAULT_SIGNATURE = ROOT / "data/ms_expression_v3.csv"
DEFAULT_EDGES = ROOT / "data/networks/string_ms_network.tsv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score", type=int, default=400, help="STRING confidence x1000")
    parser.add_argument("--add-nodes", type=int, default=150,
                        help="STRING first-shell expansion, for shortest-path connectivity")
    parser.add_argument("--species", type=int, default=9606)
    parser.add_argument("--signature", type=Path, default=DEFAULT_SIGNATURE,
                        help="signature CSV whose genes seed the query")
    parser.add_argument("--out", type=Path, default=DEFAULT_EDGES,
                        help="edge-list TSV to write; metadata goes alongside it")
    args = parser.parse_args()

    out_edges = args.out if args.out.is_absolute() else ROOT / args.out
    out_meta = out_edges.with_suffix(".meta.json")
    signature_path = args.signature if args.signature.is_absolute() else ROOT / args.signature

    with signature_path.open(newline="") as handle:
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

    out_edges.parent.mkdir(parents=True, exist_ok=True)
    with out_edges.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_a", "protein_b", "combined_score"])
        writer.writerows(sorted(edges))

    nodes = {n for a, b, _ in edges for n in (a, b)}
    missing = sorted(set(genes) - nodes)
    meta = {
        "source": "STRING v12 REST API",
        "signature": str(signature_path.relative_to(ROOT)),
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
    out_meta.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()

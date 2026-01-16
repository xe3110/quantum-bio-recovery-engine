from core.biology.ppi_network import PPINetworkBuilder

MS_PROTEINS = [
    "IL2RA",
    "HLA-DRB1",
    "MBP",
    "TNF",
    "IFNG",
    "CXCL10",
    "CD40",
    "STAT3",
    "MOG",
    "VCAM1"
]

builder = PPINetworkBuilder()
graph = builder.build_graph(MS_PROTEINS)

print("Nodes:", graph.number_of_nodes())
print("Edges:", graph.number_of_edges())

print("\nSample interactions:")
for u, v, data in list(graph.edges(data=True))[:10]:
    print(f"{u} ↔ {v} | confidence={data['weight']:.3f}")

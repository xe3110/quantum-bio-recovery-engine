"""Network-based proximity and separation metrics on a cached STRING interactome.

Implements the drug-drug separation and disease-module proximity measures used
to define the *complementary exposure* pattern for drug combinations
(Menche et al., Science 2015; Cheng et al., Nat Commun 2019, PMID: 31000720).

Two agents show complementary exposure when each target set is individually
close to the disease module while the two target sets are separated from one
another (s_AB > 0).  This is an interactome-topology criterion and is
independent of the expression-reversal score, so it supplies genuinely
orthogonal evidence rather than a restatement of the same signal.
"""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from itertools import combinations
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EDGES = ROOT / "data/networks/string_ms_network.tsv"


def load_network(path: Path | str = DEFAULT_EDGES, min_score: float = 0.4) -> nx.Graph:
    """Load the cached STRING edge list as an unweighted graph.

    Shortest-path topology, not edge confidence, drives proximity; ``min_score``
    filters which interactions are admitted as edges at all.
    """
    graph = nx.Graph()
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row["combined_score"]) >= min_score:
                graph.add_edge(row["protein_a"], row["protein_b"])
    if graph.number_of_nodes() == 0:
        raise ValueError(f"No edges loaded from {path} at min_score={min_score}")
    return graph


@lru_cache(maxsize=8)
def _degree_bins(graph_key: int, graph: nx.Graph) -> dict[int, list[str]]:
    """Bin nodes by degree so null models preserve the degree distribution.

    Bins are widened until each holds at least 20 nodes, following the binning
    strategy used for interactome proximity z-scores.
    """
    by_degree = defaultdict(list)
    for node, degree in graph.degree():
        by_degree[degree].append(node)

    bins: dict[int, list[str]] = {}
    bucket: list[str] = []
    members: list[int] = []
    for degree in sorted(by_degree):
        bucket.extend(by_degree[degree])
        members.append(degree)
        if len(bucket) >= 20:
            for d in members:
                bins[d] = list(bucket)
            bucket, members = [], []
    if bucket:  # fold the tail into the last complete bin
        leftover = list(bucket)
        for d in members:
            bins[d] = leftover
        if len(leftover) < 20:
            tail = sorted(graph.degree(), key=lambda kv: -kv[1])[:20]
            pool = leftover + [n for n, _ in tail]
            for d in members:
                bins[d] = pool
    return bins


class ProximityScorer:
    """Compute closest-distance proximity and separation with a degree-matched null."""

    def __init__(self, graph: nx.Graph, seed: int = 7, permutations: int = 200):
        self.graph = graph
        self.seed = seed
        self.permutations = permutations
        self._paths: dict[str, dict[str, int]] = {}
        self._z_cache: dict[tuple, tuple[float, float]] = {}
        self._bins = _degree_bins(id(graph), graph)

    def _distances_from(self, source: str) -> dict[str, int]:
        if source not in self._paths:
            self._paths[source] = nx.single_source_shortest_path_length(self.graph, source)
        return self._paths[source]

    def in_network(self, genes) -> list[str]:
        return sorted({g for g in genes if g in self.graph})

    def closest_distance(self, set_a, set_b) -> float:
        """Mean over ``set_a`` of the shortest distance to any node in ``set_b``."""
        a, b = self.in_network(set_a), self.in_network(set_b)
        if not a or not b:
            return float("nan")
        totals = []
        for source in a:
            dist = self._distances_from(source)
            reachable = [dist[t] for t in b if t in dist]
            if reachable:
                totals.append(min(reachable))
        return float(np.mean(totals)) if totals else float("nan")

    def _within_distance(self, genes) -> float:
        """Mean nearest-neighbour distance inside one set (0 for singletons)."""
        nodes = self.in_network(genes)
        if len(nodes) < 2:
            return 0.0
        totals = []
        for source in nodes:
            dist = self._distances_from(source)
            others = [dist[t] for t in nodes if t != source and t in dist]
            if others:
                totals.append(min(others))
        return float(np.mean(totals)) if totals else 0.0

    def separation(self, set_a, set_b) -> float:
        """Menche separation s_AB = <d_AB> - (<d_AA> + <d_BB>)/2.

        s_AB > 0 means the two target modules are topologically separated;
        s_AB < 0 means they overlap.
        """
        d_ab = self.closest_distance(set_a, set_b)
        d_ba = self.closest_distance(set_b, set_a)
        if np.isnan(d_ab) or np.isnan(d_ba):
            return float("nan")
        mean_between = (d_ab + d_ba) / 2.0
        return float(mean_between - (self._within_distance(set_a) + self._within_distance(set_b)) / 2.0)

    def _matched_sample(self, genes, rng: random.Random) -> list[str]:
        sample = []
        for gene in genes:
            pool = self._bins.get(self.graph.degree(gene), list(self.graph.nodes))
            sample.append(rng.choice(pool))
        return sample

    def proximity_z(self, targets, disease_genes) -> tuple[float, float]:
        """Return (observed closest distance, z-score vs a degree-matched null).

        A significantly negative z indicates the target set sits closer to the
        disease module than degree-matched random protein sets.
        """
        targets = self.in_network(targets)
        disease = self.in_network(disease_genes)
        if not targets or not disease:
            return float("nan"), float("nan")
        # Only a handful of distinct target sets exist across thousands of
        # pairs, so the degree-matched null is computed once per set.
        cache_key = (tuple(targets), len(disease))
        if cache_key in self._z_cache:
            return self._z_cache[cache_key]
        observed = self.closest_distance(targets, disease)
        rng = random.Random(f"{self.seed}:{','.join(targets)}")
        null = []
        for _ in range(self.permutations):
            null.append(self.closest_distance(self._matched_sample(targets, rng), disease))
        null_arr = np.array([v for v in null if not np.isnan(v)])
        if null_arr.size < 2 or null_arr.std() == 0:
            result = (round(observed, 4), 0.0)
        else:
            z = (observed - null_arr.mean()) / null_arr.std()
            result = (round(observed, 4), round(float(z), 4))
        self._z_cache[cache_key] = result
        return result


@dataclass(frozen=True)
class ComplementaryExposure:
    separation: float
    z_a: float
    z_b: float
    is_complementary: bool


def classify_exposure(scorer: ProximityScorer, targets_a, targets_b, disease_genes,
                      z_threshold: float = -0.15) -> ComplementaryExposure:
    """Label a pair by the complementary-exposure pattern.

    Both target sets must be proximal to the disease module (z below
    ``z_threshold``) while remaining topologically separated from each other.
    """
    _, z_a = scorer.proximity_z(targets_a, disease_genes)
    _, z_b = scorer.proximity_z(targets_b, disease_genes)
    sep = scorer.separation(targets_a, targets_b)
    complementary = bool(
        not np.isnan(sep) and sep > 0
        and not np.isnan(z_a) and z_a < z_threshold
        and not np.isnan(z_b) and z_b < z_threshold
    )
    return ComplementaryExposure(
        separation=round(sep, 4) if not np.isnan(sep) else float("nan"),
        z_a=z_a, z_b=z_b, is_complementary=complementary,
    )


def classify_set_exposure(scorer: ProximityScorer, target_sets, disease_genes,
                          z_threshold: float = -0.15) -> ComplementaryExposure:
    """The k-ary generalisation of :func:`classify_exposure`.

    A combination shows complementary exposure when *every* member sits close
    to the disease module while *every* pair of members remains topologically
    separated from the others. Requiring it of every pair rather than on
    average is the conservative reading: a triple in which two agents overlap
    is not exploiting three distinct neighbourhoods, whatever the mean says.

    ``separation`` is reported as the minimum pairwise separation, which is the
    binding one; ``z_a``/``z_b`` carry the best and worst member proximity.
    """
    sets = [list(s) for s in target_sets]
    zs = [scorer.proximity_z(s, disease_genes)[1] for s in sets]
    if len(sets) < 2:
        z = zs[0] if zs else float("nan")
        return ComplementaryExposure(separation=float("nan"), z_a=z, z_b=z, is_complementary=False)

    seps = [scorer.separation(a, b) for a, b in combinations(sets, 2)]
    finite = [s for s in seps if not np.isnan(s)]
    worst_sep = min(finite) if finite else float("nan")
    finite_z = [z for z in zs if not np.isnan(z)]
    complementary = bool(
        finite and len(finite) == len(seps) and worst_sep > 0
        and len(finite_z) == len(zs) and max(finite_z) < z_threshold
    )
    return ComplementaryExposure(
        separation=round(worst_sep, 4) if not np.isnan(worst_sep) else float("nan"),
        z_a=round(min(finite_z), 4) if finite_z else float("nan"),
        z_b=round(max(finite_z), 4) if finite_z else float("nan"),
        is_complementary=complementary,
    )

from __future__ import annotations
from collections import defaultdict
import random
from itertools import product
from typing import List, Dict, Iterable, Any, Optional


def kmer_composition(text: str, k: int) -> Iterable[str]:
    """Forms the k-mer composition of a string."""
    return [text[i:i+k] for i in range(len(text)-k+1)]

def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    return path[0] + ''.join([p[-1] for p in path[1:]])

def overlap_graph(patterns: List[str]) -> Dict[str, list[str]]:
    """Forms the overlap graph of a collection of patterns."""
    prefix = lambda s: s[:-1]
    suffix = lambda s: s[1:]
    graph = defaultdict(list)
    for k1 in patterns:
        for k2 in patterns:
            if suffix(k1) == prefix(k2):
                if k2 not in graph[k1]:
                    graph[k1].append(k2)
    return dict(graph)


class EdgeError(Exception):
    """Container for edge-related exceptions."""


class Node:
    """Represents a k-mer."""
    def __init__(self, val: str) -> None:
        self.val = val

    def __str__(self) -> str:
        return f"{self.val}"

    def __repr__(self):
        return f"{self.val}"

    def __len__(self) -> int:
        return len(self.val)

    def __hash__(self) -> int:
        return hash(self.val)

    def __eq__(self, other: Node) -> bool:
       return self.val == other.val

    @property
    def suffix(self) -> str:
        return self.val[1:]

    @property
    def prefix(self) -> str:
        return self.val[:-1]


class Edge:
    def __init__(self, n1: Node, n2: Node) -> None:
        assert n1.suffix == n2.prefix, EdgeError(
            "Source node suffix must overlap with destination node prefix"
        )
        self.src = n1
        self.dest = n2
        self.label = n1.val + n2.val[-1]

    def __str__(self) -> str:
        return f"{self.src}->{self.dest}"

    def __repr__(self) -> str:
        return f"{self.src}->{self.dest}"

    @property
    def prefix(self) -> str:
        return self.src.val

    @property
    def suffix(self) -> str:
        return self.dest.val


class BaseGraph:
    def __init__(self, text: str, k: int) -> None:
        self.text = text
        self.k = k
        self.nodes = [Node(kmer) for kmer in kmer_composition(text, k)]
        self.edges = []
        self.G = self._generate_graph_from_text()

    def __str__(self) -> str:
        fmt = ''
        for node, node_list in self.G.items():
            fmt += f"{node}: {' '.join([str(n) for n in node_list])}\n"
        return fmt

    def _generate_graph_from_text(self) -> dict[Node, List[Node]]:
        """Create an overlap graph represented as an adjacenct list."""
        G = defaultdict(list)
        for n1 in self.nodes:
            for n2 in self.nodes:
                try:
                    e = Edge(n1, n2)
                    if n2 not in G[n1]:
                        G[n1].append(n2)
                        self.edges.append(e)
                except (AssertionError, EdgeError) as err:
                    # print(f"{n1} does not overlap with {n2}")
                    continue
        return dict(G)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def unique_nodes(self) -> set[Node]:
        return set(self.nodes)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    def to_adjacency_list(self) -> Dict[str, List[str]]:
        """Convert the graph into an adjacency list."""
        adj = defaultdict(list)
        for src, dest_list in self.G.items():
            for dest in dest_list:
                adj[src.val].append(dest.val)
        return dict(adj)


class DeBruijn(BaseGraph):
    """A DeBruijn graph representation."""
    def __init__(self, text: Optional[str] = None, k: Optional[int] = 3, kmers: Optional[List[str]] = None) -> None:
        # Intialize graph from a string `text` and int `k`
        if text:
            assert k is not None, "k must be given as an integer"
            super().__init__(text, k)
            self.nodes = [Node(text[i:i+k-1]) for i in range(len(text)-k+2)]
            self.edges = [Edge(self.nodes[i], self.nodes[i+1]) for i in range(len(self.nodes)-1)]
            self.G = self._generate_graph_from_text()
        # Initialize graph from a list of k-mers
        elif kmers:
            substrings = [kmer[:-1] for kmer in kmers] + [kmer[1:] for kmer in kmers]
            self.unique_kmers = set(substrings)
            self.nodes = [Node(kmer) for kmer in self.unique_kmers]
            self.edges = [Edge(Node(kmer[:-1]), Node(kmer[1:])) for kmer in kmers]
            self.G = self._generate_graph_from_reads()
        else:
            raise ValueError("Provide either a single string or a list of string of generating the graph")

    def _generate_graph_from_text(self) -> dict[Node, List[Node]]:
        """Create a DeBruijn graph from a string using nodes of length k."""
        G = defaultdict(list)
        for i in range(len(self.text)-self.k+1):
            src = Node(self.text[i:i+self.k-1])
            dest = Node(self.text[i+1:i+self.k])
            G[src].append(dest)
        return dict(G)

    def _generate_graph_from_reads(self) -> dict[Node, List[Node]]:
        """Create a DeBruijn graph from a set of reads."""
        G = defaultdict(list)
        for edge in self.edges:
            for kmer in self.unique_kmers:
                if edge.prefix == kmer:
                    src = Node(kmer)
                    dest = Node(edge.suffix)
                    G[src].append(dest)
        return dict(G)


def traverse(start: str, graph: DeBruijn) -> list[str]:
    cycle = []
    G = graph.to_adjacency_list()
    curr_pos = random.choice(list(G.keys()))
    while curr_pos is not None:
        next_pos_choices = G.get(curr_pos, None)
        next_pos = random.choice(next_pos_choices)
        print(f"{curr_pos} -> {next_pos}")
        cycle.append(curr_pos)
        curr_pos = next_pos
    print("Dead end. Cycle has terminated.")
    return cycle

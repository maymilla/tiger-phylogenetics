import random
from collections import defaultdict

import numpy as np

from kimura_distance import build_distance_matrix
from neighbor_joining import neighbor_joining

def resample_alignment(sequences_dict):
    names = list(sequences_dict.keys())
    length = len(sequences_dict[names[0]])
    sampled_indices = random.choices(range(length), k=length)

    resampled = {}
    for name in names:
        seq = sequences_dict[name]
        resampled[name] = ''.join(seq[i] for i in sampled_indices)
    
    return resampled

def run_bootstrap(sequences_dict, n_bootstrap=1000):
    clade_counts = defaultdict(int)
    success = 0

    for iteration in range(n_bootstrap):
        resampled = resample_alignment(sequences_dict)
        matrix, species = build_distance_matrix(resampled)
        
        if np.any(np.isnan(matrix)):
            continue
        
        edges = neighbor_joining(matrix, species)
        clades = get_clades_from_edges(edges)
        for clade in clades:
            clade_counts[clade] += 1
        
        success += 1
        if (iteration + 1) % 100 == 0:
            print(f"Bootstrap progress: {iteration+1}/{n_bootstrap}")
    
    # convert count to support value
    if success == 0:
        print("Warning: semua iterasi bootstrap gagal")
        return {}

    support = {clade: (count / success * 100) for clade, count in clade_counts.items()}
    print(f"Bootstrap selesai: {success}/{n_bootstrap} iterasi berhasil")
    
    return support 

def get_clades_from_edges(edges):
    # helper function to extract all clades from NJ edges list
    children = defaultdict(list)
    all_children = set()
    all_parents = set()

    for parent, child, _ in edges:
        children[parent].append(child)
        all_children.add(child)
        all_parents.add(parent)

    leaf_nodes = all_children-all_parents

    def get_leaves_under(node):
        if node in leaf_nodes:
            return frozenset([node])
        result = frozenset()
        for child in children[node]:
            result = result | get_leaves_under(child)
        return result
    
    clades = []
    for node in all_parents:
        clade = get_leaves_under(node)
        if len(clade) > 1:
            clades.append(clade)
    
    return clades

def attach_bootstrap_to_tree(edges, bootstrap_support):
    children = defaultdict(list)
    all_children = set()
    all_parents = set()

    for parent, child, _ in edges:
        children[parent].append(child)
        all_children.add(child)
        all_parents.add(parent)

    leaf_nodes = all_children-all_parents

    def get_leaves_under(node):
        if node in leaf_nodes:
            return frozenset([node])
        result = frozenset()
        for child in children[node]:
            result = result | get_leaves_under(child)
        return result
    
    node_support = {}
    for node in all_parents:
        clade = get_leaves_under(node)
        support_value = bootstrap_support.get(clade, 0.0)
        node_support[node] = round(support_value, 1)
    
    return node_support

def main(sequences_dict, edges, n_bootstrap=1000):
    bootstrap_support = run_bootstrap(sequences_dict, n_bootstrap)
    node_support = attach_bootstrap_to_tree(edges, bootstrap_support)
    
    print("\nBootstrap support per node:")
    for node, support in node_support.items():
        status = "kuat" if support >= 70 else "lemah"
        print(f"  {node}: {support:.1f}% ({status})")
    
    return node_support
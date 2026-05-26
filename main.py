from fetch_sequence import main as fetch
from smith_waterman import main as align
from run_mafft import main as mafft
from kimura_distance import main as kimura
from neighbor_joining import main as nj
from bootstrap_sampling import main as bootstrap
from evolutionary_gap import main as gap
from visualizer import main as visualize

if __name__ == "__main__":
    print("=== Step 1: Fetch Sequences ===")
    fetch()

    print("\n=== Step 2: Pairwise Alignment ===")
    align()

    print("\n=== Step 3: MAFFT ===")
    mafft()

    print("\n=== Step 4: Kimura 2-Parameter Distance ===")
    sequences_dict, matrix, species = kimura()

    print("\n=== Step 5: Neighbor-Joining Tree ===")
    edges = nj(matrix, species)

    print("\n=== Step 6: Bootstrap ===") # 1000x sampling
    node_support = bootstrap(sequences_dict, edges)

    print("\n=== Step 7: Evolutionary Gap Analysis ===")
    loss, distinctiveness = gap(edges, matrix, species, node_support)

    print("\n=== Step 8: Visualisasi Pohon ===")
    visualize(edges, node_support)
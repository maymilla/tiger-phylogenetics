import math
import numpy as np

def parse_fasta(filepath):
    sequences = {}
    current_name = None

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue
            
            if line.startswith(">"):
                # header
                current_name = line[1:].strip()
                sequences[current_name] = ""
            else:
                # sequences
                sequences[current_name] += line.lower()
            
    return sequences

def is_valid_site(base1, base2):
    valid = {'a', 'c', 'g', 't'}
    return base1 in valid and base2 in valid

def classify_substitution(base1, base2):
    # classify transition if two bases come from the same group and transversion if otherwise
    if base1 == base2:
        return None

    purin = {'a', 'g'}
    pirimidin = {'c', 't'}

    if (base1 in purin and base2 in purin) or (base1 in pirimidin and base2 in pirimidin):
        return "transition"
    else:
        return "transversion"

def compute_k2p(seq1, seq2):
    # calculate kimura 2 param
    # -0.5 * ln(1 - 2P - Q) - 0.25 * ln(1 - 2Q)
    n_valid = 0
    n_transition = 0
    n_transversion = 0

    for base1, base2 in zip(seq1, seq2):
        if not is_valid_site(base1, base2):
            continue
        n_valid += 1
        if base1 != base2:
            res = classify_substitution(base1, base2)
            if res == "transition":
                n_transition += 1
            elif res == "transversion":
                n_transversion += 1
            
    if n_valid < 100: # threshold = 100, can be adjusted later
        return None

    P = n_transition / n_valid
    Q = n_transversion / n_valid

    if (1-2*P-Q <= 0) or (1-2*Q <= 0):
        return None
    
    d = -0.5*math.log(1-2*P-Q)-0.25*math.log(1-2*Q)
    return d, n_valid, n_transition, n_transversion

def build_distance_matrix(sequences_dict):
    species = list(sequences_dict.keys())
    N = len(species)
    matrix = np.zeros((N, N))

    for i in range(N):
        for j in range(i, N):
            if i == j:
                matrix[i][j] = 0.0
            else:
                res = compute_k2p(sequences_dict[species[i]], sequences_dict[species[j]])
                if res is None:
                    matrix[i][j] = np.nan
                    matrix[j][i] = np.nan
                    print(f"Warning: {species[i]} vs {species[j]} tidak bisa dihitung")
                else:
                    matrix[i][j] = res[0] # d value is stored as first element (idx = 0)
                    matrix[j][i] = res[0]

    return matrix, species

def main():
    sequences_dict = parse_fasta("data/aligned_sequences.fasta")
    matrix, species = build_distance_matrix(sequences_dict)
    
    short_names = [sp.split("|")[0].strip().replace("P_t_", "") for sp in species]
    print("\nDistance Matrix (K2P):")
    print(f"{'':15}", end="")
    for name in short_names:
        print(f"{name:>12}", end="")
    print()
    for i, name in enumerate(short_names):
        print(f"{name:15}", end="")
        for j in range(len(species)):
            val = matrix[i][j]
            print(f"{val:>12.4f}" if not np.isnan(val) else f"{'NaN':>12}", end="")
        print()
    
    return sequences_dict, matrix, species
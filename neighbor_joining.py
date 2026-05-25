import numpy as np

def compute_q_matrix(matrix, species_remaining):
    # q matrix is correction distance matrix for NJ
    N = len(species_remaining)
    q_matrix = np.zeros((N, N))
    row_sums = np.sum(matrix, axis=1)

    for i in range(N):
        for j in range(N):
            if i == j:
                q_matrix[i][j] = np.inf
                continue
            q_matrix[i][j] = (N-2) * matrix[i][j] - row_sums[i] - row_sums[j]

    return q_matrix

def find_closest_pair(q_matrix):
    MIN = np.inf
    for i in range(len(q_matrix)):
        for j in range(len(q_matrix)):
            if q_matrix[i][j] < MIN:
                MIN = q_matrix[i][j]
                min_i = i
                min_j = j
    
    return min_i, min_j
            

def compute_branch_lengths(matrix, i, j, species_remaining):
    # calculate branch length from new node to species i and j
    N = len(species_remaining)
    row_sums = np.sum(matrix, axis=1)
    
    branch_i = 0.5*matrix[i][j]+(row_sums[i]-row_sums[j])/(2*(N-2))
    branch_j = matrix[i][j]-branch_i

    return branch_i, branch_j

def update_matrix(matrix, i, j, species_remaining):
    # after i and j joined as node u, update matrix
    N = len(species_remaining)
    new_distances = []
    for k in range(N):
        if k == i or k == j:
            continue
        d_uk = 0.5*(matrix[i][k]+matrix[j][k]-matrix[i][j])
        new_distances.append(d_uk)
    
    new_matrix = np.delete(matrix, [i, j], axis=0)
    new_matrix = np.delete(new_matrix, [i, j], axis=1)

    new_row = np.array(new_distances)
    new_matrix = np.vstack([new_matrix, new_row])
    
    new_col = np.array(new_distances+[0.0]).reshape(-1, 1)
    new_matrix = np.hstack([new_matrix, new_col])

    node_name = f"Node_{N}"
    new_species = [s for idx, s in enumerate(species_remaining) if idx != i and idx != j]
    new_species.append(node_name)
    
    return new_matrix, new_species

def neighbor_joining(matrix, species):
    current_matrix = matrix.copy()
    current_species = species.copy()
    edges = []
    node_counter = 0

    while len(current_species) > 2:
        Q = compute_q_matrix(current_matrix, current_species)
        i, j = find_closest_pair(Q)
        branch_i, branch_j = compute_branch_lengths(current_matrix, i, j, current_species)
        
        node_name = f"Node_{node_counter}"
        node_counter += 1
        
        edges.append((node_name, current_species[i], branch_i))
        edges.append((node_name, current_species[j], branch_j)) 

        current_matrix, current_species = update_matrix(current_matrix, i, j, current_species)

    last_branch = current_matrix[0][1] / 2
    edges.append((current_species[0], current_species[1], last_branch))

    return edges
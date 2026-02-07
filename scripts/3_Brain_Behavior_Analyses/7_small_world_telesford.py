import numpy as np
import pandas as pd
import os

from bct.algorithms.distance import distance_wei
from bct.algorithms.distance import charpath
from bct.algorithms.clustering import clustering_coef_wu

# --- Set functions
def neg_log_transform(A):
    """
    Performs negative log transformation to generate a length matrix
    for characteristic path length calculation.
    Inputs:
        A (ndarray): Symmetric (n x n) adjacency matrix.
    Output:
        M (ndarray): Log transformed length matrix.
    """
    M = np.copy(A)
    nz_idx = np.nonzero(M)
    M[nz_idx] = -np.log10(M[nz_idx])
    M[M==0] = np.inf
    np.fill_diagonal(M,0)
    return M

def characteristic_path_length(A):
    """
    Calculates charactertic path length from an adjacency matrix where
    edge weights describe information carrying capacity.
    Inputs:
        A (ndarray): Symmetric (n x n) adjacency matrix.
    Output:
        cpl (scalar): Characteristic path length.
    """
    M = np.copy(A)
    L = neg_log_transform(M)
    D, _ = distance_wei(L)
    cpl, _, _, _, _ = charpath(D)
    return(cpl)

def regular_matrix_generator(G, r):
    """
    Generates a regular matrix with weights and # of nodes obtained from
    a given adjacency matrix.
    Inputs:
        G (ndarray): Symmetric (n x n) adjacency matrix.
        r (int): Approximate radius of the regular network
    Output:
        M (ndarray): Regular matrix with highest weights in inner radius.
    """
    n = G.shape[0]
    G_upper = np.triu(G)
    B = G_upper.flatten()
    B = np.sort(B)[::-1]
    num_els = int(np.ceil(G.size / (2 * n)))
    num_zeros = 2 * n * num_els - G.size
    B = np.concatenate([B, np.zeros(num_zeros)])
    B = B.reshape((n, -1), order='F')  # Fortran order to mimic column-wise filling
    M = np.zeros_like(G)
    for i in range(n):
        for z in range(r):
            a = np.random.randint(0, n)
            while (B[a, z] == 0 and z != r - 1) or (B[a, z] == 0 and z == r - 1 and np.any(B[:, r - 1])):
                a = np.random.randint(0, n)
            y_coor_1 = (i + z) % n
            M[i, y_coor_1] = B[a, z]
            M[y_coor_1, i] = B[a, z]
            B[a, z] = 0
    return M

def random_rewire_matrix(G):
    """
    Randomly rewires all connections of a given adjacency matrix.
    Inputs:
        G (ndarray): Symmetric (n x n) adjacency matrix.
    Ouput:
        G_rewired (ndarray): Symmetrically rewired matrix with the same weight distribution.
    """
    n = G.shape[0]
    triu_indices = np.triu_indices(n, k=1)
    upper_values = G[triu_indices]
    np.random.shuffle(upper_values)
    G_rewired = np.zeros((n, n))
    G_rewired[triu_indices] = upper_values
    G_rewired = G_rewired + G_rewired.T
    return G_rewired

def small_world_telesford(C_real, L_real, C_latt, L_rand):
    """
    Calculates the Telesford small-world metric from observed topology metrics, and
    those generated from null models.
    Inputs:
        C_real (scalar): Observed average clustering coefficient.
        L_real (scalar): Observed characteristic path length.
        C_latt (ndarray): Averaged clustering coefficient of null lattice network.
        L_rand (ndarray): Characterstic path length of null random network.
    Output:
        Sigma: (ndarray): Telesford small-world metric
        dict: (dict): Dictionary of the average clustering coefficient ratio, and characteristic path length ratio.
    """
    C_ratio = C_real / C_latt
    L_ratio = L_rand / L_real
    sigma = L_ratio - C_ratio
    return sigma, {
        "C_ratio": C_ratio,
        "L_ratio": L_ratio
    }

# --- Set working directory
rdir=os.getcwd()

# --- Get subject list
subjects=np.loadtxt(os.path.join(rdir,"data","mri_qaqc","good_subjects.csv"), skiprows=1)
subjects=subjects.tolist()

# --- Generate small-worldness values and save
results=[]
for sub in subjects:
    mat=np.genfromtxt(os.path.join(rdir,"data","jsfm",sub,"MMP1_maxed_averaged_structure_function_network.csv"),delimiter=" ")
    mat_bin = (mat > 0).astype(int)
    d = np.sum(mat_bin, axis=1)
    r = np.ceil(np.mean(d) / 4).astype(int)
    regular_matrix = regular_matrix_generator(mat, r=r)
    random_matrix = random_rewire_matrix(mat)
    C_real = np.mean(clustering_coef_wu(mat))
    C_latt = np.mean(clustering_coef_wu(regular_matrix))
    L_real = characteristic_path_length(mat)
    L_rand = characteristic_path_length(random_matrix)
    sigma, metrics = small_world_telesford(C_real, L_real, C_latt, L_rand)
    row = {
        'Subject': sub,
        'Sigma': sigma,
        'C_real': C_real,
        'L_real': L_real,
        'C_latt': C_latt,
        'L_rand': L_rand,
        'C_ratio': metrics['C_ratio'],
        'L_ratio': metrics['L_ratio']
    }
    results.append(row)

df=pd.DataFrame(results)
df.to_csv(os.path.join(rdir,"data","graph_metrics","master_small_worldness_telesford.csv"),index=False)

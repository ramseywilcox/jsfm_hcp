import os

import numpy as np
import pandas as pd

from scipy.linalg import schur

from nctpy.utils import matrix_normalization

# --- Set function
def modal_control_continuous(A_norm):
    T, U = schur(A_norm, 'real')
    eigVals = np.diag(T)
    N = A_norm.shape[0]
    phi = np.zeros(N, dtype=float)
    for i in range(N):
        Al = U[i, ] * U[i, ]
        Ar = (1.0 - np.exp(eigVals)).transpose()
        phi[i] = np.matmul(Al, Ar)
    return phi


# --- Set working directory
rdir=os.getcwd()

# --- Get subject list
subjects=np.loadtxt(os.path.join(rdir,"data","mri_qaqc","good_subjects.csv"), skiprows=1)
subjects=subjects.tolist()

# --- Generate modal control values and save
results=[]
for sub in subjects:
    mat=np.genfromtxt(os.path.join(rdir,"data","tractography",sub,"MMP1_structural_connectome_mu_scaled"),delimiter=" ")
    mat_norm=matrix_normalization(mat, system="continuous",c=1)
    phi = modal_control_continuous(mat_norm)
    row = {
        'Subject': sub,
        'Phi': phi,
    }
    results.append(row)
else:
    print("Subject has no data")

rows = []
for subject in results:
    row = {f'Node_{i+1}': val for i, val in enumerate(subject['Phi'])}
    row['Subject'] = subject['Subject']
    rows.append(row)

df = pd.DataFrame(rows)
df = df[['Subject'] + [f'Node_{i+1}' for i in range(360)]]
df.to_csv(os.path.join(rdir,"graph_metrics","modal_control_continuous_structural_connectome_mu_scaled.csv"),index=False)

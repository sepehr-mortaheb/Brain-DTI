#%% Importing 
import os 
import os.path as op 
import pandas as pd 
import matplotlib.pylab as plt 
import seaborn as sns 
import numpy as np 
import pygsp
import nibabel as nb
import nilearn as nl
from scipy.io import savemat, loadmat

# %% Parameters and Directories Initialization

data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

atlas = 'SCH100'
subjects = np.sort(os.listdir(data_dir))
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}


# %% Reading and Averaging SC matrices
TSC=[]
for sub in subjects:
    for ses in sessions[sub]:
        sc = np.loadtxt(
            op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
            delimiter=','
        )
        TSC.append(sc)
TSC = np.array(TSC)
SC = np.mean(
    TSC, 
    axis=0
)

fig, ax = plt.subplots(
    1,
    1, 
    figsize=(9,7)
)
sns.heatmap(
    SC, 
    ax=ax, 
    square=False, 
    cmap='viridis',
    xticklabels=[],
    yticklabels=[], 
)
fig.savefig(
    op.join(
        res_dir, 
        f'SC/SC_avg_{atlas}.png'
    ), 
    dpi=300
)
fig.savefig(
    op.join(
        res_dir, 
        f'SC/SC_avg_{atlas}.pdf'
    ), 
    dpi=300
)

savemat(op.join(res_dir, f'SC/SC_avg_{atlas}.mat'), {'SC':SC})

# %% Graph Creation
G = pygsp.graphs.Graph(SC)
G.compute_laplacian('normalized') 
G.compute_fourier_basis(recompute=True)
eigenValues = G.e
eigenVectors = G.U

fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.plot(
    eigenValues, 
    marker='o', 
    linewidth=0
)
ax.grid(True)
ax.axhline(
    0, 
    color='black'
)
ax.axvline(
    0, 
    color='black'
)
ax.set_ylabel(
    'SC Eigen Values', 
    size='15'
)

fig.savefig(
    op.join(
        res_dir, 
        f'SC/SC_eig_{atlas}.png'
    ), 
    dpi=300
)
fig.savefig(
    op.join(
        res_dir, 
        f'SC/SC_eig_{atlas}.pdf'
    ), 
    dpi=300
)

savemat(
    op.join(
        res_dir, 
        f'SC/SC_eig_{atlas}.mat'
    ), 
    {
        'e':eigenValues, 
        'U':eigenVectors
    }
)

#%% Zero Crossing Rate 
R = len(eigenValues)
wzc = np.zeros(R)
for r in range(R):
    u = eigenVectors[:, r];
    summ = 0; 
    for i in range(R-1):
        for j in range(i+1, R):
            if u[i]*u[j] < 0:
                summ += SC[i, j]>1
            wzc[r] = summ

fig, ax = plt.subplots(
    1,
    1, 
    figsize=(10,5)
)
ax.plot(
    eigenValues,
    wzc,  
    linewidth=2
)
ax.grid(True)
ax.axhline(
    0, 
    color='black'
)
ax.axvline(
    0, 
    color='black'
)
ax.set_ylabel(
    'Zero Crossing Rate', 
    size='15'
)
ax.set_xlabel(
    'Eigenvalue', 
    size='15'
)
fig.savefig(
    op.join(
        res_dir, 
        f'SC/SC_ZCR_{atlas}.png'
    ), 
    dpi=300
)
fig.savefig(
    op.join(
        res_dir, 
        f'SC/SC_ZCR_{atlas}.pdf'
    ), 
    dpi=300
)


#%% Creating Nifti Files of Eigenvectors 
atlas_addr = op.join(mni_dir, f'{atlas}.nii')
   

# Eigenvectors of Interest
evoi = [0, 1, 2, R-3, R-2, R-1]
for i in evoi:
    atlas_img = nb.load(atlas_addr)
    atlas_hdr = atlas_img.header
    atlas_aff = atlas_img.affine 
    atlas_mat = atlas_img.get_fdata()
    for r in range(R):
        atlas_mat[atlas_mat==r+1] = eigenVectors[r, i]
    new_img = nb.Nifti1Image(atlas_mat, atlas_aff, atlas_hdr)
    nb.save(new_img, op.join(res_dir, f'SC/eig_{i}_{atlas}.nii'))
    

np.unique(atlas_mat)
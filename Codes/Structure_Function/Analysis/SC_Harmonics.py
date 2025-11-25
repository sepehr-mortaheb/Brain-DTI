#%% Importing

import os 
import os.path as op 
import pandas as pd 
import matplotlib.pylab as plt 
import seaborn as sns 
import numpy as np 
import pygsp
from scipy.io import savemat, loadmat
from nilearn import plotting, datasets, image
from scipy.ndimage import center_of_mass

#%% Parameters and directories initialization 

data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

atlas = 'SCH100' # SCH100, SCH400, AAL

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

#%% Calculate the Average SC for Cosmonauts 

exc_list = [
  "sub-cosmonaut13", 
  "sub-cosmonaut15", 
  "sub-cosmonaut16",
  "sub-cosmonaut26"
] 

TSC=[]
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        if sub not in exc_list:
            for ses in sessions[sub]:
                flight = ses.split('-')[1][0:2]
                tp = ses.split('-')[1][2:]
                if (flight=='f1') & (tp in ['pre2', 'post']):
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

savemat(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}.mat'),
    {'SC':SC}
)

fig, ax = plt.subplots(
    1,
    1, 
    figsize=(9,7)
)
epsilon = 1e-5
SC_log = np.log10(SC+epsilon)

sns.heatmap(
    SC, 
    ax=ax, 
    square=False, 
    cmap='viridis',
    xticklabels=[],
    yticklabels=[], 
)

fig.savefig(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}.pdf'),
    dpi=300
)

fig.savefig(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}.png'),
    dpi=300
)

#%% Normalized Laplacian Eigendecomposition 

G = pygsp.graphs.Graph(SC)
G.compute_laplacian('normalized') 
G.compute_fourier_basis(recompute=True)
eigenValues = G.e
eigenVectors = G.U

savemat(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_eigenvalues.mat'),
    {'eigVals': eigenValues}
)

savemat(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_eigenvectors.mat'),
    {'eigVecs': eigenVectors}
)

#%% Plotting Eigne Values 

fig, ax = plt.subplots(1, 1, figsize=(7,4))
ax.plot(
    eigenValues, 
    'o',
    markersize=4
)
ax.grid(True)
ax.axvline(0, color='black', linewidth=1)
ax.axhline(0, color='black', linewidth=1)
ax.set_ylabel('Normalized Laplacian Eigenvalues', size=13)
fig.savefig(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_eigenvalues.pdf'),
    dpi=300
)
fig.savefig(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_eigenvalues.png'),
    dpi=300
)

#%% Zero Crossing Rate of Eigenvectors 

R = len(eigenValues)
wzc = np.zeros(R)

for r in range(R):
    u = eigenVectors[:, r]
    summ = 0; 
    for i in range(R-1):
        for j in range(i+1, R):
            if u[i]*u[j] < 0:
                summ += SC[i, j]>1
            wzc[r] = summ

fig, ax = plt.subplots(
    1,
    1, 
    figsize=(7,4)
)
ax.plot(
    eigenValues,
    wzc,  
    linewidth=2
)
ax.grid(True)
ax.axhline(
    0, 
    color='black', 
    linewidth=1
)
ax.axvline(
    0, 
    color='black',
    linewidth=1
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
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_zcr.pdf'), 
    dpi=300
)
fig.savefig(
    op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_zcr.png'), 
    dpi=300
)

#%% Plotting Eigenvectors on the brain 

if atlas == 'SCH100':
    atlas_file = datasets.fetch_atlas_schaefer_2018(
        n_rois=100, 
        yeo_networks=7, 
        resolution_mm=1
    )

atlas_img = image.load_img(atlas_file.maps)
labels = atlas_file.labels  # Labels from "1" to "R"
region_ids = [f"R{i}" for i in range(1, R+1)]
atlas_data = image.load_img(atlas_img).get_fdata()
region_coords = []
for i, region_label in enumerate(region_ids, 1):
    mask = atlas_data == i
    if np.sum(mask) == 0:
        continue
    com = center_of_mass(mask)
    affine = atlas_img.affine
    mni_coords = np.dot(affine[:3, :3], com) + affine[:3, 3]
    region_coords = region_coords + [mni_coords.tolist()]

vec = eigenVectors[:, 1]
th = np.max(np.abs(vec))
plotting.plot_markers(
    vec, 
    region_coords, 
    node_size=40,
    node_vmin=-th,
    node_vmax=th,
    node_cmap='RdBu_r'
)

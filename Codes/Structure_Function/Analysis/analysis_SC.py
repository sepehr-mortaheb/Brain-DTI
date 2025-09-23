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
import networkx as nx
from scipy.io import savemat, loadmat
from nilearn import plotting, datasets, image
from nilearn.maskers import NiftiLabelsMasker
from scipy.ndimage import center_of_mass

#%% Parameters and Directories Initialization

data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

atlas = 'SCH100'
R = 100 
subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

#%% Reading SC files for Cosmonauts and Creating the dataframe

df = pd.DataFrame([])
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        print(sub)
        for ses in sessions[sub]:
            print(f'---- {ses}')
            flight = ses.split('-')[1][0:2]
            tp = ses.split('-')[1][2:]

            sc = np.loadtxt(
                op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                delimiter=','
            )

            for row in range(R):
                for col in range(row+1,R):
                    edge = f'R{row}-R{col}'
                    val = sc[row, col]
                    dftmp = pd.DataFrame([])
                    dftmp['subject'] = [sub]
                    dftmp['flight'] = [flight]
                    dftmp['time'] = [tp]
                    dftmp['edge'] = [edge]
                    dftmp['val'] = [val]
                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'SC/SC_vals_{atlas}.xlsx'))

#%% Post vs Pre analysis plots
# The LMM analysis has been done in R and the results have been saved in an 
# Excel file that can be loaded and be visualized. 

# Reading the atlas information 
atlas = datasets.fetch_atlas_schaefer_2018(
    n_rois=100, 
    yeo_networks=7, 
    resolution_mm=1
)
atlas_img = image.load_img(atlas.maps)
labels = atlas.labels  # Labels from "1" to "100"
region_ids = [f"R{i}" for i in range(1, 101)]
atlas_data = image.load_img(atlas_img).get_fdata()
region_coords = {}
for i, region_label in enumerate(region_ids, 1):
    mask = atlas_data == i
    if np.sum(mask) == 0:
        continue
    com = center_of_mass(mask)
    affine = atlas_img.affine
    mni_coords = np.dot(affine[:3, :3], com) + affine[:3, 3]
    region_coords[region_label] = mni_coords.tolist()

# Reading the stat results 
file_addr = op.join(
    res_dir,
    'SC/MassUnivariate/SC_MassUnivariate_LMM_SCH100.xlsx'
)
df = pd.read_excel(file_addr)

# Filter significant edges
sig_edges = df[df["p_fdr"] < 0.05].copy()
sig_edges[["node1", "node2"]] = sig_edges["edge"].str.split("-", expand=True)
sig_edges = sig_edges[
    sig_edges["node1"].isin(region_coords) & sig_edges["node2"].isin(region_coords)
]

# Graph Creation
G = nx.Graph()
for _, row in sig_edges.iterrows():
    n1, n2 = row["node1"], row["node2"]
    weight = row["estimate"]
    G.add_edge(n1, n2, weight=weight)

# Plotting the Graph on the Brain Atlas 
coords = [region_coords[n] for n in G.nodes()]

f = plotting.plot_connectome(
    adjacency_matrix=nx.to_numpy_array(G, weight="weight"),
    node_coords=coords,
    node_color='black',
    node_size=5,
    edge_threshold="1%",
    edge_cmap="bwr",
    colorbar=True,
    display_mode="ortho",
    annotate=False
)

plotting.show()
f.savefig(
    op.join(
        res_dir, f'SC/MassUnivariate/SC_Post-Pre_{atlas}.png'
    ),
    dpi=300,
    bbox_inches='tight'
)
#%% Reading and Averaging SC matrices
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
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

subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

#%% Reading Cosmonauts' and Controls' SC files and dataframe creation

atlas = 'SCH100' # SCH100, SCH400, AAL

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

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

df.to_csv(op.join(res_dir, f'SC/SC_Cosmonauts_{atlas}.csv'))

df = pd.DataFrame([])
for sub in subjects:
    if sub.startswith('sub-control'):
        print(sub)
        for ses in sessions[sub]:
            print(f'---- {ses}')
            tp = ses.split('-')[1]

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
                    dftmp['time'] = [tp]
                    dftmp['edge'] = [edge]
                    dftmp['val'] = [val]
                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_csv(op.join(res_dir, f'SC/SC_Controls_{atlas}.csv'))


#%% Post vs Pre analysis plots
# The LMM analysis has been done in R and the results have been saved in an 
# Excel file that can be loaded and be visualized. 

# Reading the atlas information 

if atlas == 'SCH100':
    atlas_file = datasets.fetch_atlas_schaefer_2018(
        n_rois=100, 
        yeo_networks=7, 
        resolution_mm=1
    )
elif atlas == 'SCH400':
    atlas_file = datasets.fetch_atlas_schaefer_2018(
        n_rois=400, 
        yeo_networks=7, 
        resolution_mm=1
    )
elif atlas == 'AAL':
    atlas_file = datasets.fetch_atlas_aal()

atlas_img = image.load_img(atlas_file.maps)
labels = atlas_file.labels  # Labels from "1" to "R"
region_ids = [f"R{i}" for i in range(1, R+1)]
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
    f'SC/MassUnivariate/SC_MassUnivariate_LMM_{atlas}.xlsx'
)
df = pd.read_excel(file_addr)

# Filter significant edges
sig_edges = df[df["p_fdr"] < 0.05].copy()

if len(sig_edges) == 0:
    print('No significant edge was found!')
    df[["node1", "node2"]] = df["edge"].str.split("-", expand=True)
    df = df[df["node1"].isin(region_coords) & df["node2"].isin(region_coords)]
    # Graph Creation
    G = nx.Graph()
    for _, row in df.iterrows():
        n1, n2 = row["node1"], row["node2"]
        weight = 0
        G.add_edge(n1, n2, weight=weight)
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in G.nodes()]
    f = plotting.plot_connectome(
        adjacency_matrix=nx.to_numpy_array(G, weight="weight"),
        node_coords=coords,
        node_color='black',
        node_size=5,
        colorbar=False,
        display_mode="ortho",
        annotate=False
    )   

else:
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
    dpi=300
)
f.savefig(
    op.join(
        res_dir, f'SC/MassUnivariate/SC_Post-Pre_{atlas}.pdf'
    ),
    dpi=300
)

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
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import mne
from mne.viz import circular_layout
from mne_connectivity.viz import plot_connectivity_circle

#%% Parameters and Directories Initialization

data_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

#%% Plotting the results of Type 1 Statistical Analysis: Cosmonauts (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'

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
labels = atlas_file.labels[1:]  # Labels from "1" to "R"

R = len(labels)
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
    f'{atlas}/FC/MassUnivariate/Statistics_Cosmonauts_LMM_PostPre.xlsx'
)

df = pd.read_excel(file_addr)

# Filter significant edges
sig_edges = df[df["p_fdr"] < 0.05].copy()

# Creating node colors 
if atlas=='SCH100':
    roi_to_network = {i: net for i, net in enumerate(
        ['Visual'] * 9 + ['Somatomotor'] * 6 + ['DorsalAttention'] * 8 +
        ['VentralAttention'] * 7 + ['Limbic'] * 3 + ['Control'] * 4 + ['DMN'] * 13 + 
        ['Visual'] * 8 + ['Somatomotor'] * 8 + ['DorsalAttention'] * 7 +
        ['VentralAttention'] * 5 + ['Limbic'] * 2 + ['Control'] * 9 + ['DMN'] * 11
    )}
    network_colors = {
        'Visual': 'brown',
        'Somatomotor': 'green',
        'DorsalAttention': 'orange',
        'VentralAttention': 'olive',
        'Limbic': 'black',
        'Control': 'cyan',
        'DMN': 'magenta'
    }
    node_colors = [network_colors[roi_to_network[i]] for i in range(len(roi_to_network))]

if len(sig_edges) == 0:
    print('No significant edge was found!')
    df[["node1", "node2"]] = df["edge"].str.split("-", expand=True)
    df = df[df["node1"].isin(region_coords) & df["node2"].isin(region_coords)]
    # Graph Creation
    A = np.zeros((100, 100))
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids]
    view = plotting.view_connectome(
        A, 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser()    

else:
    sig_edges[["node1", "node2"]] = sig_edges["edge"].str.split("-", expand=True)
    sig_edges = sig_edges[
        sig_edges["node1"].isin(region_coords) & sig_edges["node2"].isin(region_coords)
    ]
    # Graph Creation
    A = np.zeros((100, 100))
    for _, row in sig_edges.iterrows():
        n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
        weight = row["beta"]
        A[n1, n2] = weight
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids] 
    view = plotting.view_connectome(
        A, 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser()

#%% Plot Circular Graph of Type 1 Statistical Analysis: Cosmonauts (Post vs Pre) using Linear Mixed Models

names = [str(labels[i])[10:] for i in range(len(labels))]
node_order = list()
node_order.extend(names)
node_angles = circular_layout(
    names, node_order, start_pos=90, group_boundaries=[0, len(names) / 2]
)
fig, ax = plt.subplots(figsize=(15, 15), facecolor='white', subplot_kw=dict(polar=True))

A = np.zeros((100, 100))
for _, row in sig_edges.iterrows():
    n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
    weight = row["beta"]
    A[n1, n2] = weight
    A[n2, n1] = weight

A[A==0] = np.nan

plot_connectivity_circle(
    A,
    names,
    node_angles=node_angles,
    node_colors=node_colors,
    title="",
    ax=ax,
    linewidth=2,
    vmin=-0.26, 
    vmax=0.26, 
    colorbar=True,
    colormap='coolwarm',
    fontsize_names = 0,
    fontsize_colorbar=12,
    facecolor='white',
    textcolor='black'
)
fig.tight_layout()
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/FC/MassUnivariate/Statistics_Cosmonauts_LMM_PostPre.pdf'
    ),
    dpi=300
)
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/FC/MassUnivariate/Statistics_Cosmonauts_LMM_PostPre.png'
    ),
    dpi=300
)

#%% Plotting the results of Type 2 Statistical Analysis: Controls (Post vs Pre) using Linear Mixed Models
atlas = 'SCH100'

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
labels = atlas_file.labels[1:]  # Labels from "1" to "R"
R = len(labels)
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
    f'{atlas}/FC/MassUnivariate/Statistics_Controls_LMM_PostPre.xlsx'
)
df = pd.read_excel(file_addr)

# Filter significant edges
sig_edges = df[df["p_fdr"] < 0.05].copy()

# Creating node colors 
if atlas=='SCH100':
    roi_to_network = {i: net for i, net in enumerate(
        ['Visual'] * 9 + ['Somatomotor'] * 6 + ['DorsalAttention'] * 8 +
        ['VentralAttention'] * 7 + ['Limbic'] * 3 + ['Control'] * 4 + ['DMN'] * 13 + 
        ['Visual'] * 8 + ['Somatomotor'] * 8 + ['DorsalAttention'] * 7 +
        ['VentralAttention'] * 5 + ['Limbic'] * 2 + ['Control'] * 9 + ['DMN'] * 11
    )}
    network_colors = {
        'Visual': 'brown',
        'Somatomotor': 'green',
        'DorsalAttention': 'orange',
        'VentralAttention': 'olive',
        'Limbic': 'black',
        'Control': 'cyan',
        'DMN': 'magenta'
    }
    node_colors = [network_colors[roi_to_network[i]] for i in range(len(roi_to_network))]

if len(sig_edges) == 0:
    print('No significant edge was found!')
    df[["node1", "node2"]] = df["edge"].str.split("-", expand=True)
    df = df[df["node1"].isin(region_coords) & df["node2"].isin(region_coords)]
    # Graph Creation
    A = np.zeros((100, 100))
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids]
    view = plotting.view_connectome(
        A, 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser()

else:
    sig_edges[["node1", "node2"]] = sig_edges["edge"].str.split("-", expand=True)
    sig_edges = sig_edges[
        sig_edges["node1"].isin(region_coords) & sig_edges["node2"].isin(region_coords)
    ]
    # Graph Creation
    A = np.zeros((100, 100))
    for _, row in sig_edges.iterrows():
        n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
        weight = row["beta"]
        A[n1, n2] = weight
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids] 
    view = plotting.view_connectome(
        A, 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser()

#%% Plot Circular Graph of Type 2 Statistical Analysis: Controls (Post vs Pre) using Linear Mixed Models

names = [str(labels[i])[10:] for i in range(len(labels))]
node_order = list()
node_order.extend(names)
node_angles = circular_layout(
    names, node_order, start_pos=90, group_boundaries=[0, len(names) / 2]
)
fig, ax = plt.subplots(figsize=(15, 15), facecolor='white', subplot_kw=dict(polar=True))

A = np.zeros((100, 100))
for _, row in sig_edges.iterrows():
    n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
    weight = row["beta"]
    A[n1, n2] = weight
    A[n2, n1] = weight

A[A==0] = np.nan

plot_connectivity_circle(
    A,
    names,
    node_angles=node_angles,
    node_colors=node_colors,
    title="",
    ax=ax,
    linewidth=2,
    vmin=-0.26, 
    vmax=0.26, 
    colorbar=True,
    colormap='coolwarm',
    fontsize_names = 0,
    fontsize_colorbar=12,
    facecolor='white',
    textcolor='black'
)
fig.tight_layout()
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/FC/MassUnivariate/Statistics_Controls_LMM_PostPre.pdf'
    ),
    dpi=300
)
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/FC/MassUnivariate/Statistics_Controls_LMM_PostPre.png'
    ),
    dpi=300
)

#%% Plotting the results of Type 3 Statistical Analysis: Cosmonauts (Post - Pre) vs Controls (Post - Pre) using Linear Mixed Models

atlas = 'SCH100'

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
labels = atlas_file.labels[1:]  # Labels from "1" to "R"
R = len(labels)
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
    f'{atlas}/FC/MassUnivariate/Statistics_Combined_LMM_PostPre.xlsx'
)
df = pd.read_excel(file_addr)

# Filter significant edges
sig_edges = df[df["p_fdr_interaction"] < 0.05].copy()

# Creating node colors 
if atlas=='SCH100':
    roi_to_network = {i: net for i, net in enumerate(
        ['Visual'] * 9 + ['Somatomotor'] * 6 + ['DorsalAttention'] * 8 +
        ['VentralAttention'] * 7 + ['Limbic'] * 3 + ['Control'] * 4 + ['DMN'] * 13 + 
        ['Visual'] * 8 + ['Somatomotor'] * 8 + ['DorsalAttention'] * 7 +
        ['VentralAttention'] * 5 + ['Limbic'] * 2 + ['Control'] * 9 + ['DMN'] * 11
    )}
    network_colors = {
        'Visual': 'brown',
        'Somatomotor': 'green',
        'DorsalAttention': 'orange',
        'VentralAttention': 'olive',
        'Limbic': 'black',
        'Control': 'cyan',
        'DMN': 'magenta'
    }
    node_colors = [network_colors[roi_to_network[i]] for i in range(len(roi_to_network))]

if len(sig_edges) == 0:
    print('No significant edge was found!')
    df[["node1", "node2"]] = df["edge"].str.split("-", expand=True)
    df = df[df["node1"].isin(region_coords) & df["node2"].isin(region_coords)]
    # Graph Creation
    A = np.zeros((100, 100))
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids]
    view = plotting.view_connectome(
        A, 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser() 

else:
    sig_edges[["node1", "node2"]] = sig_edges["edge"].str.split("-", expand=True)
    sig_edges = sig_edges[
        sig_edges["node1"].isin(region_coords) & sig_edges["node2"].isin(region_coords)
    ]
    # Graph Creation
    A = np.zeros((100, 100))
    for _, row in sig_edges.iterrows():
        n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
        weight = row["beta_interaction"]
        A[n1, n2] = weight
    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids] 
    view = plotting.view_connectome(
        A, 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser()

#%% Plot Circular Graph of Type 3 Statistical Analysis: Cosmonauts (Post - Pre) vs Controls (Post - Pre) using Linear Mixed Models

names = [str(labels[i])[10:] for i in range(len(labels))]
node_order = list()
node_order.extend(names)
node_angles = circular_layout(
    names, node_order, start_pos=90, group_boundaries=[0, len(names) / 2]
)
fig, ax = plt.subplots(figsize=(15, 15), facecolor='white', subplot_kw=dict(polar=True))

A = np.zeros((100, 100))
for _, row in sig_edges.iterrows():
    n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
    weight = row["beta_interaction"]
    A[n1, n2] = weight
    A[n2, n1] = weight

A[A==0] = np.nan

plot_connectivity_circle(
    A,
    names,
    node_angles=node_angles,
    node_colors=node_colors,
    title="",
    ax=ax,
    linewidth=2,
    vmin=-0.8, 
    vmax=0.8, 
    colorbar=True,
    colormap='coolwarm',
    fontsize_names = 0,
    fontsize_colorbar=12,
    facecolor='white',
    textcolor='black'
)
fig.tight_layout()
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/SC/MassUnivariate/Statistics_Combined_LMM_PostPre.pdf'
    ),
    dpi=300
)
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/SC/MassUnivariate/Statistics_Combined_LMM_PostPre.png'
    ),
    dpi=300
)
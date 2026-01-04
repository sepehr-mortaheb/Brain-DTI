#%% Importing

import os
import numpy as np 
import pandas as pd 
import os.path as op
import matplotlib.pylab as plt 
import seaborn as sns
import nibabel as nb
from nilearn import datasets, image
from scipy.ndimage import center_of_mass

#%% Parameters and directories initialization 
data_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

# %% Visualization of Global Metrics for Cosmonauts and Controls

atlas = 'SCH100'
norm = False

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

# Reading the data for Cosmonauts 
df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_global_Cosmonauts_auc_{atlas}_norm-{norm}.xlsx'))
df_ctrl = pd.read_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_global_Controls_auc_{atlas}_norm-{norm}.xlsx'))
df_cosm['group'] = 'cosmonaut'
df_ctrl['group'] = 'control'
df_ctrl.loc[(df_ctrl['time'] == 1), 'time'] = 'pre2'
df_ctrl.loc[(df_ctrl['time'] == 2), 'time'] = 'post'
df_ctrl['flight'] = 'f1'
df = pd.concat((df_cosm, df_ctrl), ignore_index=True)
df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])

metrics = [
    'global_efficiency',
    'char_path_length',
    'density',
    'transitivity',
    'modularity'
]

for metric in metrics:
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    sns.violinplot(data=df[df.metric==metric], 
                   x="group", 
                   y="auc", 
                   hue="time",
                   split=True, 
                   inner="quart", 
                   fill=True,
                   alpha=0.1,
                   palette={"pre2": "g", "post": "r"},
                   width=0.4,
                   ax=ax)

    ax.grid(True)
    ax.set_xlabel('', size=15)
    ax.set_ylabel('AUC', size=12)
    fig.suptitle(metric, size=15)

    fig.savefig(
        op.join(res_dir, f'{atlas}/SC/Graph/Graph_{metric}_{atlas}_norm-{norm}.pdf'),
        dpi=300,
        bbox_inches='tight'
    )
    fig.savefig(
        op.join(res_dir, f'{atlas}/SC/Graph/Graph_{metric}_{atlas}_norm-{norm}.png'),
        dpi=300,
        bbox_inches='tight'
    )

#%% Visualization of Nodal Metrics Statistical Results Type 1: Cosmonauts

atlas = 'SCH100'
norm = False
metric = 'clustering'

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
    f'{atlas}/SC/Graph/Statistics_Cosmomnauts_nodalAUC_norm-{norm}.xlsx'
)

df = pd.read_excel(file_addr)
df_metric = df[df.metric==metric]

# Filter significant nodes
sig_nodes = df_metric[df_metric["p_fdr"] < 0.05].copy()

if len(sig_nodes) == 0:
    print('No significant node was found!')

else:
    
    # image Creation
    A = np.zeros((100, 1))
    for _, row in sig_nodes.iterrows():
        n = int(row["node"][1:])-1
        beta = row["beta"]
        A[n] = beta
    img = atlas_img
    img_data = atlas_data
    for i in range(100):
        img_data[img_data==i+1] = A[i]
    
    img = nb.Nifti1Image(img_data, img.affine, img.header)

    nb.save(
        img,
        op.join(
            res_dir,
            f'{atlas}/SC/Graph/Statistics_Cosmomnauts_nodalAUC_{metric}_norm-{norm}.nii'
        )
    )

    print(sig_nodes)

#%% Visualization of Nodal Metrics Statistical Results Type 2: Controls

atlas = 'SCH100'
norm = False
metric = 'clustering'

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
    f'{atlas}/SC/Graph/Statistics_Controls_nodalAUC_norm-{norm}.xlsx'
)

df = pd.read_excel(file_addr)
df_metric = df[df.metric==metric]

# Filter significant nodes
sig_nodes = df_metric[df_metric["p_fdr"] < 0.05].copy()

if len(sig_nodes) == 0:
    print('No significant node was found!')

else:
    
    # image Creation
    A = np.zeros((100, 1))
    for _, row in sig_nodes.iterrows():
        n = int(row["node"][1:])-1
        beta = row["beta"]
        A[n] = beta
    img = atlas_img
    img_data = atlas_data
    for i in range(100):
        img_data[img_data==i+1] = A[i]
    
    img = nb.Nifti1Image(img_data, img.affine, img.header)

    nb.save(
        img,
        op.join(
            res_dir,
            f'{atlas}/SC/Graph/Statistics_Controls_nodalAUC_{metric}_norm-{norm}.nii'
        )
    )

    print(sig_nodes)

#%% Visualization of Nodal Metrics Statistical Results Type 3: Cosmonauts vs Controls interaction

atlas = 'SCH100'
norm = False
metric = 'clustering'

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
    f'{atlas}/SC/Graph/Statistics_Cosmonauts_vs_Controls_nodalAUC_norm-{norm}.xlsx'
)

df = pd.read_excel(file_addr)
df_metric = df[df.metric==metric]

# Filter significant nodes
sig_nodes = df_metric[df_metric["p_fdr_interaction"] < 0.05].copy()

if len(sig_nodes) == 0:
    print('No significant node was found!')

else:
    
    # image Creation
    A = np.zeros((100, 1))
    for _, row in sig_nodes.iterrows():
        n = int(row["node"][1:])-1
        beta = row["beta"]
        A[n] = beta
    img = atlas_img
    img_data = atlas_data
    for i in range(100):
        img_data[img_data==i+1] = A[i]
    
    img = nb.Nifti1Image(img_data, img.affine, img.header)

    nb.save(
        img,
        op.join(
            res_dir,
            f'{atlas}/SC/Graph/Statistics_Cosmonauts_vs_Controls_nodalAUC_{metric}_norm-{norm}.nii'
        )
    )

    print(sig_nodes)

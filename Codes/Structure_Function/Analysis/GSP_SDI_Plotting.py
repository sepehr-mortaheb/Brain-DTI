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
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import nibabel as nb

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

#%% Visualization of the SDI results for Cosmonauts

atlas = 'SCH100' # SCH100, SCH400, AAL
metric = 'logsdi'

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

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
    f'{atlas}/GSP/SDI/Statistics_SDI_Cosmomnauts_{metric}.xlsx'
)

df = pd.read_excel(file_addr)

# Filter significant nodes
sig_nodes = df[df["p_fdr"] < 0.05].copy()

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
            f'{atlas}/GSP/SDI/Statistics_SDI_Cosmomnauts_{metric}.nii'
        )
    )

    print(sig_nodes)

#%% Visualization of the SDI results for Cosntrols

atlas = 'SCH100' # SCH100, SCH400, AAL
metric = 'logsdi'  

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

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
    f'{atlas}/GSP/SDI/Statistics_SDI_Controls_{metric}.xlsx'
)
df = pd.read_excel(file_addr)

# Filter significant nodes
sig_nodes = df[df["p_fdr"] < 0.05].copy()
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
            f'{atlas}/GSP/SDI/Statistics_SDI_Controls_{metric}.nii'
        )
    )

    print(sig_nodes)


#%% Visualization of the SDI results for Cosmonauts vs Controls

atlas = 'SCH100' # SCH100, SCH400, AAL
metric = 'logsdi'

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

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
    f'{atlas}/GSP/SDI/Statistics_SDI_Cosmonauts_vs_Controls_{metric}.xlsx'
)
df = pd.read_excel(file_addr)

# Filter significant nodes
sig_nodes = df[df["p_fdr_interaction"] < 0.05].copy()
if len(sig_nodes) == 0:
    print('No significant node was found!')
else:
    # image Creation
    A = np.zeros((100, 1))
    for _, row in sig_nodes.iterrows():
        n = int(row["node"][1:])-1
        beta = row["beta_interaction"]
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
            f'{atlas}/GSP/SDI/Statistics_SDI_Cosmonauts_vs_Controls_{metric}.nii'
        )
    )
    print(sig_nodes)



#%% 

rs_file = op.join(data_dir, 'sub-cosmonaut01', 'ses-f1pre2', f'ts_{atlas}.mat')
rs = loadmat(rs_file)[f'ts_{atlas}']
rs = rs.T

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

vec = rs[:, 40]
th = np.max(np.abs(vec))
fig, ax = plt.subplots(1, 1, figsize=(12,4))
plotting.plot_markers(
    vec, 
    region_coords, 
    node_size=40,
    node_vmin=-th,
    node_vmax=th,
    node_cmap='RdBu_r',
    axes=ax
)
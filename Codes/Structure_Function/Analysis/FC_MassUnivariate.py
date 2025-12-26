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

data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

#%% Reading Cosmonauts' and Controls' FC files and dataframe creation

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

            fc = loadmat(
                op.join(data_dir, sub, ses, f'FC_{atlas}.mat'), 
            )['FC']

            for row in range(R):
                for col in range(row+1,R):
                    edge = f'R{row+1}-R{col+1}'
                    val = fc[row, col]
                    dftmp = pd.DataFrame([])
                    dftmp['subject'] = [sub]
                    dftmp['flight'] = [flight]
                    dftmp['time'] = [tp]
                    dftmp['edge'] = [edge]
                    dftmp['val'] = [val]
                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_csv(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_Cosmonauts.csv'))

df = pd.DataFrame([])
for sub in subjects:
    if sub.startswith('sub-control'):
        print(sub)
        for ses in sessions[sub]:
            print(f'---- {ses}')
            tp = ses.split('-')[1]

            fc = loadmat(
                op.join(data_dir, sub, ses, f'FC_{atlas}.mat'), 
            )['FC']

            for row in range(R):
                for col in range(row+1,R):
                    edge = f'R{row+1}-R{col+1}'
                    val = fc[row, col]
                    dftmp = pd.DataFrame([])
                    dftmp['subject'] = [sub]
                    dftmp['time'] = [tp]
                    dftmp['edge'] = [edge]
                    dftmp['val'] = [val]
                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_csv(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_Controls.csv'))

#%% Statistical Analysis 

# For cosmonauts 

atlas = 'SCH100'

df = pd.read_csv(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_Cosmonauts.csv'))

df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])
df['fish_val'] = np.arctanh(df['val'])

results = []
for edge in tqdm(df['edge'].unique(), desc="Running LMMs per edge"):

    df_edge = df[df['edge'] == edge]

    model = smf.mixedlm("fish_val ~ time", df_edge, groups=df_edge["subject"])
    result = model.fit(reml=False)

    coef = result.params.get("time[T.post]", float("nan"))
    pval = result.pvalues.get("time[T.post]", float("nan"))
    tval = result.tvalues.get("time[T.post]", float("nan"))

    results.append({
        "edge": edge,
        "beta": coef,
        "t_value": tval,
        "p_value": pval
    })

# ---- Compile results ----
df_results = pd.DataFrame(results)

valid_mask = df_results['p_value'].notna()
df_results['p_fdr'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr'] = corrected_pvals
df_results['significant'] = df_results['p_fdr'] < 0.05


df_results.to_excel(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_MassUnivariate_LMM.xlsx'), index=False)


#%% Plotting Results

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
    f'{atlas}/FC/MassUnivariate/FC_MassUnivariate_LMM.xlsx'
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
    for _, row in sig_edges.iterrows():
        n1, n2 = int(row["node1"][1:])-1, int(row["node2"][1:])-1
        weight = 0
        A[n1, n2] = weight

    # Plotting the Graph on the Brain Atlas 
    coords = [region_coords[n] for n in region_ids]
    view = plotting.view_connectome(
        nx.to_numpy_array(G, weight="weight"), 
        coords, 
        node_size=7,
        edge_cmap="bwr",
        node_color=node_colors,
        linewidth=5
    )
    view.open_in_browser()   

else:
    sig_edges[["node1", "node2"]] = sig_edges["edge"].str.split("-", expand=True)
    sig_edges['node11'] = [f'R{int(i[1:])+1}' for i in sig_edges['node1']]
    sig_edges['node22'] = [f'R{int(i[1:])+1}' for i in sig_edges['node2']]
    sig_edges = sig_edges[
        sig_edges["node11"].isin(region_coords) & sig_edges["node22"].isin(region_coords)
    ]
    # Graph Creation
    A = np.zeros((100, 100))
    for _, row in sig_edges.iterrows():
        n1, n2 = int(row["node11"][1:])-1, int(row["node22"][1:])-1
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

#%% Plot Circular Graph

names = [str(labels[i])[12:] for i in range(len(labels))]
node_order = list()
node_order.extend(names)
node_angles = circular_layout(
    names, node_order, start_pos=90, group_boundaries=[0, len(names) / 2]
)
fig, ax = plt.subplots(figsize=(15, 15), facecolor='white', subplot_kw=dict(polar=True))

A = np.zeros((100, 100))
for _, row in sig_edges.iterrows():
    n1, n2 = int(row["node11"][1:])-1, int(row["node22"][1:])-1
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
        f'{atlas}/FC/MassUnivariate/Statistics_Cosmonauts_LMM.pdf'
    ),
    dpi=300
)
fig.savefig(
    op.join(
        res_dir,
        f'{atlas}/FC/MassUnivariate/Statistics_Cosmonauts_LMM.png'
    ),
    dpi=300
)

#%% Network Analysis 

atlas = 'SCH100'

df = pd.read_csv(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_Cosmonauts.csv'))

df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['fish_val'] = np.arctanh(df['val'])
df[['node1', 'node2']] = df['edge'].str.split("-", expand=True)
df['node11'] = [f'R{int(i[1:])+1}' for i in df['node1']]
df['node22'] = [f'R{int(i[1:])+1}' for i in df['node2']]
    

networks = {
        'R1':'Vis',
        'R2':'Vis',
        'R3':'Vis',
        'R4':'Vis',
        'R5':'Vis',
        'R6':'Vis',
        'R7':'Vis',
        'R8':'Vis',
        'R9':'Vis',
        'R51':'Vis',
        'R52':'Vis',
        'R53':'Vis',
        'R54':'Vis',
        'R55':'Vis',
        'R56':'Vis',
        'R57':'Vis',
        'R58':'Vis',
        'R10':'SomMot',
        'R11':'SomMot',
        'R12':'SomMot',
        'R13':'SomMot',
        'R14':'SomMot',
        'R15':'SomMot',
        'R59':'SomMot',
        'R60':'SomMot',
        'R61':'SomMot',
        'R62':'SomMot',
        'R63':'SomMot',
        'R64':'SomMot',
        'R65':'SomMot',
        'R66':'SomMot',
        'R16':'DAN',
        'R17':'DAN',
        'R18':'DAN',
        'R19':'DAN',
        'R20':'DAN',
        'R21':'DAN',
        'R22':'DAN',
        'R23':'DAN',
        'R67':'DAN',
        'R68':'DAN',
        'R69':'DAN',
        'R70':'DAN',
        'R71':'DAN',
        'R72':'DAN',
        'R73':'DAN',
        'R24':'Sal',
        'R25':'Sal',
        'R26':'Sal',
        'R27':'Sal',
        'R28':'Sal',
        'R29':'Sal',
        'R30':'Sal',
        'R74':'Sal',
        'R75':'Sal',
        'R76':'Sal',
        'R77':'Sal',
        'R78':'Sal',
        'R31':'Lim',
        'R32':'Lim',
        'R33':'Lim',
        'R79':'Lim',
        'R80':'Lim',
        'R34':'Cont',
        'R35':'Cont',
        'R36':'Cont',
        'R37':'Cont',
        'R81':'Cont',
        'R82':'Cont',
        'R83':'Cont',
        'R84':'Cont',
        'R85':'Cont',
        'R86':'Cont',
        'R87':'Cont',
        'R88':'Cont',
        'R89':'Cont',
        'R38':'DMN',
        'R39':'DMN',
        'R40':'DMN',
        'R41':'DMN',
        'R42':'DMN',
        'R43':'DMN',
        'R44':'DMN',
        'R45':'DMN',
        'R46':'DMN',
        'R47':'DMN',
        'R48':'DMN',
        'R49':'DMN',
        'R50':'DMN',
        'R90':'DMN',
        'R91':'DMN',
        'R92':'DMN',
        'R93':'DMN',
        'R94':'DMN',
        'R95':'DMN',
        'R96':'DMN',
        'R97':'DMN',
        'R98':'DMN',
        'R99':'DMN',
        'R100':'DMN'
}

net1 = [networks[i] for i in list(df['node11'])]
net2 = [networks[i] for i in list(df['node22'])]

df['net1'] = net1
df['net2'] = net2

final = []

for subj in df['subject'].unique():
    print(subj)
    for time in df['time'].unique(): 
        tmp = []
        for net1 in df['net1'].unique():
            for net2 in df['net2'].unique():
                if (not f'{net1}-{net2}' in tmp) | (not f'{net2}-{net1}' in tmp):
                    df_net = df[
                        (df['subject']==subj) &
                        (df['time']==time) &
                        (df['net1']==net1) & 
                        (df['net2']==net2)
                    ]
                    tmp.append(f'{net1}-{net2}')
                    tmp.append(f'{net2}-{net1}')
                    if len(df_net)>0:
                        mean_val = np.mean(df_net['fish_val'])
                        final.append({
                            "subject": subj,
                            "time": time,
                            "net1": net1,
                            "net2": net2,
                            "val": mean_val
                        })
df_final = pd.DataFrame(final)

df_final.to_excel(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_Cosmonauts_networks.xlsx'))

#%% Statistics for networks

df_final = pd.read_excel(op.join(res_dir, f'{atlas}/FC/MassUnivariate/FC_Cosmonauts_networks.xlsx'))
df_final['time'] = pd.Categorical(df_final['time'], categories=['pre2', 'post'])
results = []
for n in tmp[0:len(tmp)+1:2]:
    dftmp = df_final[(df_final['net1']==n.split('-')[0]) & (df_final['net2']==n.split('-')[1])]
    model = smf.mixedlm("val ~ time", dftmp, groups=dftmp["subject"])
    result = model.fit(reml=False)

    coef = result.params.get("time[T.post]", float("nan"))
    pval = result.pvalues.get("time[T.post]", float("nan"))
    tval = result.tvalues.get("time[T.post]", float("nan"))

    results.append({
    "net1": n.split('-')[0],
    "net2": n.split('-')[1],
    "beta": coef,
    "t_value": tval,
    "p_value": pval
})
df_results = pd.DataFrame(results)
df_results['p_fdr'] = multipletests(df_results['p_value'], method='fdr_bh')[1]
df_results['significance'] = df_results['p_fdr']<0.05

convert = {
    'Vis':0,
    'SomMot':1,
    'DAN':2,
    'Sal':3,
    'Lim':4,
    'Cont':5,
    'DMN':6
}

res = np.zeros((len(convert), len(convert)))
for n in tmp[0:len(tmp)+1:2]:
    dftmp = df_results[(df_results['net1']==n.split('-')[0]) & (df_results['net2']==n.split('-')[1])]
    if list(dftmp['p_value'] < 0.05)[0]:
        res[convert[n.split('-')[0]], convert[n.split('-')[1]]] = dftmp['beta']
vmin = np.min(res)
vmax = np.max(res)
lim = np.max(np.abs([vmin, vmax]))
betas = df_results.pivot(index="net1", columns="net2", values="beta")
ps = df_results.pivot(index="net1", columns="net2", values="p_value")
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
sns.heatmap(
    betas, 
    vmin=-lim, 
    vmax=lim,
    cmap='jet',
    annot=True,
    ax=ax
)
sns.heatmap(
    betas[ps>=0.05], 
    annot=True, 
    cbar=False, 
    cmap=sns.color_palette("Greys", n_colors=1, desat=1),
    ax=ax
)
ax.set_xlabel('')
ax.set_ylabel('')
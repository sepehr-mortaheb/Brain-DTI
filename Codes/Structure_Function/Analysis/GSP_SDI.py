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

# %% Preparing the structural harmonics info
norm = False
e = loadmat(op.join(res_dir, f'{atlas}/SC/Harmonics/SC_avg_cosm_{atlas}_norm-{norm}_eigenvalues.mat'))['eigVals'][0]
U = loadmat(op.join(res_dir, f'{atlas}/SC/Harmonics/SC_avg_cosm_{atlas}_norm-{norm}_eigenvectors.mat'))['eigVecs']
R = len(e)

#%% Calculate Filtering Cutoff Frequency 

TPSD = np.zeros((R,1))
count=0
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        for ses in sessions[sub]:
            flight = ses.split('-')[1][0:2]
            tp = ses.split('-')[1][2:]
            if (flight=='f1') & (tp in ['pre2', 'post']):
                count +=1
                rs_file = op.join(data_dir, sub, ses, f'ts_{atlas}.mat')
                rs = loadmat(rs_file)[f'ts_{atlas}']
                rs = rs.T
                Xhat = np.matmul(np.transpose(U), rs)
                Xhat = np.array(Xhat)
                PSD = Xhat**2 
                PSD = np.mean(PSD, axis=1)
                TPSD = np.concatenate((TPSD, PSD.reshape((R,1))), axis=1)
TPSD = TPSD[:, 1:]
MPSD = np.mean(TPSD, axis=1) 
SPSD = np.std(TPSD, axis=1)

auctot = np.trapz(MPSD)
auc = 0
c = 0
while auc<auctot/2: 
    print(f'{c}  {e[c]}  {auc} {auctot/2}')
    auc = np.trapz(MPSD[0:c])
    c = c+1
freq_cut = e[c-1]

#%% Plotting Filtering Cutoff Frequency

fig, ax = plt.subplots(1,1, figsize=(9,5))
ax.semilogy(e, MPSD, color='black', linewidth=2)
ax.fill_between(e, MPSD-SPSD, MPSD+SPSD, alpha=0.1, color='black')
ax.grid(True, 'minor')

ax.axvspan(0, freq_cut, alpha=0.2)
ax.axvspan(freq_cut, max(e), color='red', alpha=0.2)

ax.set_xlim([min(e), max(e)])
ax.set_xlabel('Eigenvalues', size=15)
ax.set_ylabel('Power Spectrum Density', size=15)

fig.savefig(op.join(res_dir, f'{atlas}/GSP/freq_cut_{atlas}_norm-{norm}.png'), dpi=300)
fig.savefig(op.join(res_dir, f'{atlas}/GSP/freq_cut_{atlas}_norm-{norm}.pdf'), dpi=300)

# %% GFT-based LPF and HPF
c = 41
Ulow = np.zeros(U.shape)
Ulow[:, 0:c+1] = U[:, 0:c+1]
Uhigh = np.zeros(U.shape)
Uhigh[:, c+1:R] = U[:, c+1:R]

for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        for ses in sessions[sub]:
            flight = ses.split('-')[1][0:2]
            tp = ses.split('-')[1][2:]
            if (flight=='f1') & (tp in ['pre2', 'post']):
                rs_file = op.join(data_dir, sub, ses, f'ts_{atlas}.mat')
                rs = loadmat(rs_file)[f'ts_{atlas}']
                rs = rs.T
                rshat = np.matmul(U.T, rs)
                xc = np.matmul(Ulow, rshat)
                xd = np.matmul(Uhigh, rshat)
                savemat(op.join(data_dir, sub, ses, f'GFT_xlow_{atlas}_c-{c}.mat'), {'xc':xc})
                savemat(op.join(data_dir, sub, ses, f'GFT_xhigh_{atlas}_c-{c}.mat'), {'xd':xd})

#%% SDI calculation

c=41

# for cosmonauts
df_sdi_cosmonauts = pd.DataFrame([])
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        for ses in sessions[sub]:
            flight = ses.split('-')[1][0:2]
            tp = ses.split('-')[1][2:]
            if (flight=='f1') & (tp in ['pre2', 'post']):
                xc = loadmat(op.join(data_dir, sub, ses, f'GFT_xlow_{atlas}_c-{c}.mat'))['xc']
                xd = loadmat(op.join(data_dir, sub, ses, f'GFT_xhigh_{atlas}_c-{c}.mat'))['xd']
                cind = np.linalg.norm(xc, axis=1)
                dind = np.linalg.norm(xd, axis=1)
                sdi = dind/cind
                for i in range(len(sdi)):
                    r=i+1
                    tmpdf = pd.DataFrame([])
                    tmpdf['subject'] = [sub]
                    tmpdf['session'] = [ses]
                    tmpdf['flight'] = [flight]
                    tmpdf['time'] = [tp]
                    tmpdf['region'] = [f'R{r}']
                    tmpdf['coupling'] = [cind[i]]
                    tmpdf['decoupling'] = [dind[i]]
                    tmpdf['sdi'] = [sdi[i]]
                    tmpdf['logsdi'] = [np.log2(sdi[i])]
                    df_sdi_cosmonauts = pd.concat((df_sdi_cosmonauts, tmpdf), ignore_index=True)

df_sdi_cosmonauts.to_excel(op.join(res_dir, f'{atlas}/GSP/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))


#%% Statistical Analysis

c=41
df = pd.read_excel(op.join(res_dir, f'{atlas}/GSP/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])

metric = 'logsdi'

results = []
for r in range(1, R+1):
    node = f'R{r}'
    df_node = df[df['region'] == node]
    model = smf.mixedlm(f"{metric} ~ time", df_node, groups=df_node["subject"])
    result = model.fit(reml=False)

    coef = result.params.get("time[T.post]", float("nan"))
    pval = result.pvalues.get("time[T.post]", float("nan"))
    tval = result.tvalues.get("time[T.post]", float("nan"))

    results.append({
        "node": node,
        "beta": coef,
        "t_value": tval,
        "p_value": pval,
        "metric": metric
    })

df_results = pd.DataFrame(results)

valid_mask = df_results['p_value'].notna()
df_results['p_fdr'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr'] = corrected_pvals
df_results['significant'] = df_results['p_fdr'] < 0.05

df_results.to_excel(op.join(res_dir, f'{atlas}/GSP/Statistics_Cosmomnauts_{metric}.xlsx'), index=False)

#%% Visualization of the results 

atlas = 'SCH100'
metric = 'sdi'

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
    f'{atlas}/GSP/Statistics_Cosmomnauts_{metric}.xlsx'
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
            f'{atlas}/GSP/Statistics_Cosmomnauts_{metric}.nii'
        )
    )

    print(sig_nodes)

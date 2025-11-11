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
# %%
data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'
res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results'
mni_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/MNI'

atlas = 'SCH100'
subjects = np.sort(os.listdir(data_dir))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

# %% Preparing the structural harmonic info

e = loadmat(op.join(res_dir, f'SC/SC_eig_{atlas}.mat'))['e'][0]
U = loadmat(op.join(res_dir, f'SC/SC_eig_{atlas}.mat'))['U']
R = len(e)

# %% Cutoff frequency calculation

TPSD = np.zeros((R,1))
for sub in subjects:
    for ses in sessions[sub]:
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
    auc = np.trapz(MPSD[0:c])
    c = c+1
freq_cut = e[c-1]

# Plotting
fig, ax = plt.subplots(1,1, figsize=(9,5))
ax.semilogy(e, MPSD, color='black', linewidth=2)
ax.fill_between(e, MPSD-SPSD, MPSD+SPSD, alpha=0.1, color='black')
ax.grid(True, 'minor')

ax.axvspan(0, freq_cut, alpha=0.2)
ax.axvspan(freq_cut, max(e), color='red', alpha=0.2)

ax.set_xlim([min(e), max(e)])
ax.set_xlabel('Eigenvalues', size=15)
ax.set_ylabel('Power Spectrum Density', size=15)

fig.savefig(op.join(res_dir, f'GSP/freq_cut_{atlas}.png'), dpi=300)
fig.savefig(op.join(res_dir, f'GSP/freq_cut_{atlas}.pdf'), dpi=300)

# %% GFT-based LPF and HPF
c = 21
Ulow = np.zeros(U.shape)
Ulow[:, 0:c-1] = U[:, 0:c-1]
Uhigh = np.zeros(U.shape)
Uhigh[:, c-1:R] = U[:, c-1:R]

for sub in subjects:
    for ses in sessions[sub]:
        rs_file = op.join(data_dir, sub, ses, f'ts_{atlas}.mat')
        rs = loadmat(rs_file)[f'ts_{atlas}']
        rs = rs.T
        rshat = np.matmul(U.T, rs)
        xc = np.matmul(Ulow, rshat)
        xd = np.matmul(Uhigh, rshat)
        savemat(op.join(data_dir, sub, ses, f'GFT_xlow_{atlas}_c-{c}.mat'), {'xc':xc})
        savemat(op.join(data_dir, sub, ses, f'GFT_xhigh_{atlas}_c-{c}.mat'), {'xd':xd})

#%% SDI calculation

c=21

# for cosmonauts
df_sdi_cosmonauts = pd.DataFrame([])
for sub in subjects:
    if ('cosmonaut' in sub):
        for ses in sessions[sub]:
            xc = loadmat(op.join(data_dir, sub, ses, f'GFT_xlow_{atlas}_c-{c}.mat'))['xc']
            xd = loadmat(op.join(data_dir, sub, ses, f'GFT_xhigh_{atlas}_c-{c}.mat'))['xd']
            cind = np.linalg.norm(xc, axis=1)
            dind = np.linalg.norm(xd, axis=1)
            sdi = dind/cind
            fnum = ses.split('-')[1][0:2]
            tp = ses.split('-')[1][2:]
            for i in range(len(sdi)):
                r=i+1
                tmpdf = pd.DataFrame([])
                tmpdf['subject'] = [sub]
                tmpdf['session'] = [ses]
                tmpdf['flight'] = [fnum]
                tmpdf['time'] = [tp]
                tmpdf['region'] = [f'R{r}']
                tmpdf['coupling'] = [cind[i]]
                tmpdf['decoupling'] = [dind[i]]
                tmpdf['sdi'] = [sdi[i]]
                tmpdf['logsdi'] = [np.log2(sdi[i])]
                df_sdi_cosmonauts = pd.concat((df_sdi_cosmonauts, tmpdf), ignore_index=True)

df_sdi_cosmonauts.to_excel(op.join(res_dir, f'GSP/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))

# for controls 
df_sdi_controls = pd.DataFrame([])
for sub in subjects:
    if ('control' in sub):
        for ses in sessions[sub]:
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
                tmpdf['region'] = [f'R{r}']
                tmpdf['coupling'] = [cind[i]]
                tmpdf['decoupling'] = [dind[i]]
                tmpdf['sdi'] = [sdi[i]]
                tmpdf['logsdi'] = [np.log2(sdi[i])]
                df_sdi_controls = pd.concat((df_sdi_controls, tmpdf), ignore_index=True)

df_sdi_controls.to_excel(op.join(res_dir, f'GSP/df_sdi_controls_{atlas}_c-{c}.xlsx'))

# %%

c=21
R=100
param='logsdi'
df = pd.read_excel(op.join(res_dir, f'GSP/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))


dftmp = df[(df.flight=='f1') & (df.time!='pre1')]

for i in range(R): 
    r=i+1

    fig, ax = plt.subplots(
        1,
        1,
        figsize=(6,4)
    )
    sns.stripplot(
        data=dftmp[dftmp.region==f'R{r}'],
        x='time',
        y=param,
        order=['pre2', 'post', 'foll'],
        ax=ax
    )
    sns.pointplot(
        data=dftmp[dftmp.region==f'R{r}'],
        x='time',
        y=param,
        order=['pre2', 'post', 'foll'],
        ax=ax,
        color='black'
    )
    ax.grid(True)

    fig.savefig(op.join(res_dir, f'GSP/region_sdi_figs/{atlas}/{param}/{param}_R{r}_c-{c}.png'), dpi=300)
    fig.savefig(op.join(res_dir, f'GSP/region_sdi_figs/{atlas}/{param}/{param}_R{r}_c-{c}.pdf'), dpi=300)

# %%

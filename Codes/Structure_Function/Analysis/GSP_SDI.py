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

# %% Preparing the structural harmonics info

e = loadmat(op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_eigenvalues.mat'))['eigVals'][0]
U = loadmat(op.join(res_dir, f'SC/Harmonics/SC_avg_cosm_{atlas}_eigenvectors.mat'))['eigVecs']
R = len(e)

#%% Calculate Filtering Cutoff Frequency 

exc_list = [
  "sub-cosmonaut13", 
  "sub-cosmonaut15", 
  "sub-cosmonaut16",
  "sub-cosmonaut26"
]

TPSD = np.zeros((R,1))
count=0
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        if sub not in exc_list:
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

fig.savefig(op.join(res_dir, f'GSP/freq_cut_{atlas}.png'), dpi=300)
fig.savefig(op.join(res_dir, f'GSP/freq_cut_{atlas}.pdf'), dpi=300)

#%% GFT-base rs-fMRI filtering

# %% GFT-based LPF and HPF
c = 20
Ulow = np.zeros(U.shape)
Ulow[:, 0:c+1] = U[:, 0:c+1]
Uhigh = np.zeros(U.shape)
Uhigh[:, c+1:R] = U[:, c+1:R]

for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        if sub not in exc_list:
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

c=20

# for cosmonauts
df_sdi_cosmonauts = pd.DataFrame([])
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        if sub not in exc_list:
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

df_sdi_cosmonauts.to_excel(op.join(res_dir, f'GSP/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))

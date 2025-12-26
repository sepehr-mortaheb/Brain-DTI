#%% Importing 

import os.path as op 
import numpy as np 
from scipy.io import loadmat, savemat
from nilearn.connectome import ConnectivityMeasure
import os
import seaborn as sns
import matplotlib.pylab as plt 
from scipy.stats import zscore

#%% Parameters and Directories Initialization

data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc/ROS'

subjects = np.sort(os.listdir(op.join(data_dir)))
# to remove .DS files in mac
subjects = [sub for sub in subjects if sub.startswith('sub')] 
sessions = {sub:np.sort(os.listdir(op.join(data_dir, sub))) for sub in subjects}
# to remove .DS files in mac
sessions = {sub:[ses for ses in sessions[sub] if ses.startswith('ses')] for sub in subjects}

#%% Reading Cosmonauts' and Controls' RS time series and creating FC matrices 

atlas = 'SCH100' # SCH100, SCH400, AAL

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

# Cosmonauts 
conn_estimator = ConnectivityMeasure(kind='correlation')
for sub in subjects:
    if sub.startswith('sub-cosmonaut'):
        print(sub)
        for ses in sessions[sub]:
            print(f'---- {ses}')
            flight = ses.split('-')[1][0:2]
            tp = ses.split('-')[1][2:]

            ts = loadmat(
                op.join(data_dir, f'{sub}/{ses}/ts_{atlas}.mat')
            )[f'ts_{atlas}']

            fc_matrix = conn_estimator.fit_transform([ts])[0]

            savemat(
                op.join(data_dir, f'{sub}/{ses}/FC_{atlas}.mat'), 
                {'FC':fc_matrix}
            )

# Controls 
conn_estimator = ConnectivityMeasure(kind='correlation')
for sub in subjects:
    if sub.startswith('sub-control'):
        print(sub)
        for ses in sessions[sub]:
            print(f'---- {ses}')
            tp = ses.split('-')[1]

            ts = loadmat(
                op.join(data_dir, f'{sub}/{ses}/ts_{atlas}.mat')
            )[f'ts_{atlas}']

            fc_matrix = conn_estimator.fit_transform([ts])[0]

            savemat(
                op.join(data_dir, f'{sub}/{ses}/FC_{atlas}.mat'), 
                {'FC':fc_matrix}
            )

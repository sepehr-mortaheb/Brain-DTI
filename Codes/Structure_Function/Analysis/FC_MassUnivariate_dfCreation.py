#%% Importing 
import os 
import os.path as op 
import pandas as pd 
import numpy as np 
from scipy.io import loadmat

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

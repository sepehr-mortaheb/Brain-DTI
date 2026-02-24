#%% importing libraries
import numpy as np
import matplotlib.pyplot as plt
import os.path as op
import pandas as pd
import seaborn as sns

#%% loading data

data_dir = '/Users/sepehrleia/Library/CloudStorage/GoogleDrive-sepmori2023@gmail.com/My Drive/Academic/LEIA/Projects/BRAIN-DTI/Joystick/dataframe.xlsx'
df = pd.read_excel(data_dir)
df['abs_tilt'] = np.abs(df['tilt'])

#%% plotting
sns.set_style('whitegrid')
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

df = df[df.flight == 'f1']

sns.pointplot(data=df[df.direction == 'CW'], x='session', y='abs_tilt', hue='experience', ax=ax[0])
sns.pointplot(data=df[df.direction == 'CCW'], x='session', y='abs_tilt', hue='experience', ax=ax[1])

ax[0].set_title('Clockwise (CW) Rotation')
ax[1].set_title('Counter-Clockwise (CCW) Rotation')

ax[0].set_xlabel('Session')
ax[1].set_xlabel('Session')
ax[0].set_ylabel('Absolute Tilt (degrees)')
ax[1].set_ylabel('Absolute Tilt (degrees)')
ax[0].legend().set_visible(False)
ax[1].legend(loc='upper left', bbox_to_anchor=(1, 1), title='Experience')
plt.tight_layout()
plt.show()

#%% plotting
sns.set_style('whitegrid')
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

df = df[df.flight == 'f1']

sns.pointplot(data=df[df.direction == 'CW'], x='session', y='abs_tilt', ax=ax[0])
sns.pointplot(data=df[df.direction == 'CCW'], x='session', y='abs_tilt', ax=ax[1])

ax[0].set_title('Clockwise (CW) Rotation')
ax[1].set_title('Counter-Clockwise (CCW) Rotation')

ax[0].set_xlabel('Session')
ax[1].set_xlabel('Session')
ax[0].set_ylabel('Absolute Tilt (degrees)')
ax[1].set_ylabel('Absolute Tilt (degrees)')
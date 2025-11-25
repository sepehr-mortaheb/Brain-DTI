#%% Importing 

import os 
import os.path as op 
import pandas as pd 
import seaborn as sns 
import matplotlib.pylab as plt

from scipy import stats

# %% Reading the demographic file 

data_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/data_ts_sc'
file_name = 'subjects_demographics.xlsx'

res_dir = '/Users/sepehrmortaheb/MyDrive/Academic/LEIA/Projects/BRAIN-DTI/Cosmonauts_StructFunc/Analysis/results/demog'

file_addr = op.join(
    data_dir,
    file_name
)

#%% Age Comparison (All Subjects)

df = pd.read_excel(file_addr, sheet_name='all')

fig, ax = plt.subplots(1, 1, figsize=(6,4))
sns.violinplot(
    x='group',
    y='age', 
    data=df,
    color='white',
    inner='quart',
    ax=ax,
)

sns.swarmplot(
    x='group',
    y='age', 
    data=df,
    ax=ax,
    hue='group',
    hue_order=['cosm', 'cont']
)
ax.grid(True)
ax.set_xlabel('Group', size=15)
ax.set_ylabel('Age at the first Scan', size=15)
fig.savefig(op.join(res_dir, f'age_dist.png'), dpi=300)
fig.savefig(op.join(res_dir, f'age_dist.pdf'), dpi=300)

# Stats 

pd.set_option('display.precision', 3)
pd.set_option('display.max_columns', 100)
display(df.groupby('group').describe())

stats.ttest_ind(
    df[df.group=='cosm']['age'], 
    df[df.group=='cont']['age']
)

#%% Cosmonauts Data Summary

df_cosm = pd.read_excel(file_addr, sheet_name='cosmonauts')
#df_cosm['f1post']=df_cosm['f1post']/365.25
#df_cosm['f1foll']=df_cosm['f1foll']/365.25
#df_cosm['f2pre']=df_cosm['f2pre']/365.25
#df_cosm['f2post']=df_cosm['f2post']/365.25
#df_cosm['f2foll']=df_cosm['f2foll']/365.25
#df_cosm['f3pre']=df_cosm['f3pre']/365.25
#df_cosm['f3post']=df_cosm['f3post']/365.25
#df_cosm['f3foll']=df_cosm['f3foll']/365.25

fig, ax = plt.subplots(1, 1, figsize=(9, 7))
sns.stripplot(
    x='f1pre', 
    y='subject',
    data=df_cosm, 
    color='blue',
    ax=ax,
    marker='o',
    size=4
)
sns.stripplot(
    x='f1post', 
    y='subject',
    data=df_cosm, 
    color='blue',
    ax=ax,
    marker='*',
    size=10
)
sns.stripplot(
    x='f1foll', 
    y='subject',
    data=df_cosm, 
    color='blue',
    ax=ax,
    marker='d',
    size=8
)

sns.stripplot(
    x='f2pre', 
    y='subject',
    data=df_cosm, 
    color='red',
    ax=ax,
    marker='o',
    size=4
)
sns.stripplot(
    x='f2post', 
    y='subject',
    data=df_cosm, 
    color='red',
    ax=ax,
    marker='*',
    size=10
)
sns.stripplot(
    x='f2foll', 
    y='subject',
    data=df_cosm, 
    color='red',
    ax=ax,
    marker='d',
    size=8
)

sns.stripplot(
    x='f3pre', 
    y='subject',
    data=df_cosm, 
    color='green',
    marker='o',
    size=4
)
sns.stripplot(
    x='f3post', 
    y='subject',
    data=df_cosm, 
    color='green',
    ax=ax,
    marker='*',
    size=10
)
sns.stripplot(
    x='f3foll', 
    y='subject',
    data=df_cosm, 
    color='green',
    ax=ax,
    marker='d',
    size=8
)

ax.grid(True)
ax.set_ylabel('Subjects', size=15)
ax.set_xlabel('Days after the first scan', size=15)

fig.savefig(op.join(res_dir, f'available_data_cosm.png'), dpi=300, bbox_inches='tight')
fig.savefig(op.join(res_dir, f'available_data_cosm.pdf'), dpi=300, bbox_inches='tight')
#%% Controls Data Summary

df_cont = pd.read_excel(file_addr, sheet_name='controls')
#df_cont['MRI1']=df_cont['MRI1']/365.25
#df_cont['MRI2']=df_cont['MRI2']/365.25
#df_cont['MRI3']=df_cont['MRI3']/365.25
#df_cont['MRI4']=df_cont['MRI4']/365.25
#df_cont['MRI5']=df_cont['MRI5']/365.25
#df_cont['MRI6']=df_cont['MRI6']/365.25
#df_cont['MRI7']=df_cont['MRI7']/365.25
#df_cont['MRI8']=df_cont['MRI8']/365.25

fig, ax = plt.subplots(1, 1, figsize=(9, 7))
colors = sns.color_palette('tab10', 8)
 
sns.stripplot(
    x='MRI1', 
    y='subject',
    data=df_cont, 
    color=colors[0],
    ax=ax
)
sns.stripplot(
    x='MRI2', 
    y='subject',
    data=df_cont, 
    color=colors[1],
    ax=ax
)
sns.stripplot(
    x='MRI3', 
    y='subject',
    data=df_cont, 
    color=colors[2],
    ax=ax
)
sns.stripplot(
    x='MRI4', 
    y='subject',
    data=df_cont, 
    color=colors[3],
    ax=ax
)
sns.stripplot(
    x='MRI5', 
    y='subject',
    data=df_cont, 
    color=colors[4],
    ax=ax
)
sns.stripplot(
    x='MRI6', 
    y='subject',
    data=df_cont, 
    color=colors[5],
    ax=ax
)
sns.stripplot(
    x='MRI7', 
    y='subject',
    data=df_cont, 
    color=colors[6],
    ax=ax
)
sns.stripplot(
    x='MRI8', 
    y='subject',
    data=df_cont, 
    color=colors[7],
    ax=ax
)

ax.grid(True)
ax.set_ylabel('Subjects', size=15)
ax.set_xlabel('Days after the first scan', size=15)

fig.savefig(op.join(res_dir, f'available_data_cont.png'), dpi=300, bbox_inches='tight')
fig.savefig(op.join(res_dir, f'available_data_cont.pdf'), dpi=300, bbox_inches='tight')
#%% Age Comparison (Selected Subjects)

df = pd.read_excel(file_addr, sheet_name='all')
exc_list = [
    'sub-cosmonaut13',
    'sub-cosmonaut15',
    'sub-cosmonaut16',
    'sub-cosmonaut26',
    'sub-cosmonaut19',
    'sub-cosmonaut20',
    'sub-control16',
    'sub-control18',
    'sub-control19',
    'sub-control22'
]
for sub in exc_list:
    df = df[df.subject!=sub]

df = df.reset_index()
df = df.drop('index', axis=1)

fig, ax = plt.subplots(1, 1, figsize=(6,4))
sns.violinplot(
    x='group',
    y='age', 
    data=df,
    color='white',
    inner='quart',
    ax=ax,
)

sns.swarmplot(
    x='group',
    y='age', 
    data=df,
    ax=ax,
    hue='group',
    hue_order=['cosm', 'cont']
)
ax.grid(True)
ax.set_xlabel('Group', size=15)
ax.set_ylabel('Age at the first Scan', size=15)
fig.savefig(op.join(res_dir, f'age_dist_sel.png'), dpi=300)
fig.savefig(op.join(res_dir, f'age_dist_sel.pdf'), dpi=300)

# Stats 

pd.set_option('display.precision', 3)
pd.set_option('display.max_columns', 100)
display(df.groupby('group').describe())

stats.ttest_ind(
    df[df.group=='cosm']['age'], 
    df[df.group=='cont']['age']
)

# %%

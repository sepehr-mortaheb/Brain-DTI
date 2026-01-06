#%% Importing

import os 
import os.path as op 
import pandas as pd 
import numpy as np 
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

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

#%% Statistical Analysis Type 1: Cosmonauts (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100' # SCH100, SCH400, AAL
c=41

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

df = pd.read_excel(op.join(res_dir, f'{atlas}/GSP/SDI/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))
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

df_results.to_excel(op.join(res_dir, f'{atlas}/GSP/SDI/Statistics_SDI_Cosmomnauts_{metric}.xlsx'), index=False)

#%% Statistical Analysis Type 2: Controls (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100' # SCH100, SCH400, AAL
c=40

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

df = pd.read_excel(op.join(res_dir, f'{atlas}/GSP/SDI/df_sdi_controls_{atlas}_c-{c}.xlsx'))
df['time'] = pd.Categorical(df['time'], categories=[1, 2])

metric = 'logsdi'

results = []
for r in range(1, R+1):
    node = f'R{r}'
    df_node = df[df['region'] == node]
    model = smf.mixedlm(f"{metric} ~ time", df_node, groups=df_node["subject"])
    result = model.fit(reml=False)

    coef = result.params.get("time[T.2]", float("nan"))
    pval = result.pvalues.get("time[T.2]", float("nan"))
    tval = result.tvalues.get("time[T.2]", float("nan"))

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

df_results.to_excel(op.join(res_dir, f'{atlas}/GSP/SDI/Statistics_SDI_Controls_{metric}.xlsx'), index=False)

#%% Global Metrics Statistical Analysis Type 3: Cosmonauts (Post - Pre) vs Controls (Post - Pre) using Linear Mixed Models

atlas = 'SCH100' # SCH100, SCH400, AAL
if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

c=41
df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/GSP/SDI/df_sdi_cosmonauts_{atlas}_c-{c}.xlsx'))

c=40
df_ctrl = pd.read_excel(op.join(res_dir, f'{atlas}/GSP/SDI/df_sdi_controls_{atlas}_c-{c}.xlsx'))

df_cosm['group'] = 'cosmonauts'
df_ctrl['group'] = 'controls'
df_ctrl.loc[(df_ctrl['time'] == 1), 'time'] = 'pre2'
df_ctrl.loc[(df_ctrl['time'] == 2), 'time'] = 'post'
df_ctrl['flight'] = 'f1'
df = pd.concat((df_cosm, df_ctrl), ignore_index=True)
df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])

metric = 'logsdi'

results = []
for r in range(1, R+1):
    node = f'R{r}'
    df_node = df[df['region'] == node]
    model = smf.mixedlm(f"{metric} ~ time * group", df_node, groups=df_node["subject"])
    result = model.fit(reml=False)

    coef_time = result.params.get("time[T.post]", float("nan"))
    coef_group = result.params.get("group[T.cosmonauts]", float("nan"))
    coef_interaction = result.params.get("time[T.post]:group[T.cosmonauts]", float("nan"))
    pval_time = result.pvalues.get("time[T.post]", float("nan"))
    pval_group = result.pvalues.get("group[T.cosmonauts]", float("nan"))
    pval_interaction = result.pvalues.get("time[T.post]:group[T.cosmonauts]", float("nan"))
    tval_time = result.tvalues.get("time[T.post]", float("nan"))
    tval_group = result.tvalues.get("group[T.cosmonauts]", float("nan"))
    tval_interaction = result.tvalues.get("time[T.post]:group[T.cosmonauts]", float("nan"))

    results.append({
        "node": node,
        "beta_time": coef_time,
        "beta_group": coef_group,
        "beta_interaction": coef_interaction,
        "t_value_time": tval_time,
        "t_value_group": tval_group,
        "t_value_interaction": tval_interaction,
        "p_value_time": pval_time,
        "p_value_group": pval_group,
        "p_value_interaction": pval_interaction,
        "metric": metric
    })

df_results = pd.DataFrame(results)

valid_mask = df_results['p_value_interaction'].notna()
df_results['p_fdr_interaction'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value_interaction'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr_interaction'] = corrected_pvals
df_results['significant_interaction'] = df_results['p_fdr_interaction'] < 0.05

df_results.to_excel(op.join(res_dir, f'{atlas}/GSP/SDI/Statistics_SDI_Cosmonauts_vs_Controls_{metric}.xlsx'), index=False)
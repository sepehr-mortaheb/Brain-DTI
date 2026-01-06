#%% Importing

import os
import numpy as np 
import pandas as pd 
import os.path as op
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

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

#%% Global Metrics Statistical Analysis Type 1: Cosmonauts (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'
graph = 'pos'

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))

df_cosm_filt = df_cosm
df_cosm_filt = df_cosm_filt[df_cosm_filt.flight=='f1']
df_cosm_filt = df_cosm_filt[(df_cosm_filt.time=='pre2') | (df_cosm_filt.time=='post')]
df_cosm_filt = df_cosm_filt.reset_index()
df_cosm_filt = df_cosm_filt.drop(axis='columns', labels=['index'])
df_cosm_filt['time'] = pd.Categorical(df_cosm_filt['time'], categories=['pre2', 'post'])

metric = 'modularity' # Possible metrics: 
                             # 'global_efficiency'
                             # 'char_path_length'
                             # 'modularity


df_used = df_cosm_filt[df_cosm_filt.metric==metric]
model = smf.mixedlm(f" auc ~ time", df_used, groups=df_used["subject"])
result = model.fit(reml=False)

coef = result.params.get("time[T.post]")
pval = result.pvalues.get("time[T.post]")
tval = result.tvalues.get("time[T.post]")

df_results = pd.DataFrame({
    'beta':[coef],
    't_value': [tval],
    'p_value': [pval]
})

df_results

#%% Global Metrics Statistical Analysis Type 2: Controls (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'
graph = 'pos'

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

df_ctrl = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Controls_auc_{atlas}_graph-{graph}.xlsx'))

df_ctrl_filt = df_ctrl
df_ctrl_filt = df_ctrl_filt[(df_ctrl_filt.time==1) | (df_ctrl_filt.time==2)]
df_ctrl_filt = df_ctrl_filt.reset_index()
df_ctrl_filt = df_ctrl_filt.drop(axis='columns', labels=['index'])
df_ctrl_filt['time'] = pd.Categorical(df_ctrl_filt['time'], categories=[1, 2])

metric = 'char_path_length' # Possible metrics: 
                             # 'global_efficiency'
                             # 'char_path_length'
                             # 'modularity


df_used = df_ctrl_filt[df_ctrl_filt.metric==metric]
model = smf.mixedlm(f" auc ~ time", df_used, groups=df_used["subject"])
result = model.fit(reml=False)

coef = result.params.get("time[T.2]")
pval = result.pvalues.get("time[T.2]")
tval = result.tvalues.get("time[T.2]")

df_results = pd.DataFrame({
    'beta':[coef],
    't_value': [tval],
    'p_value': [pval]
})

df_results

#%% Global Metrics Statistical Analysis Type 3: Cosmonauts (Post - Pre) vs Controls (Post - Pre) using Linear Mixed Models

atlas = 'SCH100'
graph = 'pos'

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':   
    R = 400 
elif atlas=='AAL':
    R = 170

df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))
df_ctrl = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Controls_auc_{atlas}_graph-{graph}.xlsx'))
df_cosm['group'] = 'cosmonaut'
df_ctrl['group'] = 'control'
df_ctrl.loc[(df_ctrl['time'] == 1), 'time'] = 'pre2'
df_ctrl.loc[(df_ctrl['time'] == 2), 'time'] = 'post'
df_ctrl['flight'] = 'f1'
df = pd.concat((df_cosm, df_ctrl), ignore_index=True)
df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])

metric = 'modularity' # Possible metrics: 
                             # 'global_efficiency'
                             # 'char_path_length'
                             # 'modularity'


df_used = df[df.metric==metric]
model = smf.mixedlm(f" auc ~ time*group", df_used, groups=df_used["subject"])
result = model.fit(reml=False)

coef_time = result.params.get("time[T.post]", float("nan"))
coef_group = result.params.get("group[T.cosmonaut]", float("nan"))
coef_interaction = result.params.get("time[T.post]:group[T.cosmonaut]", float("nan"))
pval_time = result.pvalues.get("time[T.post]", float("nan"))
pval_group = result.pvalues.get("group[T.cosmonaut]", float("nan"))
pval_interaction = result.pvalues.get("time[T.post]:group[T.cosmonaut]", float("nan"))
tval_time = result.tvalues.get("time[T.post]", float("nan"))
tval_group = result.tvalues.get("group[T.cosmonaut]", float("nan"))
tval_interaction = result.tvalues.get("time[T.post]:group[T.cosmonaut]", float("nan"))

df_results = pd.DataFrame({
    "beta_time": [coef_time],
    "beta_group": [coef_group],
    "beta_interaction": [coef_interaction],
    "t_value_time": [tval_time],
    "t_value_group": [tval_group],
    "t_value_interaction": [tval_interaction],
    "p_value_time": [pval_time],
    "p_value_group": [pval_group],
    "p_value_interaction": [pval_interaction]
})

df_results

#%% Nodal Metrics Statistical Analysis Type 1: Cosmonauts (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'
R = 100
graph = 'pos'

df = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))
df = df[df.flight=='f1']
df = df[(df.time=='pre2') | (df.time=='post')]
df = df.reset_index()
df = df.drop(axis='columns', labels=['index'])
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])

metrics = [
    'clustering',
    'betweenness',
    'closeness'
]

df_final = pd.DataFrame([])
for metric in metrics:
    results = []
    df_metric = df[df.metric==metric]
    for r in range(1, R+1):
        node = f'R{r}'
        df_node = df_metric[df_metric['region'] == node]
        model = smf.mixedlm("auc ~ time", df_node, groups=df_node["subject"])
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

    df_final = pd.concat((df_final, df_results), ignore_index=True)


df_final.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Statistics_Cosmomnauts_nodalAUC_graph-{graph}.xlsx'), index=False)

#%% Nodal Metrics Statistical Analysis Type 2: Controls (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'
R = 100
graph = 'pos'

df = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Controls_auc_{atlas}_graph-{graph}.xlsx'))
df = df[(df.time==1) | (df.time==2)]
df = df.reset_index()
df = df.drop(axis='columns', labels=['index'])
df['time'] = pd.Categorical(df['time'], categories=[1, 2])

metrics = [
    'clustering',
    'betweenness',
    'closeness'
]

df_final = pd.DataFrame([])
for metric in metrics:
    results = []
    df_metric = df[df.metric==metric]
    for r in range(1, R+1):
        node = f'R{r}'
        df_node = df_metric[df_metric['region'] == node]
        model = smf.mixedlm("auc ~ time", df_node, groups=df_node["subject"])
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

    df_final = pd.concat((df_final, df_results), ignore_index=True)

df_final.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Statistics_Controls_nodalAUC_graph-{graph}.xlsx'), index=False)

#%% Nodal Metrics Statistical Analysis Type 3: Cosmonauts (Post - Pre) vs Controls (Post - Pre) using Linear Mixed Models

atlas = 'SCH100'
R = 100
graph = 'pos'

df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))
df_ctrl = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Controls_auc_{atlas}_graph-{graph}.xlsx'))
df_cosm['group'] = 'cosmonaut'
df_ctrl['group'] = 'control'
df_ctrl.loc[(df_ctrl['time'] == 1), 'time'] = 'pre2'
df_ctrl.loc[(df_ctrl['time'] == 2), 'time'] = 'post'
df_ctrl['flight'] = 'f1'
df = pd.concat((df_cosm, df_ctrl), ignore_index=True)
df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])

metrics = [
    'clustering',
    'betweenness',
    'closeness'
]

df_final = pd.DataFrame([])
for metric in metrics:
    results = []
    df_metric = df[df.metric==metric]
    for r in range(1, R+1):
        node = f'R{r}'
        df_node = df_metric[df_metric['region'] == node]
        model = smf.mixedlm("auc ~ time*group", df_node, groups=df_node["subject"])
        result = model.fit(reml=False)

        coef_time = result.params.get("time[T.post]", float("nan"))
        coef_group = result.params.get("group[T.cosmonaut]", float("nan"))
        coef_interaction = result.params.get("time[T.post]:group[T.cosmonaut]", float("nan"))
        pval_time = result.pvalues.get("time[T.post]", float("nan"))
        pval_group = result.pvalues.get("group[T.cosmonaut]", float("nan"))
        pval_interaction = result.pvalues.get("time[T.post]:group[T.cosmonaut]", float("nan"))
        tval_time = result.tvalues.get("time[T.post]", float("nan"))
        tval_group = result.tvalues.get("group[T.cosmonaut]", float("nan"))
        tval_interaction = result.tvalues.get("time[T.post]:group[T.cosmonaut]", float("nan"))

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

    df_final = pd.concat((df_final, df_results), ignore_index=True)

df_final.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Statistics_Cosmonauts_vs_Controls_nodalAUC_graph-{graph}.xlsx'), index=False)
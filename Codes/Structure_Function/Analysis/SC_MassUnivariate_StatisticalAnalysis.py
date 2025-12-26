#%% Importing 
import os 
import os.path as op 
import numpy as np 
import pandas as pd 
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

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

#%% Statistical Analysis Type 1: Cosmonauts (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'
norm = False

df = pd.read_csv(op.join(res_dir, f'{atlas}/SC/MassUnivariate/SC_Cosmonauts_norm-{norm}.csv'))

df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])
df['log_val'] = np.log2(df['val']+1)

results = []
for edge in tqdm(df['edge'].unique(), desc="Running LMMs per edge"):
    df_edge = df[df['edge'] == edge]
    model = smf.mixedlm("log_val ~ time", df_edge, groups=df_edge["subject"])
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

df_results = pd.DataFrame(results)
valid_mask = df_results['p_value'].notna()
df_results['p_fdr'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr'] = corrected_pvals
df_results['significant'] = df_results['p_fdr'] < 0.05


df_results.to_excel(op.join(res_dir, f'{atlas}/SC/MassUnivariate/Statistics_Cosmonauts_LMM_PostPre_norm-{norm}.xlsx'), index=False)

#%% Statistical Analysis Type 2: Controls (Post vs Pre) using Linear Mixed Models

atlas = 'SCH100'
norm = False

df = pd.read_csv(op.join(res_dir, f'{atlas}/SC/MassUnivariate/SC_Controls_norm-{norm}.csv'))

df = df[df['time'].isin([1, 2])].copy()
df['time'] = pd.Categorical(df['time'], categories=[1, 2])
df['log_val'] = np.log2(df['val']+1)

results = []
for edge in tqdm(df['edge'].unique(), desc="Running LMMs per edge"):
    df_edge = df[df['edge'] == edge]
    model = smf.mixedlm("log_val ~ time", df_edge, groups=df_edge["subject"])
    result = model.fit(reml=False)
    coef = result.params.get("time[T.2]", float("nan"))
    pval = result.pvalues.get("time[T.2]", float("nan"))
    tval = result.tvalues.get("time[T.2]", float("nan"))

    results.append({
        "edge": edge,
        "beta": coef,
        "t_value": tval,
        "p_value": pval
    })

df_results = pd.DataFrame(results)
valid_mask = df_results['p_value'].notna()
df_results['p_fdr'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr'] = corrected_pvals
df_results['significant'] = df_results['p_fdr'] < 0.05


df_results.to_excel(op.join(res_dir, f'{atlas}/SC/MassUnivariate/Statistics_Controls_LMM_PostPre_norm-{norm}.xlsx'), index=False)

#%% Statistical Analysis Type 3: Cosmonauts (Post - Pre) vs Controls (Post - Pre) using Linear Mixed Models

atlas = 'SCH100'
norm = False

df_cosm = pd.read_csv(op.join(res_dir, f'{atlas}/SC/MassUnivariate/SC_Cosmonauts_norm-{norm}.csv'))
df_ctrl = pd.read_csv(op.join(res_dir, f'{atlas}/SC/MassUnivariate/SC_Controls_norm-{norm}.csv'))

df_cosm['group'] = 'cosmonaut'
df_ctrl['group'] = 'control'

df_ctrl.loc[(df_ctrl['time'] == 1), 'time'] = 'pre2'
df_ctrl.loc[(df_ctrl['time'] == 2), 'time'] = 'post'
df_ctrl['flight'] = 'f1'

df = pd.concat((df_cosm, df_ctrl), ignore_index=True)

df = df[(df['flight'] == 'f1') & (df['time'].isin(['pre2', 'post']))].copy()
df['time'] = pd.Categorical(df['time'], categories=['pre2', 'post'])
df['log_val'] = np.log2(df['val']+1)

results = []
for edge in tqdm(df['edge'].unique(), desc="Running LMMs per edge"):
    df_edge = df[df['edge'] == edge]
    model = smf.mixedlm("log_val ~ time*group", df_edge, groups=df_edge["subject"])
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
        "edge": edge,
        "beta_time": coef_time,
        "beta_group": coef_group,
        "beta_interaction": coef_interaction,
        "t_value_time": tval_time,
        "t_value_group": tval_group,
        "t_value_interaction": tval_interaction,
        "p_value_time": pval_time,
        "p_value_group": pval_group,
        "p_value_interaction": pval_interaction
    })

df_results = pd.DataFrame(results)
valid_mask = df_results['p_value_time'].notna()
df_results['p_fdr_time'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value_time'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr_time'] = corrected_pvals
df_results['significant_time'] = df_results['p_fdr_time'] < 0.05

valid_mask = df_results['p_value_group'].notna()
df_results['p_fdr_group'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value_group'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr_group'] = corrected_pvals
df_results['significant_group'] = df_results['p_fdr_group'] < 0.05

valid_mask = df_results['p_value_interaction'].notna()
df_results['p_fdr_interaction'] = np.nan
corrected_pvals = multipletests(df_results.loc[valid_mask, 'p_value_interaction'], method='fdr_bh')[1]
df_results.loc[valid_mask, 'p_fdr_interaction'] = corrected_pvals
df_results['significant_interaction'] = df_results['p_fdr_interaction'] < 0.05

df_results.to_excel(op.join(res_dir, f'{atlas}/SC/MassUnivariate/Statistics_Combined_LMM_PostPre_norm-{norm}.xlsx'), index=False)

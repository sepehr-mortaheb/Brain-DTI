#%% Importing

import numpy as np 
import pandas as pd 
import networkx as nx 
import community as community_louvain
import os.path as op
import matplotlib.pylab as plt 
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import auc
import os
from nilearn import plotting, datasets, image
from scipy.ndimage import center_of_mass
import nibabel as nb
from scipy.io import loadmat

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

#%% Graph Global Metrics Estimation 

atlas = 'SCH100' 

df = pd.DataFrame([])
graph = 'absolute' # pos or absolute

thresholds = [
    0.25,
    0.275,
    0.3,
    0.325,
    0.35,
    0.375,
    0.4,
    0.425,
    0.45,
    0.475,
    0.5, 
    0.525, 
    0.55, 
    0.575, 
    0.6, 
    0.625, 
    0.65, 
    0.675, 
    0.7, 
    0.725, 
    0.75, 
    0.775, 
    0.8, 
]

for threshold in thresholds:
    print(f'Cosmonauts: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-cosmonaut'):           
            for ses in sessions[sub]:
                flight = ses.split('-')[1][0:2]
                tp = ses.split('-')[1][2:]

                adj_matrix = loadmat(
                    op.join(data_dir, sub, ses, f'FC_{atlas}.mat')
                )['FC']
                np.fill_diagonal(adj_matrix, 0)

                if graph=='pos':
                    adj_matrix[adj_matrix < 0] = 0
                elif graph=='absolute':
                    adj_matrix = np.abs(adj_matrix)

                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                # Global metrics
                dftmp = pd.DataFrame({})
                dftmp['subject'] = [sub]
                dftmp['flight'] = [flight]
                dftmp['time'] = [tp]

                dftmp['threshold'] = [threshold]

                dftmp["global_efficiency"] = [nx.global_efficiency(G)]
                try:
                    dftmp["char_path_length"] = [nx.average_shortest_path_length(G, weight="weight")]
                except:
                    dftmp["char_path_length"] = [np.nan]
                dftmp["density"] = [nx.density(G)]
                dftmp["transitivity"] = [nx.transitivity(G)]
                partition = community_louvain.best_partition(G, weight='weight')
                dftmp["modularity"] = [community_louvain.modularity(partition, G, weight="weight")]
                #dftmp['small_world'] = [nx.sigma(G)]

                df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_{atlas}_graph-{graph}.xlsx'))

df = pd.DataFrame([])
for threshold in thresholds:
    print(f'Controls: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-control'):
            for ses in sessions[sub]:
                tp = ses.split('-')[1]

                adj_matrix = loadmat(
                    op.join(data_dir, sub, ses, f'FC_{atlas}.mat')
                )['FC']
                np.fill_diagonal(adj_matrix, 0)

                if graph=='pos':
                    adj_matrix[adj_matrix < 0] = 0
                elif graph=='absolute':
                    adj_matrix = np.abs(adj_matrix)

                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                # Global metrics
                dftmp = pd.DataFrame({})
                dftmp['subject'] = [sub]
                dftmp['time'] = [tp]

                dftmp['threshold'] = [threshold]

                dftmp["global_efficiency"] = [nx.global_efficiency(G)]
                try:
                    dftmp["char_path_length"] = [nx.average_shortest_path_length(G, weight="weight")]
                except:
                    dftmp["char_path_length"] = [np.nan]               
                dftmp["density"] = [nx.density(G)]
                dftmp["transitivity"] = [nx.transitivity(G)]
                partition = community_louvain.best_partition(G, weight='weight')
                dftmp["modularity"] = [community_louvain.modularity(partition, G, weight="weight")]
                #dftmp['small_world'] = [nx.sigma(G)]

                df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Controls_{atlas}_graph-{graph}.xlsx'))


# %% Visualization of Global Metrics 

atlas = 'SCH100'
graph = 'absolute'

# For Cosmonauts 
df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_{atlas}_graph-{graph}.xlsx'))

df_cosm_filt = df_cosm
df_cosm_filt = df_cosm_filt[df_cosm_filt.flight=='f1']
df_cosm_filt = df_cosm_filt[(df_cosm_filt.time=='pre2') | (df_cosm_filt.time=='post')]
df_cosm_filt = df_cosm_filt.reset_index()
df_cosm_filt = df_cosm_filt.drop(axis='columns', labels=['index'])

metrics = [
    'global_efficiency',
    'char_path_length',
    'density',
    'transitivity',
    'modularity'
]

df_cosm_filt_auc = pd.DataFrame([])

for subject in np.unique(df_cosm_filt.subject):
    for time in np.unique(df_cosm_filt.time):
        for metric in metrics:
            x = df_cosm_filt[(df_cosm_filt.time == time) & (df_cosm_filt.subject == subject)]['threshold']
            y = df_cosm_filt[(df_cosm_filt.time == time) & (df_cosm_filt.subject == subject)][metric]
            if len(x) > 1:
                auc_val = auc(x, y)
                dftmp = pd.DataFrame([])
                dftmp['subject'] = [subject]
                dftmp['time'] = [time]
                dftmp['metric'] = [metric]
                dftmp['auc'] = [auc_val]
                df_cosm_filt_auc = pd.concat((df_cosm_filt_auc, dftmp), ignore_index=True)
            
df_cosm_filt_auc.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))

for metric in metrics:
    fig, ax = plt.subplots(1, 2, figsize=(17, 5))

    sns.pointplot(
        x='threshold', 
        y=metric,
        data=df_cosm_filt, 
        hue='time',
        hue_order=['pre2', 'post'],
        ax=ax[0],
    )

    thresholds = [
        0.25,
        0.275,
        0.3,
        0.325,
        0.35,
        0.375,
        0.4,
        0.425,
        0.45,
        0.475,
        0.5, 
        0.525, 
        0.55, 
        0.575, 
        0.6, 
        0.625, 
        0.65, 
        0.675, 
        0.7, 
        0.725, 
        0.75, 
        0.775, 
        0.8, 
    ]
    ax[0].set_xticklabels(thresholds, rotation=45)
    ax[0].set_xlabel('Threshold', size=15)
    ax[0].set_ylabel(metric, size=12)
    ax[0].grid(True)

    sns.violinplot(
        x='time',
        y='auc',
        data=df_cosm_filt_auc[df_cosm_filt_auc.metric==metric],
        order=['pre2', 'post'],
        ax=ax[1],
        color='white'
    )

    sns.stripplot(
        x='time',
        y='auc',
        data=df_cosm_filt_auc[df_cosm_filt_auc.metric==metric],
        order=['pre2', 'post'],
        ax=ax[1], 
        hue='time',
        hue_order=['pre2', 'post'],
        jitter=0
    )

    sns.lineplot(
        x='time',
        y='auc',
        data=df_cosm_filt_auc[df_cosm_filt_auc.metric==metric],
        ax=ax[1],
        hue='subject',
        legend=False,
        alpha=0.4
    )

    ax[1].grid(True)
    ax[1].set_title('Area Under the Curve', size=12)
    ax[1].set_xlabel('', size=15)
    ax[1].set_ylabel('AUC', size=12)

    fig.suptitle(metric, size=15)

    fig.savefig(
        op.join(res_dir, f'{atlas}/FC/Graph/Graph_{metric}_Cosmonauts_{atlas}_graph-{graph}.pdf'),
        dpi=300,
        bbox_inches='tight'
    )
    fig.savefig(
        op.join(res_dir, f'{atlas}/FC/Graph/Graph_{metric}_Cosmonauts_{atlas}_graph-{graph}.png'),
        dpi=300,
        bbox_inches='tight'
    )

#%% Statistics of Global Metrics 

# Statistics per threshold value 
atlas = 'SCH100'
graph = 'absolute'

df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_{atlas}_graph-{graph}.xlsx'))

df_cosm_filt = df_cosm
df_cosm_filt = df_cosm_filt[df_cosm_filt.flight=='f1']
df_cosm_filt = df_cosm_filt[(df_cosm_filt.time=='pre2') | (df_cosm_filt.time=='post')]
df_cosm_filt = df_cosm_filt.reset_index()
df_cosm_filt = df_cosm_filt.drop(axis='columns', labels=['index'])
df_cosm_filt['time'] = pd.Categorical(df_cosm_filt['time'], categories=['pre2', 'post'])

metric = 'global_efficiency' # Possible metrics: 
                             # 'global_efficiency'
                             # 'char_path_length'
                             # 'density'
                             # 'transitivity'
                             # 'modularity

thresholds = [
    0.25,
    0.275,
    0.3,
    0.325,
    0.35,
    0.375,
    0.4,
    0.425,
    0.45,
    0.475,
    0.5, 
    0.525, 
    0.55, 
    0.575, 
    0.6, 
    0.625, 
    0.65, 
    0.675, 
    0.7, 
    0.725, 
    0.75, 
    0.775, 
    0.8, 
]
results = []
for th in thresholds:
    df_th = df_cosm_filt[df_cosm_filt['threshold'] == th]
    model = smf.mixedlm(f" {metric} ~ time", df_th, groups=df_th["subject"])
    result = model.fit(reml=False)

    coef = result.params.get("time[T.post]")
    pval = result.pvalues.get("time[T.post]")
    tval = result.tvalues.get("time[T.post]")

    results.append({
        "threshold": th,
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

df_results

#%% statistics on the AUC value
#   
atlas = 'SCH100'
graph = 'absolute'

df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_global_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))

df_cosm_filt = df_cosm
df_cosm_filt = df_cosm_filt[(df_cosm_filt.time=='pre2') | (df_cosm_filt.time=='post')]
df_cosm_filt = df_cosm_filt.reset_index()
df_cosm_filt = df_cosm_filt.drop(axis='columns', labels=['index'])
df_cosm_filt['time'] = pd.Categorical(df_cosm_filt['time'], categories=['pre2', 'post'])

metric = 'modularity' # Possible metrics: 
                             # 'global_efficiency'
                             # 'char_path_length'
                             # 'density'
                             # 'transitivity'
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

#%% Nodal Metrics 

df = pd.DataFrame([])
graph='absolute'

thresholds = [
    0.25,
    0.275,
    0.3,
    0.325,
    0.35,
    0.375,
    0.4,
    0.425,
    0.45,
    0.475,
    0.5, 
    0.525, 
    0.55, 
    0.575, 
    0.6, 
    0.625, 
    0.65, 
    0.675, 
    0.7, 
    0.725, 
    0.75, 
    0.775, 
    0.8, 
]
for threshold in thresholds:
    print(f'Cosmonauts: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-cosmonaut'):           
            for ses in sessions[sub]:
                flight = ses.split('-')[1][0:2]
                tp = ses.split('-')[1][2:]

                adj_matrix = loadmat(
                    op.join(data_dir, sub, ses, f'FC_{atlas}.mat')
                )['FC']
                np.fill_diagonal(adj_matrix, 0)

                if graph=='pos':
                    adj_matrix[adj_matrix < 0] = 0
                elif graph=='absolute':
                    adj_matrix = np.abs(adj_matrix)

                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                clustering = nx.clustering(G, weight="weight")
                between = nx.betweenness_centrality(G, weight="weight")
                close = nx.closeness_centrality(G, distance="weight")               
                
                for r in range(G.number_of_nodes()):
                    dftmp = pd.DataFrame({})
                    dftmp['subject'] = [sub]
                    dftmp['time'] = [tp]
                    dftmp['flight'] = [flight]
                    dftmp['region'] = [f'R{r+1}']

                    dftmp['threshold'] = [threshold]

                    dftmp['clustering'] = [list(clustering.values())[r]]
                    dftmp['betweenness'] = [list(between.values())[r]]
                    dftmp['closeness'] = [list(close.values())[r]]

                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Cosmonauts_{atlas}_graph-{graph}.xlsx'))

df = pd.DataFrame([])
for threshold in thresholds:
    print(f'Controls: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-control'):
            for ses in sessions[sub]:
                tp = ses.split('-')[1]

                adj_matrix = loadmat(
                    op.join(data_dir, sub, ses, f'FC_{atlas}.mat')
                )['FC']
                np.fill_diagonal(adj_matrix, 0)

                if graph=='pos':
                    adj_matrix[adj_matrix < 0] = 0
                elif graph=='absolute':
                    adj_matrix = np.abs(adj_matrix)

                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                # Node-level metrics (summarized by mean)
                clustering = nx.clustering(G, weight="weight")
                between = nx.betweenness_centrality(G, weight="weight")
                close = nx.closeness_centrality(G, distance="weight")              
                
                for r in range(G.number_of_nodes()):
                    dftmp = pd.DataFrame({})
                    dftmp['subject'] = [sub]
                    dftmp['time'] = [tp]
                    dftmp['region'] = [f'R{r+1}']

                    dftmp['threshold'] = [threshold]

                    dftmp['clustering'] = [list(clustering.values())[r]]
                    dftmp['betweenness'] = [list(between.values())[r]]
                    dftmp['closeness'] = [list(close.values())[r]]

                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Controls_{atlas}_graph-{graph}.xlsx'))

#%% Calculation of Nodal Metrics AUC 

atlas = 'SCH100'
graph = 'absolute'

# For Cosmonauts 
df_cosm = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Cosmonauts_{atlas}_graph-{graph}.xlsx'))

df_cosm_filt = df_cosm
df_cosm_filt = df_cosm_filt[df_cosm_filt.flight=='f1']
df_cosm_filt = df_cosm_filt[(df_cosm_filt.time=='pre2') | (df_cosm_filt.time=='post')]
df_cosm_filt = df_cosm_filt.reset_index()
df_cosm_filt = df_cosm_filt.drop(axis='columns', labels=['index'])

metrics = [
    'clustering',
    'betweenness',
    'closeness'
]

df_cosm_filt_auc = pd.DataFrame([])

for subject in np.unique(df_cosm_filt.subject):
    print(subject)
    for reg in np.unique(df_cosm_filt.region):
        for time in np.unique(df_cosm_filt.time):
            for metric in metrics:
                x = df_cosm_filt[(df_cosm_filt.time == time) & (df_cosm_filt.subject == subject) & (df_cosm_filt.region == reg)]['threshold']
                y = df_cosm_filt[(df_cosm_filt.time == time) & (df_cosm_filt.subject == subject) & (df_cosm_filt.region == reg)][metric]
                if len(x) > 1:
                    auc_val = auc(x, y)
                    dftmp = pd.DataFrame([])
                    dftmp['subject'] = [subject]
                    dftmp['time'] = [time]
                    dftmp['region'] = [reg]
                    dftmp['metric'] = [metric]
                    dftmp['auc'] = [auc_val]
                    df_cosm_filt_auc = pd.concat((df_cosm_filt_auc, dftmp), ignore_index=True)
                

df_cosm_filt_auc.to_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))


#%% Statistics on the nodal metrics AUC 

atlas = 'SCH100'
R = 100
graph = 'absolute'

# For Cosmonauts 
df = pd.read_excel(op.join(res_dir, f'{atlas}/FC/Graph/Graph_metrics_nodal_Cosmonauts_auc_{atlas}_graph-{graph}.xlsx'))
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

#%% Visualization of the results 

atlas = 'SCH100'
graph = 'absolute'
metric = 'closeness'

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
    f'{atlas}/FC/Graph/Statistics_Cosmomnauts_nodalAUC_graph-{graph}.xlsx'
)

df = pd.read_excel(file_addr)
df_metric = df[df.metric==metric]

# Filter significant nodes
sig_nodes = df_metric[df_metric["p_fdr"] < 0.05].copy()

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
            f'{atlas}/FC/Graph/Statistics_Cosmomnauts_nodalAUC_{metric}_graph-{graph}.nii'
        )
    )

    print(sig_nodes)
#%% Importing

import numpy as np 
import pandas as pd 
import networkx as nx 
import community as community_louvain
import os.path as op
import matplotlib.pylab as plt 
import seaborn as sns

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

binary = False # Bool, whether to binarize the matrix 

#%% Graph parameters Estimation 
df = pd.DataFrame([])

for threshold in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]:
    print(f'Cosmonauts: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-cosmonaut'):           
            for ses in sessions[sub]:
                flight = ses.split('-')[1][0:2]
                tp = ses.split('-')[1][2:]

                adj_matrix = np.loadtxt(
                    op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                    delimiter=','
                )
                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Binarize (optional)
                if binary:
                    adj_matrix = (adj_matrix > 0).astype(int)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                # Global metrics
                dftmp = pd.DataFrame({})
                dftmp['subject'] = [sub]
                dftmp['flight'] = [flight]
                dftmp['time'] = [tp]
                dftmp["global_efficiency"] = [nx.global_efficiency(G)]
                dftmp["average_clustering"] = [nx.average_clustering(G, weight=None if binary else "weight")]
                dftmp["density"] = [nx.density(G)]
                dftmp["transitivity"] = [nx.transitivity(G)]

                # Shortest path metrics
                try:
                    dftmp["char_path_length"] = [nx.average_shortest_path_length(G, weight=None if binary else "weight")]
                except:
                    dftmp["char_path_length"] = [np.nan]

                # Modularity (Louvain)
                partition = community_louvain.best_partition(G, weight=None if binary else "weight")
                dftmp["modularity"] = [community_louvain.modularity(partition, G, weight=None if binary else "weight")]

                # Node-level metrics (summarized by mean)
                deg = dict(G.degree(weight=None if binary else None))
                strength = dict(G.degree(weight="weight")) if not binary else deg
                between = nx.betweenness_centrality(G, weight=None if binary else "weight")
                close = nx.closeness_centrality(G, distance=None if binary else "weight")
                eigen = nx.eigenvector_centrality_numpy(G, weight=None if binary else "weight")

                # Add average node-level metrics
                dftmp["mean_degree"] = [np.mean(list(deg.values()))]
                dftmp["mean_strength"] = [np.mean(list(strength.values()))]
                dftmp["mean_betweenness"] = [np.mean(list(between.values()))]
                dftmp["mean_closeness"] = [np.mean(list(close.values()))]
                dftmp["mean_eigenvector"] = [np.mean(list(eigen.values()))]
                dftmp['threshold'] = [threshold]

                df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'SC/Graph/Graph_metrics_Cosmonauts_{atlas}_bin-{binary}.xlsx'))

df = pd.DataFrame([])
for threshold in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]:
    print(f'Controls: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-control'):
            for ses in sessions[sub]:
                tp = ses.split('-')[1]
                adj_matrix = np.loadtxt(
                    op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                    delimiter=','
                )
                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Binarize (optional)
                if binary:
                    adj_matrix = (adj_matrix > 0).astype(int)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                # Global metrics
                dftmp = pd.DataFrame({})
                dftmp['subject'] = [sub]
                dftmp['time'] = [tp]
                dftmp["global_efficiency"] = [nx.global_efficiency(G)]
                dftmp["average_clustering"] = [nx.average_clustering(G, weight=None if binary else "weight")]
                dftmp["density"] = [nx.density(G)]
                dftmp["transitivity"] = [nx.transitivity(G)]

                # Shortest path metrics
                try:
                    dftmp["char_path_length"] = [nx.average_shortest_path_length(G, weight=None if binary else "weight")]
                except:
                    dftmp["char_path_length"] = [np.nan]

                # Modularity (Louvain)
                partition = community_louvain.best_partition(G, weight=None if binary else "weight")
                dftmp["modularity"] = [community_louvain.modularity(partition, G, weight=None if binary else "weight")]

                # Node-level metrics (summarized by mean)
                deg = dict(G.degree(weight=None if binary else None))
                strength = dict(G.degree(weight="weight")) if not binary else deg
                between = nx.betweenness_centrality(G, weight=None if binary else "weight")
                close = nx.closeness_centrality(G, distance=None if binary else "weight")
                eigen = nx.eigenvector_centrality_numpy(G, weight=None if binary else "weight")

                # Add average node-level metrics
                dftmp["mean_degree"] = [np.mean(list(deg.values()))]
                dftmp["mean_strength"] = [np.mean(list(strength.values()))]
                dftmp["mean_betweenness"] = [np.mean(list(between.values()))]
                dftmp["mean_closeness"] = [np.mean(list(close.values()))]
                dftmp["mean_eigenvector"] = [np.mean(list(eigen.values()))]
                dftmp['threshold'] = [threshold]

                df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'SC/Graph/Graph_metrics_Controls_{atlas}_bin-{binary}.xlsx'))

# %% Visualization 

atlas = 'SCH100'
binary = False 

# For Cosmonauts 
df_cosm = pd.read_excel(op.join(res_dir, f'SC/Graph/Graph_metrics_Cosmonauts_{atlas}_bin-{binary}.xlsx'))

exc_list = [
  "sub-cosmonaut13", 
  "sub-cosmonaut15", 
  "sub-cosmonaut16",
  "sub-cosmonaut26"
] 

df_cosm_filt = df_cosm
for sub in exc_list:
    df_cosm_filt = df_cosm_filt[df_cosm_filt.subject!=sub]
df_cosm_filt = df_cosm_filt.reset_index()
df_cosm_filt = df_cosm_filt.drop('index', axis=1)

df_cosm_filt = df_cosm_filt[df_cosm_filt.flight=='f1']
df_cosm_filt = df_cosm_filt[(df_cosm_filt.time=='pre2') | (df_cosm_filt.time=='post')]


metrics = [
    'global_efficiency',
    'average_clustering',
    'density',
    'transitivity',
    'char_path_length',
    'modularity',
    'mean_degree',
    'mean_strength',
    'mean_betweenness',
    'mean_closeness',
    'mean_eigenvector',
]

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

    ax[0].grid(True)
    fig.suptitle(metric, size=15)

#%%
df_cosm_filt
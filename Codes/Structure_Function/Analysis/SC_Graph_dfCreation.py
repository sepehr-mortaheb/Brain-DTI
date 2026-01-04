#%% Importing

import os
import os.path as op
import numpy as np 
import pandas as pd 
import networkx as nx 
import community as community_louvain
from sklearn.metrics import auc

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

#%% Graph Global Metrics Estimation 

atlas = 'SCH100' # SCH100, SCH400, AAL
norm = False

if atlas=='SCH100':
    R = 100 
elif atlas=='SCH400':
    R = 400 
elif atlas=='AAL':
    R = 170

# For Cosmonauts 

df = pd.DataFrame([])
for threshold in [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]:
    print(f'Cosmonauts: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-cosmonaut'):           
            for ses in sessions[sub]:
                flight = ses.split('-')[1][0:2]
                tp = ses.split('-')[1][2:]
                if norm:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                        delimiter=','
                    )
                else:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}_nonorm.csv'), 
                        delimiter=','
                    )
                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                largest_cc = max(nx.connected_components(G), key=len)
                G = G.subgraph(largest_cc).copy()

                # Global metrics
                dftmp = pd.DataFrame({})
                dftmp['subject'] = [sub]
                dftmp['flight'] = [flight]
                dftmp['time'] = [tp]

                dftmp['threshold'] = [threshold]

                dftmp["global_efficiency"] = [nx.global_efficiency(G)]
                dftmp["char_path_length"] = [nx.average_shortest_path_length(G, weight="weight")]
                dftmp["density"] = [nx.density(G)]
                dftmp["transitivity"] = [nx.transitivity(G)]
                partition = community_louvain.best_partition(G, weight='weight')
                dftmp["modularity"] = [community_louvain.modularity(partition, G, weight="weight")]

                df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_global_Cosmonauts_{atlas}_norm-{norm}.xlsx'))

# AUC calculation for Cosmonauts

metrics = [
    'global_efficiency',
    'char_path_length',
    'density',
    'transitivity',
    'modularity'
]

df_auc = pd.DataFrame([])

for subject in np.unique(df.subject):
    for flight in np.unique(df.flight):
        for time in np.unique(df.time):
            for metric in metrics:
                x = df[(df.flight == flight) & (df.time == time) & (df.subject == subject)]['threshold']
                y = df[(df.flight == flight) & (df.time == time) & (df.subject == subject)][metric]
                if len(x) > 1:
                    auc_val = auc(x, y)
                    dftmp = pd.DataFrame([])
                    dftmp['subject'] = [subject]
                    dftmp['flight'] = [flight]
                    dftmp['time'] = [time]
                    dftmp['metric'] = [metric]
                    dftmp['auc'] = [auc_val]
                    df_auc = pd.concat((df_auc, dftmp), ignore_index=True)
                

df_auc.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_global_Cosmonauts_auc_{atlas}_norm-{norm}.xlsx'))

# For Controls 

df = pd.DataFrame([])
for threshold in [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]:
    print(f'Controls: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-control'):
            for ses in sessions[sub]:
                tp = ses.split('-')[1]
                if norm:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                        delimiter=','
                    )
                else:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}_nonorm.csv'), 
                        delimiter=','
                    )
                # Keep top X% of weights
                flat = adj_matrix[np.triu_indices_from(adj_matrix, k=1)]
                thresh_val = np.percentile(flat, 100 - threshold * 100)
                adj_matrix = np.where(adj_matrix >= thresh_val, adj_matrix, 0)

                # Create undirected graph
                G = nx.from_numpy_array(adj_matrix)

                # Remove isolated nodes
                G.remove_nodes_from(list(nx.isolates(G)))

                largest_cc = max(nx.connected_components(G), key=len)
                G = G.subgraph(largest_cc).copy()

                # Global metrics
                dftmp = pd.DataFrame({})
                dftmp['subject'] = [sub]
                dftmp['time'] = [tp]

                dftmp['threshold'] = [threshold]

                dftmp["global_efficiency"] = [nx.global_efficiency(G)]
                dftmp["char_path_length"] = [nx.average_shortest_path_length(G, weight="weight")]
                dftmp["density"] = [nx.density(G)]
                dftmp["transitivity"] = [nx.transitivity(G)]
                partition = community_louvain.best_partition(G, weight='weight')
                dftmp["modularity"] = [community_louvain.modularity(partition, G, weight="weight")]

                df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_global_Controls_{atlas}_norm-{norm}.xlsx'))

# AUC calculation for Controls

metrics = [
    'global_efficiency',
    'char_path_length',
    'density',
    'transitivity',
    'modularity'
]

df_auc = pd.DataFrame([])

for subject in np.unique(df.subject):
    for time in np.unique(df.time):
        for metric in metrics:
            x = df[(df.time == time) & (df.subject == subject)]['threshold']
            y = df[(df.time == time) & (df.subject == subject)][metric]
            if len(x) > 1:
                auc_val = auc(x, y)
                dftmp = pd.DataFrame([])
                dftmp['subject'] = [subject]
                dftmp['time'] = [time]
                dftmp['metric'] = [metric]
                dftmp['auc'] = [auc_val]
                df_auc = pd.concat((df_auc, dftmp), ignore_index=True)
            

df_auc.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_global_Controls_auc_{atlas}_norm-{norm}.xlsx'))

#%% Graph Nodal Metrics Estimation 

atlas = 'SCH100'
norm = False

# For Cosmonauts

df = pd.DataFrame([])

for threshold in [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]:
    print(f'Cosmonauts: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-cosmonaut'):           
            for ses in sessions[sub]:
                flight = ses.split('-')[1][0:2]
                tp = ses.split('-')[1][2:]

                if norm:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                        delimiter=','
                    )
                else:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}_nonorm.csv'), 
                        delimiter=','
                    )
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
                    dftmp['flight'] = [flight]
                    dftmp['region'] = [f'R{r+1}']

                    dftmp['threshold'] = [threshold]

                    dftmp['clustering'] = [list(clustering.values())[r]]
                    dftmp['betweenness'] = [list(between.values())[r]]
                    dftmp['closeness'] = [list(close.values())[r]]

                    df = pd.concat((df, dftmp), ignore_index=True)

df.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_nodal_Cosmonauts_{atlas}_norm-{norm}.xlsx'))

# AUC calculation for Cosmonauts
metrics = [
    'clustering',
    'betweenness',
    'closeness'
]

df_auc = pd.DataFrame([])

for subject in np.unique(df.subject):
    for reg in np.unique(df.region):
        for flight in np.unique(df.flight):
            for time in np.unique(df.time):
                for metric in metrics:
                    x = df[(df.flight == flight) & (df.time == time) & (df.subject == subject) & (df.region == reg)]['threshold']
                    y = df[(df.flight == flight) & (df.time == time) & (df.subject == subject) & (df.region == reg)][metric]
                    if len(x) > 1:
                        auc_val = auc(x, y)
                        dftmp = pd.DataFrame([])
                        dftmp['subject'] = [subject]
                        dftmp['time'] = [time]
                        dftmp['flight'] = [flight]
                        dftmp['region'] = [reg]
                        dftmp['metric'] = [metric]
                        dftmp['auc'] = [auc_val]
                        df_auc = pd.concat((df_auc, dftmp), ignore_index=True)          

df_auc.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_nodal_Cosmonauts_auc_{atlas}_norm-{norm}.xlsx'))

#%%
# For Controls

df = pd.DataFrame([])
for threshold in [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]:
    print(f'Controls: th = {threshold}')
    for sub in subjects:
        if sub.startswith('sub-control'):
            for ses in sessions[sub]:
                tp = ses.split('-')[1]

                if norm:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}.csv'), 
                        delimiter=','
                    )
                else:
                    adj_matrix = np.loadtxt(
                        op.join(data_dir, sub, ses, f'SC_{atlas}_nonorm.csv'), 
                        delimiter=','
                    )

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

df.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_nodal_Controls_{atlas}_norm-{norm}.xlsx'))

# AUC calculation for Controls
metrics = [
    'clustering',
    'betweenness',
    'closeness'
]

df_auc = pd.DataFrame([])

for subject in np.unique(df.subject):
    for reg in np.unique(df.region):
        for time in np.unique(df.time):
            for metric in metrics:
                x = df[(df.time == time) & (df.subject == subject) & (df.region == reg)]['threshold']
                y = df[(df.time == time) & (df.subject == subject) & (df.region == reg)][metric]
                if len(x) > 1:
                    auc_val = auc(x, y)
                    dftmp = pd.DataFrame([])
                    dftmp['subject'] = [subject]
                    dftmp['time'] = [time]
                    dftmp['region'] = [reg]
                    dftmp['metric'] = [metric]
                    dftmp['auc'] = [auc_val]
                    df_auc = pd.concat((df_auc, dftmp), ignore_index=True)          

df_auc.to_excel(op.join(res_dir, f'{atlas}/SC/Graph/Graph_metrics_nodal_Controls_auc_{atlas}_norm-{norm}.xlsx'))

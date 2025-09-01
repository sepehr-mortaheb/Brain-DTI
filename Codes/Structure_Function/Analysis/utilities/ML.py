import numpy as np
from scipy.io import loadmat  
from sklearn.cluster import KMeans
#import alphashape

######################################################################

def K_means(X, n_clusters=4):
    data = np.transpose(X)
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(data)
    return np.transpose(kmeans.cluster_centers_), kmeans.labels_

def convert_conn(conn, nec_files_addr):
    ind_file = loadmat(nec_files_addr+'match_Yeo.mat')
    indT = np.zeros(360)
    length = []
    c=0
    for net in ['DA', 'DMN', 'FP', 'LIM', 'SM', 'VA', 'VIS']:
        ind = ind_file['ind_'+net]
        ind = ind - 1 
        length.append(int(len(ind)))
        indT[c:c+len(ind)] = ind.reshape(len(ind))
        c = c+len(ind)
    ind = []
    for i in range(len(indT)):
        ind.append(int(indT[i]))
    conn_net = conn[ind,:][:, ind]
    return conn_net, length


import pygsp 
from pygsp import graphs
import numpy as np
import matplotlib.pylab as plt 
import seaborn as sns 
from scipy.io import loadmat 
from os import listdir 
from sklearn.cluster import KMeans
#import alphashape
#from descartes import PolygonPatch
import sys
from numpy import linalg as LA

def get_subjects_id():
    drct = '../../Data/Structural'
    allfiles = listdir(drct)
    files = [f for f in allfiles if not f.startswith('.')]
    subj_id = [f.split('_')[-1].split('.')[0] for f in files]
    return subj_id

def read_SC(subj_id):
    sc_file = 'Connectome_sb_'+subj_id+'.mat'
    sc = loadmat('../../Data/Structural/'+sc_file)['Cnorm']
    return sc 

def read_rs(subj_id):
    rs_file = 'Rest1LR_Sub'+subj_id+'_Glasser.mat'
    rs = loadmat('../../Data/Functional/Rest/Rest1LR/'+rs_file)['TCS']
    return rs 

def Graph_creation(mat, kind): 
    G = graphs.Graph(mat)
    if kind == 'combinatorial':
        G.compute_laplacian()
        G.compute_fourier_basis(recompute=True)
        eigenValues = G.e
        eigenVectors = G.U
    elif kind == 'normalized':
        G.compute_laplacian('normalized') 
        G.compute_fourier_basis(recompute=True)
        eigenValues = G.e
        eigenVectors = G.U
    elif kind == 'adjacency':
        eigenValues , eigenVectors = LA.eig(mat)
        idx = eigenValues.argsort()   
        eigenValues = eigenValues[idx]
        eigenVectors = eigenVectors[:,idx]
    return eigenValues, eigenVectors 

def LPF(X, eigval, eigvec, low, kind): 
    G_low = _LPF_matrix(len(eigval), low, kind)
    U = eigvec
    X_low = np.matmul(np.matmul(np.matmul(U, G_low), np.transpose(U)), X)
    return X_low

def HPF(X, eigval, eigvec, high, kind): 
    G_high = _HPF_matrix(len(eigval), high, kind)
    U = eigvec
    X_high = np.matmul(np.matmul(np.matmul(U, G_high), np.transpose(U)), X)
    return X_high

def Alignment(X, eigval, eigvec, low, kind): 
    X_low = LPF(X, eigval, eigvec, low, kind)
    return LA.norm(X_low, axis=0)

def Liberality(X, eigval, eigvec, high, kind): 
    X_high = HPF(X, eigval, eigvec, high, kind)
    return LA.norm(X_high, axis=0)

def SDI_calculation(X, eigval, eigvec):
    U = eigvec
    X_hat = np.matmul(np.transpose(U), X)
    cut_lambda = _freq_selector(X_hat, eigval, eigvec)
    Xc = LPF(X, eigval, eigvec, cut_lambda)
    Xd = HPF(X, eigval, eigvec, len(eigval)-cut_lambda)
    NXc = LA.norm(Xc, axis=1)
    NXd = LA.norm(Xd, axis=1)
    return NXd/NXc, cut_lambda

def Group_SDI_calculation(subjects):
    cut_lambda = _Group_freq_selector(subjects)
    sc = _Group_avg_SC(subjects)
    eigval, eigvec = Graph_creation(sc, kind='normalized')
    U = eigvec
    #TSDI = []
    NXC = []
    NXD = []
    for subj in subjects:
        rest = read_rs(subj)
        Xc = LPF(rest, eigval, eigvec, cut_lambda)
        Xd = HPF(rest, eigval, eigvec, len(eigval)-cut_lambda)
        NXc = np.linalg.norm(Xc, axis=1)
        NXd = np.linalg.norm(Xd, axis=1)
        NXC.append(NXc)
        NXD.append(NXd)
        #SDI = NXd/NXc
        #TSDI.append(SDI)
    return np.mean(NXD, 0)/np.mean(NXC, 0)

def dSDI_calculation(X, eigval, eigvec): 
    U = eigvec
    X_hat = np.matmul(np.transpose(U), X)
    cut_lambda = _freq_selector(X_hat, eigval, eigvec)
    Xc = LPF(X, eigval, eigvec, cut_lambda)
    Xd = HPF(X, eigval, eigvec, len(eigval)-cut_lambda)
    return Xd/Xc, cut_lambda

def K_means(X, n_clusters=4):
    data = np.transpose(X)
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(data)
    return np.transpose(kmeans.cluster_centers_), kmeans.labels_

def convert_conn(conn):
    ind_file = loadmat('match_Yeo.mat')
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

def plot_conn_patterns(conn_patt, length):
    if len(conn_patt.shape)>2:
        n_patterns = conn_patt.shape[2]
    else:
        n_patterns = 1
    cumlength = np.cumsum(length)
    cumlength = cumlength[0:-1]
    fig, ax = plt.subplots(1,n_patterns, figsize=(6*n_patterns,6))
    if n_patterns>1:
        for patt in range(n_patterns):
            ax[patt].imshow(conn_patt[:,:,patt])
            for i in range(len(cumlength)):
                ax[patt].axvline(cumlength[i], color='red')
                ax[patt].axhline(cumlength[i], color='red')
                ax[patt].set_xticks([])
                ax[patt].set_yticks([])
            loc = np.cumsum(length)
            ax[patt].text(0, -5, 'DA')
            ax[patt].text(loc[0], -5, 'DMN')
            ax[patt].text(loc[1], -5, 'FP')
            ax[patt].text(loc[2], -5, 'LIM')
            ax[patt].text(loc[3], -5, 'SM')
            ax[patt].text(loc[4], -5, 'VA')
            ax[patt].text(loc[5], -5, 'VIS')
    else:
        ax.imshow(conn_patt)
        for i in range(len(cumlength)):
            ax.axvline(cumlength[i], color='red')
            ax.axhline(cumlength[i], color='red')
            ax.set_xticks([])
            ax.set_yticks([])
        loc = np.cumsum(length)
        ax.text(0, -5, 'DA')
        ax.text(loc[0], -5, 'DMN')
        ax.text(loc[1], -5, 'FP')
        ax.text(loc[2], -5, 'LIM')
        ax.text(loc[3], -5, 'SM')
        ax.text(loc[4], -5, 'VA')
        ax.text(loc[5], -5, 'VIS')
    return ax 

def boundary(data, alpha):
    points =[]
    for i in range(len(data)):
        x = (data[i,0], data[i,1])
        points.append(x)
    alpha_shape = alphashape.alphashape(points, alpha)
    return alpha_shape, alpha_shape.area

#def traj_length

###############################################################################

def _LPF_matrix(N,low, kind): 
    if kind in ['normalized', 'combinatorial']:
        G_low = np.zeros([N,N])
        for i in range(low):
            G_low[i,i] = 1
    elif kind == 'adjacency':
        G_low = np.zeros([N,N])
        for i in range(N-low,N):
            G_low[i,i] = 1
    return G_low

def _HPF_matrix(N,high, kind):
    if kind in ['normalized', 'combinatorial']:
        G_high = np.zeros([N,N])
        for i in range(N-high,N):
            G_high[i,i] = 1
    elif kind == 'adjacency':
        G_high = np.zeros([N,N])
        for i in range(high):
            G_high[i,i] = 1
    return G_high

def _freq_selector(Xhat, eigval, eigvec):
    Xhat = np.array(Xhat)
    PSD = Xhat**2 
    PSD = np.mean(PSD, axis=1)
    AUCTOT = np.trapz(PSD, eigval)
    i=1
    AUC=0
    while AUC < AUCTOT/2:
        AUC = np.trapz(PSD[0:i], eigval[0:i])
        i += 1 
    return i-1 

def _Group_avg_SC(subjects):
    TSC=[]
    for subj in subjects: 
        sc = read_SC(subj)
        TSC.append(sc)
    TSC = np.array(TSC)
    return np.mean(TSC, axis=0)
    

def _Group_freq_selector(subjects):
    TPSD = np.zeros((360, len(subjects)))
    sc = _Group_avg_SC(subjects)
    eigval, eigvec = Graph_creation(sc, kind='normalized')
    U = eigvec
    c = 0
    for subj in subjects:
        rest = read_rs(subj)
        Xhat = np.matmul(np.transpose(U), rest)
        Xhat = np.array(Xhat)
        PSD = Xhat**2 
        PSD = np.mean(PSD, axis=1)
        TPSD[:,c] = PSD
        c=c+1
    MPSD = np.mean(TPSD, axis=1)
    AUCTOT = np.trapz(MPSD, eigval)
    i=1
    AUC=0
    while AUC < AUCTOT/2:
        AUC = np.trapz(MPSD[0:i], eigval[0:i])
        i += 1 
    return i-1

def _PSD_spliter(X, subjects, N):
    res=np.zeros(N+1)
    res[0]=0
    res[-1]=359
    sc = _Group_avg_SC(subjects)
    eigval, eigvec = Graph_creation(sc, kind='normalized')
    U = eigvec
    c = 0
    TPSD = np.zeros((360, len(subjects)))
    for i in range(len(subjects)): 
        tmp = X[:,i*1190:(i+1)*1190]
        Xhat = np.matmul(np.transpose(U), tmp)
        Xhat = np.array(Xhat)
        PSD = Xhat**2 
        PSD = np.mean(PSD, axis=1)
        TPSD[:,i] = PSD
    MPSD = np.mean(TPSD, axis=1)
    SPSD = np.std(TPSD, axis=1)
    AUCTOT = np.trapz(MPSD, eigval)
    cum_power = []
    for i in range(1,len(MPSD)+1):
        cum_power.append(np.trapz(MPSD[0:i], eigval[0:i]))
    cum_power = np.array(cum_power)
    for i in range(1,N):
        res[i] = np.max(np.where((cum_power>(i-1)*AUCTOT/N) & (cum_power<i*AUCTOT/N)))
    
    return res, MPSD, SPSD

import pygsp 
from pygsp import graphs
import numpy as np
import matplotlib.pylab as plt 
import seaborn as sns 
from scipy.io import loadmat 
#import alphashape
#from descartes import PolygonPatch
import sys
from numpy import linalg as LA
from io import *

################################################################################

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

def SDI_calculation(X, eigval, eigvec, kind):
    U = eigvec
    X_hat = np.matmul(np.transpose(U), X)
    cut_lambda = _freq_selector(X_hat, eigval, eigvec)
    Xc = LPF(X, eigval, eigvec, cut_lambda, kind)
    Xd = HPF(X, eigval, eigvec, len(eigval)-cut_lambda, kind)
    NXc = LA.norm(Xc, axis=1)
    NXd = LA.norm(Xd, axis=1)
    return NXd/NXc, cut_lambda

def Group_SDI_calculation(subjects, data_address, kind):
    cut_lambda = _Group_freq_selector(subjects, data_address)
    sc = Group_avg_SC(subjects, data_address)
    eigval, eigvec = Graph_creation(sc, kind='normalized')
    U = eigvec
    #TSDI = []
    NXC = []
    NXD = []
    for subj in subjects:
        rest = read_rs(subj, data_address)
        Xc = LPF(rest, eigval, eigvec, cut_lambda, kind)
        Xd = HPF(rest, eigval, eigvec, len(eigval)-cut_lambda,kind)
        NXc = np.linalg.norm(Xc, axis=1)
        NXd = np.linalg.norm(Xd, axis=1)
        NXC.append(NXc)
        NXD.append(NXd)
        #SDI = NXd/NXc
        #TSDI.append(SDI)
    return np.mean(NXD, 0)/np.mean(NXC, 0)

def dSDI_calculation(X, eigval, eigvec, kind): 
    U = eigvec
    X_hat = np.matmul(np.transpose(U), X)
    cut_lambda = _freq_selector(X_hat, eigval, eigvec)
    Xc = LPF(X, eigval, eigvec, cut_lambda, kind)
    Xd = HPF(X, eigval, eigvec, len(eigval)-cut_lambda, kind)
    return Xd/Xc, cut_lambda

def Dir_Energy(sig, sc):
    R = sc.shape[0]
    E = np.zeros(R)
    for r in range(R):
        tmp = 0
        for s in range(R):
            tmp += sc[r,s]*((sig[r] - sig[s])**2)
        E[r] = tmp
    return E.reshape(R,1)

def MDE(sig, sc, network, nec_files_addr):
    R = sc.shape[0]
    ind_file = loadmat(nec_files_addr+'match_Yeo.mat')
    M = ind_file['ind_'+network]
    M = M-1
    mde = 0
    for r in range(len(M)):
        for s in range(R):
            mde += sc[M[r],s]*((sig[M[r]] - sig[s])**2)
    return mde

def BMDE(sig, sc, network1, network2, nec_files_addr):
    R = sc.shape[0]
    ind_file = loadmat(nec_files_addr+'match_Yeo.mat')
    M1 = ind_file['ind_'+network1]
    M1 = M1 -1
    M2 = ind_file['ind_'+network2]
    M2 = M2 - 1 
    bmde = 0 
    for r in range(len(M1)):
        for s in range(len(M2)):
            bmde += sc[M1[r],M2[s]]*((sig[M1[r]] - sig[M2[s]])**2)
    return bmde 
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
    

def _Group_freq_selector(subjects, data_address):
    TPSD = np.zeros((360, len(subjects)))
    sc = Group_avg_SC(subjects, data_address)
    eigval, eigvec = Graph_creation(sc, kind='normalized')
    U = eigvec
    c = 0
    for subj in subjects:
        rest = read_rs(subj,data_address)
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

def _PSD_spliter(X, subjects, N, data_address):
    res=np.zeros(N+1)
    res[0]=0
    res[-1]=359
    sc = Group_avg_SC(subjects, data_address)
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

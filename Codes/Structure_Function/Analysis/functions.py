import os
import pickle
import os.path as op 
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns 
from scipy import signal
from scipy.io import savemat
from scipy.stats import zscore
from sklearn.cluster import KMeans
from os.path import join
from itertools import groupby
    
def _bp_filter(sig, fmin, fmax, TR, order): 
    ## bandpass filtering of BOLD signal 
    fnyq = 1/(2*TR)  
    wn = np.array([fmin/fnyq, fmax/fnyq]) 
    [b, a] = signal.butter(order, wn, btype='band')
    return signal.lfilter(b, a, sig, axis=1)

def _adif(phi1, phi2):
    ## Performs phase wrapping 
    if abs(phi2-phi1)>np.pi:
        return 2*np.pi-abs(phi2-phi1)
    else:
        return abs(phi2-phi1)

def pattern_generation(subj_list, flights, tps, data_dir, out_dir, atlas):
    if atlas == 'AAL3':
        R = 170
    elif atlas == 'HO':
        R = 63
    elif atlas == 'SCH_100':
        R = 119
    idx = np.where(np.tril(np.ones((R, R)), k=-1)) # Indices of the lower triangle of connectivity matrices 
    Patterns = np.zeros((int(R*(R-1)/2), 1)) # Variable to store connectivity patterns of all time points
    c=0 # Subject Counter 
    subj_pointer = []
    group_pointer = []
    for subj in subj_list:
        for f in flights: 
            for tp in tps:
                if op.isdir(op.join(data_dir, subj, f, tp)):
                    if f'{subj}-{tp}' != 'cosm-03-post':
                        print(f'{subj} - {f} - {tp}')
                        sig_dir = op.join(data_dir, subj, f, tp, f'Sig_{atlas}.txt')
                        sig = np.loadtxt(sig_dir)
                        sig = zscore(sig, axis=1)
                        sig = np.float16(sig)
                        # Bandpass Filtering 
                        sigc = np.concatenate((np.fliplr(sig), sig, np.fliplr(sig)), axis=1)
                        sigf = _bp_filter(sigc, fmin=0.01, fmax=0.04, TR=2, order=4)    
                        sig = sigf[:,sig.shape[1]:sig.shape[1]*2]
                        # Hilbert transform of whole signal and calculating instantaneous phase 
                        phi = np.angle(signal.hilbert(sig, axis=1))
                        T = sig.shape[1]
                        # Calculating connectivity patterns at each time point 
                        for t in range(T): 
                            patt = np.zeros((R, R))
                            for i in range(R): 
                                for j in range(i): 
                                    patt[i, j] = np.cos(_adif(phi[i,t], phi[j,t]))
                            Patterns = np.concatenate((Patterns, patt[idx].reshape((len(patt[idx]),1))), axis=1)
                            subj_pointer.append(subj)
                            gp = subj.split('-')[0]
                            group_pointer.append(f'{gp}-{tp}')
    Patterns = Patterns[:,1:]
    Patterns = np.float16(Patterns)
    my_dict = {
        'Patterns':Patterns,
        'Subject_Pointer':subj_pointer, 
        'Group_Pointer': group_pointer
    } 
    res_dir = op.join(out_dir, f'conn_matrices')
    if op.isdir(res_dir)==0: 
        os.makedirs(res_dir)
    savemat(os.path.join(res_dir, f'Patterns_{atlas}_all.mat'), my_dict)
    return Patterns

def cluster_generation(patterns, out_dir, cluster_num, atlas):
    if atlas == 'AAL3':
        R = 170
    elif atlas == 'HO':
        R = 63
    elif atlas == 'SCH_100':
        R = 119
    kmeans = KMeans(n_clusters=cluster_num, n_init=500, max_iter=200).fit(patterns)
    res_dir = op.join(out_dir, f'clustering_results')
    if op.isdir(res_dir)==0: 
        os.makedirs(res_dir)
    pickle.dump(kmeans, open(os.path.join(res_dir, f'kmeans_model_{atlas}.sav'), 'wb'))
    return kmeans 

def label_corrector(main_dir, order, atlas):
    if atlas == 'AAL3':
        R = 170
    elif atlas == 'HO':
        R = 63
    elif atlas == 'SCH_100':
        R = 119
    kmeans_dir = join(main_dir, f'clustering_results')
    model = pickle.load(open(join(kmeans_dir, f'kmeans_model_{atlas}.sav'), 'rb'))     
    corr_main_patterns = model.cluster_centers_.T[:, order-1]
    labels = model.labels_
    for i in range(len(labels)):
        labels[i] = np.where(order==labels[i]+1)[0][0]           
    corr_labels = labels
    return corr_main_patterns, corr_labels


def schaefer_conn_plotting(mat):
    f, a = plt.subplots(1,1, figsize=(7,5))
    sns.heatmap(mat, cmap='RdBu_r', vmin=-1, vmax=1, ax=a)
    a.set_xticks([])
    a.set_yticks([])
    a.hlines(0, 0, 24, color='k', linewidth=2)
    a.hlines(24, 0, 24, color='k', linewidth=2)
    a.vlines(0, 0, 24, color='k', linewidth=2)
    a.vlines(24, 0, 24, color='k', linewidth=2)
    a.text(6,-2,'DMN', size=14)
    a.text(6,-2,'DMN', size=14)
    a.hlines(37, 24, 37, color='k', linewidth=2)
    a.hlines(24, 24, 37, color='k', linewidth=2)
    a.vlines(37, 24, 37, color='k', linewidth=2)
    a.vlines(24, 24, 37, color='k', linewidth=2)
    a.text(24,-2,'Cont', size=14)
    a.text(24,-2,'Cont', size=14)
    a.hlines(51, 37, 51, color='k', linewidth=2)
    a.hlines(37, 37, 51, color='k', linewidth=2)
    a.vlines(51, 37, 51, color='k', linewidth=2)
    a.vlines(37, 37, 51, color='k', linewidth=2)
    a.text(40,-2,'SM', size=14)
    a.text(40,-2,'SM', size=14)
    a.hlines(57, 51, 57, color='k', linewidth=2)
    a.hlines(51, 51, 57, color='k', linewidth=2)
    a.vlines(57, 51, 57, color='k', linewidth=2)    
    a.vlines(51, 51, 57, color='k', linewidth=2)
    a.text(51,-2,'Lm', size=14)
    a.text(51,-2,'Lm', size=14)
    a.hlines(68, 57, 68, color='k', linewidth=2)
    a.hlines(57, 57, 68, color='k', linewidth=2)
    a.vlines(68, 57, 68, color='k', linewidth=2)
    a.vlines(57, 57, 68, color='k', linewidth=2)
    a.text(61,-2,'VA', size=14)
    a.text(61,-2,'VA', size=14) 
    a.hlines(84, 68, 84, color='k', linewidth=2)
    a.hlines(68, 68, 84, color='k', linewidth=2)
    a.vlines(84, 68, 84, color='k', linewidth=2)
    a.vlines(68, 68, 84, color='k', linewidth=2)
    a.text(72,-2,'DA', size=14)
    a.text(72,-2,'DA', size=14)
    a.hlines(100, 84, 100, color='k', linewidth=2)
    a.hlines(84, 84, 100, color='k', linewidth=2)
    a.vlines(100, 84, 100, color='k', linewidth=2)
    a.vlines(84, 84, 100, color='k', linewidth=2)
    a.text(86,-2,'Vis', size=14)
    a.text(86,-2,'Vis', size=14)

    a.text(-6,18,'DMN', size=14, rotation=90)
    a.text(-6,35,'Cont', size=14, rotation=90)
    a.text(-6,45,'SM', size=14, rotation=90)
    a.text(-6,55,'Lm', size=14, rotation=90)
    a.text(-6,65,'VA', size=14, rotation=90)    
    a.text(-6,77,'DA', size=14, rotation=90)
    a.text(-6,94,'Vis', size=14, rotation=90)
    
    return f 

def mdt_calc(labs, element):
    labs2 = [labs[0]]
    for i in range(1,len(labs)):
        if labs[i] != labs2[-1]:
            labs2 = labs2 + [labs[i]]

    count_dups = [sum(1 for _ in group) for _, group in groupby(labs)]
    idx = np.where(np.array(labs2)==element)
    dt = np.array(count_dups)[idx]
    mdt = np.mean(dt)

    return mdt


def sc_average(group, flight, time, atlas, data_dir):
    if atlas == 'AAL3':
        R = 170
    elif atlas == 'HO':
        R = 63
    elif atlas == 'SCH_100':
        R = 119

    all_subj = os.listdir(data_dir)
    subjects = [all_subj[i] for i in range(len(all_subj)) if all_subj[i].startswith(group)]
    
    SCT = np.zeros((R,R,1))
    for sub in subjects: 
        path = op.join(data_dir, sub, flight, time)
        if op.isdir(path):
            sc = np.loadtxt(op.join(data_dir, sub, flight, time, f'SC_Norm_{atlas}.csv')
                       , dtype=str, delimiter=',')
            tmpsc = np.zeros(sc.shape)
            for i in range(sc.shape[0]): 
                for j in range(sc.shape[1]):
                    tmpsc[i,j] = float(sc[i,j])
            SCT = np.concatenate((SCT,tmpsc.reshape(R,R,1)), axis=2)
    SCT = SCT[:,:,1:]
    avg_sc = np.mean(SCT, axis=2)
    
    return avg_sc

def cutoff_estimation(group, flight, time, atlas, data_dir, e, U):
    if atlas == 'AAL3':
        R = 170
    elif atlas == 'HO':
        R = 63
    elif atlas == 'SCH_100':
        R = 119
        
    all_subj = os.listdir(data_dir)
    subjects = [all_subj[i] for i in range(len(all_subj)) if all_subj[i].startswith(group)]
    TPSD = np.zeros((R,1))
    for sub in subjects:
        path = op.join(data_dir, sub, flight, time)
        if op.isdir(path):
            rs_file = op.join(data_dir, sub, flight, time, f'Sig_{atlas}.txt')
            rs = np.loadtxt(rs_file)
            nrs = zscore(rs, axis=1)
            nrs = np.nan_to_num(nrs)
            Xhat = np.matmul(np.transpose(U), nrs)
            Xhat = np.array(Xhat)
            PSD = Xhat**2 
            PSD = np.mean(PSD, axis=1)
            TPSD = np.concatenate((TPSD, PSD.reshape((R,1))), axis=1)
    TPSD = TPSD[:, 1:]
    MPSD = np.mean(TPSD, axis=1) 
    SPSD = np.std(TPSD, axis=1)
    
    auctot = np.trapz(MPSD)
    auc = 0
    c = 0
    while auc<auctot/2: 
        auc = np.trapz(MPSD[0:c])
        c = c+1
    freq_cut = e[c-1]
    
    return c-1, freq_cut, MPSD, SPSD
import numpy as np
from scipy.io import loadmat  
from sklearn.cluster import KMeans
#import alphashape
import numpy.linalg as LA
#import entropy
#import nolds

def boundary(data, alpha):
    points =[]
    for i in range(len(data)):
        x = (data[i,0], data[i,1])
        points.append(x)
    alpha_shape = alphashape.alphashape(points, alpha)
    return alpha_shape, alpha_shape.area

def transition_probability(labels):
    n_states = len(np.unique(labels))
    P = np.zeros((n_states, n_states))
    for i in range(n_states): 
        for j in range(n_states):
            denum = np.sum(labels[0:len(labels)-1]==i)
            num = 0
            for k in range(len(labels)-1):
                if labels[k] == i and labels[k+1] ==j:
                    num += 1 
            P[i,j] = num/denum
    return P 

def transition_probability_sec(labels): 
    n_states = len(np.unique(labels))
    P = np.zeros((n_states,n_states,n_states))
    tmp = transition_probability(labels)
    for i in range(n_states):
        for j in range(n_states):
            for k in range(n_states):
                #denum = 0
                #for c in range(len(labels)-2):
                #    if labels[c] == i and labels[c+1] == j:
                #        denum += 1 
                #num = 0 
                #for c in range(len(labels)-2):
                #    if labels[c] == i and labels[c+1] == j and labels[c+2]==k:
                #        num  += 1
                #if denum ==0:
                #    P[i,j,k] = 0
                #else:
                P[i,j,k] = tmp[i,j]*tmp[j,k]
    return P 

def path_length(data):
    L = np.shape(data)[0]
    D = 0
    for i in range(L-1):
        x1 = data[i,:]
        x2 = data[i+1,:]
        dist = np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2)
        D = D + dist
    return D 

def mean_path_length(data):
    L = np.shape(data)[0]
    D = []
    for i in range(L-1):
        x1 = data[i,:]
        x2 = data[i+1,:]
        dist = np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2)
        D.append(dist)
    return np.mean(np.array(D))

def conv_hull_area(data):
    alpha = 0 
    sub_boundary, area = boundary(data, alpha)
    return area 

def state_prob(labels): 
    n_states = len(np.unique(labels))
    prob = np.zeros(n_states)
    for i in range(n_states):
        prob[i] = np.sum(labels == i)/len(labels)
    return prob 

# Entropy Measures 

def shannon_entropy(labels):
    prob = state_prob(labels)
    n_states = len(prob)
    rep = np.zeros(n_states)
    for i in range(n_states):
        rep[i] = np.sum(labels == i)
    H = 0 
    for i in range(n_states): 
        H = H - rep[i]*prob[i]*np.log2(prob[i])
    return H

def permutation_entropy(labels, order):
    return entropy.perm_entropy(labels, order)

def spectral_entropy(labels, order):
    return entropy.spectral_entropy(labels, order)

def svd_entropy(labels, order, delay):
    return entropy.svd_entropy(labels, order, delay)

def approximate_entropy(labels, order):
    return entropy.app_entropy(labels, order)

def sample_entropy(labels, order): 
    return entropy.sample_entropy(labels, order)

def lempelzif_entropy(labels):
    seq = ''
    for i in range(len(labels)):
        if labels[i] == 0:
            seq = seq + '00'
        elif labels[i] == 1:
            seq = seq + '01'
        elif labels[i] == 2: 
            seq = seq + '10'
        else: 
            seq = seq + '11'
    return entropy.lziv_complexity(seq)

# Fractal Measures 

def fd_petrosian(labels):
    return entropy.petrosian_fd(labels)

def fd_katz(labels):
    return entropy.katz_fd(labels)

def fd_higuchi(labels, kmax):
    return entropy.higuchi_fd(labels, kmax)

def fd_det_fluc(labels):
    return entropy.detrended_fluctuation(labels)

def hurst_exponent(labels):
    return nolds.hurst_rs(labels)

def corr_dimension(labels, embdim):
    return nolds.corr_dim(labels, embdim)


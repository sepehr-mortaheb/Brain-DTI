
import numpy as np
import matplotlib.pylab as plt 

#############################################################################

def plot_conn_patterns(conn_patt, length, ax):
    cumlength = np.cumsum(length)
    cumlength = cumlength[0:-1]
    ax.imshow(conn_patt, cmap='jet', vmin=-1, vmax=1)
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
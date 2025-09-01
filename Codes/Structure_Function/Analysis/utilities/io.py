from scipy.io import loadmat 
from os import listdir 
import numpy as np

###########################################################################

def get_subjects_id(data_address):
    drct = data_address + 'Structural'
    allfiles = listdir(drct)
    files = [f for f in allfiles if not f.startswith('.')]
    subj_id = [f.split('_')[-1].split('.')[0] for f in files]
    return subj_id

def read_SC(subj_id, data_address):
    sc_file = 'Connectome_sb_'+subj_id+'.mat'
    sc = loadmat(data_address+'Structural/'+sc_file)['Cnorm']
    return sc 

def read_rs(subj_id, data_address):
    rs_file = 'Rest1LR_Sub'+subj_id+'_Glasser.mat'
    rs = loadmat(data_address+'Functional/Rest/Rest1LR/'+rs_file)['TCS']
    return rs 

def Group_avg_SC(subjects, data_address):
    TSC=[]
    for subj in subjects: 
        sc = read_SC(subj, data_address)
        TSC.append(sc)
    TSC = np.array(TSC)
    return np.mean(TSC, axis=0)

def read_task(subj_id, task, data_addr):
    par_file_addr = data_addr+'Functional/Task/'+task+'/Paradigm_'+task+'.mat'
    task_addr = data_addr+'Functional/Task/'+task+'/Task'+task[0:2]+'_Sub'+\
        subj_id+'_Glasser.mat'
    tasksig = loadmat(task_addr)
    tasksig = tasksig['TCS']
    raw_par = loadmat(par_file_addr)
    par = []
    for i in range(tasksig.shape[1]):
        par.append(raw_par[task][0][0][1][0][i][0])
    return tasksig, par 

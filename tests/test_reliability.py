# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
import numpy as np
import Functional_Fusion.reliability as rel
import matplotlib.pyplot as plt



def sim_data(n_subjects=10,
             n_part=4,
             n_conditions=5,
             n_vox=200,
             var_comp = [1,1,1]):
    """Simulate data with group, subject and error variance.
    Args:
        n_subjects (int): Number of subjects.
        n_conditions (int): Number of conditions.
        n_part (int): Number of partitions.
        n_vox (int): Number of voxels.
        var_comp (list): List of variance components.
    Returns:
        (nsubj,npart,ncond,nvox) array of simulated data.
    """
    # generate a 4-d tensor of data
    group = np.random.normal(0, np.sqrt(var_comp[0]), size=(1,1,n_conditions,n_vox))
    subject = np.random.normal(0, np.sqrt(var_comp[1]), size=(n_subjects,1,n_conditions,n_vox))
    error = np.random.normal(0, np.sqrt(var_comp[2]), size=(n_subjects,n_part,n_conditions,n_vox))
    data = group + subject + error
    # Flatten the tensor to 3-d
    data = data.reshape(n_subjects,n_part*n_conditions,n_vox)
    cond_vec = np.tile(np.arange(n_conditions),n_part)
    part_vec = np.repeat(np.arange(n_part),n_conditions)
    return data,cond_vec,part_vec

def compare_within():
    data,cond_vec,part_vec = sim_data(var_comp=[1,0.2,1],n_subjects=10,n_part=3,n_conditions=5,n_vox=1000)
    # data=data+1
    var_est1 = rel.decompose_subj_group(data,cond_vec,part_vec,separate='subject_wise',subtract_mean=False)
    print((var_est1[:,0]+var_est1[:,1])/var_est1.sum(axis=1))
    var_est2 = rel.within_subj(data,cond_vec,part_vec,separate='none',subtract_mean=False)
    print(var_est2)
    pass

def compare_between():
    data,cond_vec,part_vec = sim_data(var_comp=[1,0,1],n_subjects=10,n_part=1,n_conditions=5,n_vox=1000)
    # data=data+1
    var_est1 = rel.decompose_subj_group(data,cond_vec,part_vec,separate='subject_wise',subtract_mean=True)
    print(var_est1[:,0]/(var_est1[:,1]+var_est1[:,2]))
    var_est2 = rel.between_subj(data,cond_vec,separate='condition_wise',subtract_mean=True)
    print(var_est2)


def test_size():
    data,cond_vec,part_vec = sim_data(var_comp=[1,0,1],n_subjects=10,n_part=8,n_conditions=5,n_vox=200)
    var_est1 = rel.decompose_subj_group(data,cond_vec,part_vec,separate='subject_wise',subtract_mean=True)
    print(var_est1.shape)


    var_est2 = rel.within_subj(data,cond_vec,part_vec,separate='condition_wise',subtract_mean=True)
    print(var_est2.shape)

    var_est3 = rel.within_subj(data,cond_vec,part_vec,separate='voxel_wise',subtract_mean=True)
    print(var_est3.shape)

    var_est4 = rel.between_subj(data,cond_vec,separate='condition_wise',subtract_mean=True)
    print(var_est4.shape)

    var_est5 = rel.between_subj(data,cond_vec,separate='voxel_wise',subtract_mean=True)
    print(var_est5.shape)

    var_est6 = rel.within_subj_loo(data,cond_vec,part_vec,separate='voxel_wise',subtract_mean=True)
    print(var_est6.shape)

    var_est7 = rel.between_subj_loo(data,cond_vec,separate='none',subtract_mean=True)
    print(var_est7.shape)
    pass

def test_between_subj_loo(var=0.5, N=10):
    data,cond_vec,part_vec = sim_data(var_comp=[var,1-var,0],n_subjects=N,n_part=1,n_conditions=1000,n_vox=1)
    rbs = rel.between_subj(data,cond_vec,separate='none',subtract_mean=True)
    low = rel.between_subj_loo(data,cond_vec,separate='none',subtract_mean=True)
    high = rel.between_subj_avrg(data,cond_vec,separate='none',subtract_mean=True)
    nc = np.sqrt(rbs.mean())
    ncl = low.mean()
    ncu =  high.mean()
    return nc, ncl, ncu

def plot_noiseceil_relationship():
    N = np.array([2,6,10,30], dtype=int)
    nc = np.zeros(len(N))
    ncl = np.zeros(len(N))
    ncu = np.zeros(len(N))
    for i,n in enumerate(N):
        print(f"Testing N={n}")
        nc[i], ncl[i], ncu[i] = test_between_subj_loo(var=0.5,N=n)
    plt.plot(N, nc, '-', color='black')
    plt.plot(N, ncu, ':', color='red')
    plt.plot(N, ncl, ':', color='blue')
    pass

if __name__ == "__main__":
    plot_noiseceil_relationship()

    # test_between_subj_loo()


    pass
# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
import numpy as np
import Functional_Fusion.reliability as rel
import matplotlib.pyplot as plt



def sim_data(n_subjects=10,
             n_part=4,
             n_conditions=5,
             n_vox=20,
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


if __name__ == "__main__":
    data,cond_vec,part_vec = sim_data(var_comp=[1,0.2,1],n_subjects=10,n_part=8,n_conditions=5,n_vox=20)
    data=data+1
    var_est1 = rel.decompose_subj_group(data,cond_vec,part_vec,separate='none',subtract_mean=False)
    var_est2 = rel.decompose_subj_group(data,cond_vec,part_vec,separate='none',subtract_mean=True)
    var_est3 = rel.decompose_subj_group(data,cond_vec,part_vec,separate='voxel_wise',subtract_mean=True)



    pass
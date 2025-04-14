# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
from pathlib import Path
import numpy as np
import Functional_Fusion.dataset as ds
import paths as paths
import Functional_Fusion.connectivity as conn
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
    group = np.random.normal(0, np.sqrt(var_comp[0]), size=(1,1,n_conditions,n_vox))
    subject = np.random.normal(0, np.sqrt(var_comp[1]), size=(n_subjects,1,n_conditions,n_vox))
    error = np.random.normal(0, np.sqrt(var_comp[2]), size=(n_subjects,n_part,n_conditions,n_vox))
    data = group + subject + error
    return data






if __name__ == "__main__":
    data = sim_data(var_comp=[1,1,1])
    var_est1 = ds.decompose_pattern_into_group_indiv_noise(data,criterion='global')
    var_subj = ds.decompose_pattern_into_group_indiv_noise(data,criterion='subject_wise')
    var_est2 = ds.decompose_pattern_into_group_indiv_noise(data,criterion='voxel_wise')
    var_est3 = ds.decompose_pattern_into_group_indiv_noise(data,criterion='condition_wise')
    pass
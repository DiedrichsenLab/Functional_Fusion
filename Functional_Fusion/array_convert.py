import numpy as np
import pandas as pd
import nibabel as nb 

def flat2ndarray(flat_data, cond_vec, part_vec):
    """
    convert flat data (n_subjects x n_trials x n_voxels) into a 4d ndarray (n_subjects x n_partitions x n_conditions x n_voxels)

    Args:
        flat_data:
        part_vec:
        cond_vec:

    Returns:
        data

    """

    [n_subjects, n_trials, n_voxels] = flat_data.shape

    unique_conditions = np.unique(cond_vec)
    unique_partitions = np.unique(part_vec)

    n_conditions = unique_conditions.size
    n_partitions = unique_partitions.size
   
    data = np.zeros((n_subjects, n_partitions, n_conditions, n_voxels))

    for p,partI in enumerate(unique_partitions):
        for c,condI in enumerate(unique_conditions):
            trial_inds = np.where(np.logical_and(cond_vec == condI, part_vec == partI))
            data[:, p, c, :] = flat_data[:, trial_inds, :].squeeze()
            
    return data

def get_structure_data(structure='pontine', data_dir='/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/RegionOfInterest_BOLDMNI/data/group'):
    T = pd.read_csv('/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/participants.tsv', sep='\t')
    A = []
    for i, good_value in zip(T.participant_id, T.good):
        if good_value==1:
            file_path = f'{data_dir}/beta_glm2_{structure}_{i}.dscalar.nii'
            cifti = nb.load(file_path)
            A.append(cifti.get_fdata())
        to_tensor = np.array(A)
    return to_tensor
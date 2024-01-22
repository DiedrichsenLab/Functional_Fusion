# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetMDTB
import Functional_Fusion.dataset as ds
import nibabel as nb
import subprocess
import paths as paths
import Functional_Fusion.connectivity as conn
import matplotlib.pyplot as plt


base_dir = paths.set_base_dir()
atlas_dir = paths.set_atlas_dir(base_dir)
dname = 'MDTB'
data_dir = base_dir + '/' + dname
atlas_dir = base_dir + '/Atlases'

def flat2ndarray(flat_data, part_vec, cond_vec):
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

    unique_partitions = np.unique(part_vec)
    n_partitions = unique_partitions.size

    unique_conditions = np.unique(cond_vec)
    n_conditions = unique_conditions.size

    data = np.zeros((n_subjects, n_partitions, n_conditions, n_voxels))

    for partI in np.arange(n_partitions):
        for condI in np.arange(n_conditions):
            trial_inds = np.where(np.logical_and(cond_vec == unique_conditions[condI], part_vec == unique_partitions[partI]))
            data[:, partI, condI, :] = np.nanmean(flat_data[:, trial_inds, :], axis=1).squeeze()

    return data

def reshape_data(data, info, cond_column='cond_num_uni', part_column='run', mean_centering=False):
    """Reshape data from (n_subjects, n_trials, n_voxels) to (n_subjects, n_runs, n_conditions, n_voxels) to comply with the decompose_pattern_into_group_indiv_noise function."""
    # Extract each run and concatenate in third dimension
    data_reshaped = np.zeros((data.shape[0], info[part_column].max(), info[cond_column].max(), data.shape[-1]))
    for i in range(1, info.run.max()+1):
        run_data = data[:, info.run == i, :]
        data_reshaped[:, i-1, :, :] = run_data

    # Set nans to number and print number of nans
    print(f'Setting {np.sum(np.isnan(data_reshaped))} nan values to zero.')
    data_reshaped = np.nan_to_num(data_reshaped)
    
    # Mean centering
    if mean_centering:
        # Subtract the mean across subjects and runs from each voxel
        mean_across_conditions = np.mean(data_reshaped, axis=(2))
        # Repeat the mean across conditions for each condition
        mean_across_conditions = np.repeat(mean_across_conditions[:, :, np.newaxis, :], data_reshaped.shape[2], axis=2)
        data_reshaped = data_reshaped - mean_across_conditions
    
    return data_reshaped

if __name__ == "__main__":
    mdtb_dataset = DataSetMDTB(data_dir)

    # -- Get resting-state subjects --
    T = pd.read_csv(
        data_dir + '/participants.tsv', delimiter='\t')
    subject_subset = T.participant_id[T['ses-rest'] == 1].tolist()

    # -- Variance decomposition --
    mean_centering = True
    # --- Resting-state ---
    data_rest, info_rest = mdtb_dataset.get_data(ses_id='ses-rest', type='Net69Run', space='MNISymC2', subj=subject_subset)
    data_reshaped_rest = reshape_data(data_rest, info_rest, cond_column='net_id', mean_centering=mean_centering)
    vars_rest = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_rest, criterion='global')
    print(f'Rest variance\nGroup: {vars_rest[0][0]:.2f}\nSubject: {vars_rest[0][1]:.2f}\nError: {vars_rest[0][2]:.2f}')

    # --- Task ---
    data_task, info_task = mdtb_dataset.get_data(ses_id='ses-s1', type='CondRun', space='MNISymC2', subj=subject_subset)
    data_reshaped_task = reshape_data(data_task, info_task, cond_column=mdtb_dataset.cond_ind, mean_centering=mean_centering)
    vars_task = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_task, criterion='global')
    print(f'Task variance\nGroup: {vars_task[0][0]:.2f}\nSubject: {vars_task[0][1]:.2f}\nError: {vars_task[0][2]:.2f}')
    
    # --- Task Neocortex ---
    data_task_neocortex, info_task_neocortex = mdtb_dataset.get_data(ses_id='ses-s1', type='CondRun', space='fs32k')
    data_reshaped_task_neocortex = reshape_data(data_task_neocortex, info_task, cond_column=mdtb_dataset.cond_ind, mean_centering=mean_centering)
    data_reshaped_task_neocortex = reshape_data(data_task_neocortex, info_task, cond_column=mdtb_dataset.cond_ind, mean_centering=False)
    
    # Create dataset all runs are data from the first run for each subject (to maximize subject variance)
    data_max_subject_variance = np.zeros_like(data_reshaped_task)
    for i in range(data_reshaped_task.shape[1]):
        data_max_subject_variance[:, i, :, :] = data_reshaped_task[:, 0, :, :]
    vars_subject = ds.decompose_pattern_into_group_indiv_noise(data_max_subject_variance, criterion='global')
    print(f'Same runs variance\nGroup: {vars_subject[0][0]:.2f}\nSubject: {vars_subject[0][1]:.2f}\nError: {vars_subject[0][2]:.2f}')

    # Create dataset where all 24 subjects are data from the first subject (to maximize group variance)
    data_max_group_variance = np.zeros_like(data_reshaped_task)
    for i in range(data_reshaped_task.shape[0]):
        data_max_group_variance[i, :, :, :] = data_reshaped_task[0, :, :, :]
    vars_group = ds.decompose_pattern_into_group_indiv_noise(data_max_group_variance, criterion='global')
    print(f'Same subjects variance\nGroup: {vars_group[0][0]:.2f}\nSubject: {vars_group[0][1]:.2f}\nError: {vars_group[0][2]:.2f}')

    # Create dataset where error variance is maximized
    data_max_error_variance = np.random.rand(*data_reshaped_task.shape)
    vars_error = ds.decompose_pattern_into_group_indiv_noise(data_max_error_variance, criterion='global')
    print(f'Random numbers variance\nGroup: {vars_error[0][0]:.2f}\nSubject: {vars_error[0][1]:.2f}\nError: {vars_error[0][2]:.2f}')


    vars_task_neo = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_task_neocortex, criterion='global')
    print(f'Neocortical task variance\nGroup: {vars_task_neo[0][0]:.2f}\nSubject: {vars_task_neo[0][1]:.2f}\nError: {vars_task_neo[0][2]:.2f}')

    pass
    
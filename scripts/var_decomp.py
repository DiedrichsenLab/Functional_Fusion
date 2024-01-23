# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetMDTB, DataSetPontine
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
    # Extract each partition and concatenate in third dimension
    data_reshaped = np.zeros((data.shape[0], info[part_column].max(), info[cond_column].max(), data.shape[-1]))
    for i in range(1, info[part_column].max()+1):
        print(i)
        part_data = data[:, info[part_column] == i, :]
        data_reshaped[:, i-1, :, :] = part_data

    # Set nans to number and print number of nans
    print(f'Setting {np.sum(np.isnan(data_reshaped))} nan values to zero.')
    data_reshaped = np.nan_to_num(data_reshaped)
    
    # Mean centering
    if mean_centering:
        # Subtract the mean across subjects and partition from each voxel
        mean_across_conditions = np.mean(data_reshaped, axis=(2))
        # Repeat the mean across conditions for each condition
        mean_across_conditions = np.repeat(mean_across_conditions[:, :, np.newaxis, :], data_reshaped.shape[2], axis=2)
        data_reshaped = data_reshaped - mean_across_conditions
    
    return data_reshaped



if __name__ == "__main__":
    

    # -- Get resting-state subjects --
    T = pd.read_csv(
        data_dir + '/participants.tsv', delimiter='\t')
    subject_subset = T.participant_id[T['ses-rest'] == 1].tolist()

    # Reproduce derricks results
    X_individuals_task, info_individuals_task, dataset_obj_individuals = ds.get_dataset(base_dir,
                                                                              dataset='MDTB',
                                                                              atlas='MNISymC2',
                                                                              subj=None,
                                                                              sess='all',
                                                                              type='CondHalf')
    X_individuals_rest, info_individuals_rest, dataset_obj_individuals = ds.get_dataset(base_dir,
                                                                              dataset='MDTB',
                                                                              atlas='MNISymC2',
                                                                              subj=subject_subset,
                                                                              sess='ses-rest',
                                                                              type='Net69Run')

    #  -- Task s1 --
    # X_individuals_task_ses1 = X_individuals_task[:, info_individuals_task['sess'] == 'ses-s1', :]
    # data = X_individuals_task_ses1
    # task_conds = list(info_individuals_task[info_individuals_task['sess'] == 'ses-s1'].names)

    # half1_inds = np.array([x for x in np.arange(len(task_conds)) if task_conds[x].__contains__('half1')])
    # half2_inds = np.setdiff1d(np.arange(len(task_conds)), half1_inds)

    #  -- Task all --
    data = X_individuals_task
    task_conds = list(info_individuals_task.names)

    half1_inds = np.array([x for x in np.arange(len(task_conds)) if task_conds[x].__contains__('half1')])
    half2_inds = np.setdiff1d(np.arange(len(task_conds)), half1_inds)
    # -- Rest --
    # data = X_individuals_rest
    # task_conds = list(info_individuals_rest.names)

    # half1_inds = np.array([x for x in np.arange(len(task_conds)) if info_individuals_rest.iloc[x].run == 1])
    # half2_inds = np.setdiff1d(np.arange(len(task_conds)), half1_inds)



    X_individuals_half_1 = data[:, half1_inds, :]
    X_individuals_half_2 = data[:, half2_inds, :]

    data = np.array([X_individuals_half_1, X_individuals_half_2])
    data = data.transpose([1, 0, 2, 3])


    # fill nans with 0
    data[np.isnan(data)] = 0

    criterion = 'global'
    variances = ds.decompose_pattern_into_group_indiv_noise(data, criterion=criterion)
    output = dict()
    output[criterion] = variances
    print(f'Task variance all task sessions\nGroup: {output[criterion][0][0]:.2f}\nSubject: {output[criterion][0][1]:.2f}\nError: {output[criterion][0][2]:.2f}')

    

    # -- MDTB Variance decomposition --
    mean_centering = True
    type='CondHalf'
    part_column='half'
    mdtb_dataset = DataSetMDTB(data_dir)
    # --- Resting-state ---
    data_rest, info_rest = mdtb_dataset.get_data(ses_id='ses-rest', type='Net69Run', space='MNISymC2', subj=subject_subset)
    data_reshaped_rest = reshape_data(data_rest, info_rest, cond_column='net_id', mean_centering=mean_centering)
    vars_rest = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_rest, criterion='global')
    print(f'Rest variance\nGroup: {vars_rest[0][0]:.2f}\nSubject: {vars_rest[0][1]:.2f}\nError: {vars_rest[0][2]:.2f}')

    # --- Task ---
    data_task, info_task = mdtb_dataset.get_data(ses_id='ses-s1', type='CondHalf', space='MNISymC2', subj=subject_subset)
    data_reshaped_task = reshape_data(data_task, info_task, cond_column=mdtb_dataset.cond_ind, part_column=part_column, mean_centering=mean_centering)
    vars_task = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_task, criterion='global')
    print(f'Task variance\nGroup: {vars_task[0][0]:.2f}\nSubject: {vars_task[0][1]:.2f}\nError: {vars_task[0][2]:.2f}')
    
    # --- Task Neocortex ---
    data_task_neocortex, info_task_neocortex = mdtb_dataset.get_data(ses_id='ses-s1', type='CondHalf', space='fs32k')
    data_reshaped_task_neocortex = reshape_data(data_task_neocortex, info_task, cond_column=mdtb_dataset.cond_ind, part_column=part_column, mean_centering=mean_centering)
    data_reshaped_task_neocortex = reshape_data(data_task_neocortex, info_task, cond_column=mdtb_dataset.cond_ind, part_column=part_column, mean_centering=False)
    
    vars_task_neo = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_task_neocortex, criterion='global')
    print(f'Neocortical task variance\nGroup: {vars_task_neo[0][0]:.2f}\nSubject: {vars_task_neo[0][1]:.2f}\nError: {vars_task_neo[0][2]:.2f}')

    # Create dataset all runs are data from the first run for each subject (to maximize subject variance)
    # data_max_subject_variance = np.zeros_like(data_reshaped_task)
    # for i in range(data_reshaped_task.shape[1]):
    #     data_max_subject_variance[:, i, :, :] = data_reshaped_task[:, 0, :, :]
    # vars_subject = ds.decompose_pattern_into_group_indiv_noise(data_max_subject_variance, criterion='global')
    # print(f'Same runs variance\nGroup: {vars_subject[0][0]:.2f}\nSubject: {vars_subject[0][1]:.2f}\nError: {vars_subject[0][2]:.2f}')

    # # Create dataset where all 24 subjects are data from the first subject (to maximize group variance)
    # data_max_group_variance = np.zeros_like(data_reshaped_task)
    # for i in range(data_reshaped_task.shape[0]):
    #     data_max_group_variance[i, :, :, :] = data_reshaped_task[0, :, :, :]
    # vars_group = ds.decompose_pattern_into_group_indiv_noise(data_max_group_variance, criterion='global')
    # print(f'Same subjects variance\nGroup: {vars_group[0][0]:.2f}\nSubject: {vars_group[0][1]:.2f}\nError: {vars_group[0][2]:.2f}')

    # # Create dataset where error variance is maximized
    # data_max_error_variance = np.random.rand(*data_reshaped_task.shape)
    # vars_error = ds.decompose_pattern_into_group_indiv_noise(data_max_error_variance, criterion='global')
    # print(f'Random numbers variance\nGroup: {vars_error[0][0]:.2f}\nSubject: {vars_error[0][1]:.2f}\nError: {vars_error[0][2]:.2f}')


    
    # -- Pontine Variance decomposition --
    mean_centering = False
    type='CondHalf'
    part_column='half'
    pontine = DataSetPontine(data_dir)

    # --- Task ---
    data_task, info_task = pontine.get_data(type='CondHalf', space='MNISymC2')
    data_reshaped_task = reshape_data(data_task, info_task, cond_column=pontine.cond_ind, part_column=part_column, mean_centering=mean_centering)
    vars_task = ds.decompose_pattern_into_group_indiv_noise(data_reshaped_task, criterion='global')
    print(f'Task variance\nGroup: {vars_task[0][0]:.2f}\nSubject: {vars_task[0][1]:.2f}\nError: {vars_task[0][2]:.2f}')

    pass
    
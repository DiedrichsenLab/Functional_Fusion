import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import os, pickle


def setProjectPath():
    """
    setProjectPath specifies the main directory of the project and the results folder

    usage
    projectPath, mainResultsPath = setProjectPath()

    last modified: 2023.11.06
    """
    pythonAnalysisPath = os.getcwd()
    projectPath = os.path.join(pythonAnalysisPath, '..')    # exploration path
    mainResultsPath = os.path.join(projectPath, 'results')

    return projectPath, mainResultsPath


def decompose_pattern_into_group_indv_noise(data, criterion='global'):
    """
    this function decompose a collection of (across subjects and partitions) activity patterns (N condition x P voxels)
    into group, individual and noise components, return the variance estimates of each component

    input
        data
            * either a list of numpy ndarrays [sub-01: (n_partitions_01 x n_conditions x n_voxels), sub-02: (n_partitions_02 x n_conditions x n_voxels), ...]
            * or an ndarray of shape n_subjects x n_partitions x n_conditions x n_voxels, i.e., S x R x N x P
        criterion (str)
            * 'voxel_wise':     partition variance components for each voxel separately
            * 'condition_wise': partition variance components for each condition separately
            * 'global':         partition variance components for the whole pattern (N x P)
    output
        v_g, v_s, v_e (variance for group, subject, and noise)

    last modified: 2023.11.17
    """

    if isinstance(data, list):
        n_subjects = len(data)
        n_conditions = data[0].shape[1]
        n_voxels = data[0].shape[2]
        n_partitions_each_subject = [x.shape[0] for x in data]
        n_partitions = np.max(n_partitions_each_subject)

        # X = np.full((n_subjects, n_partitions, n_conditions, n_voxels), np.nan)
        X = np.full((n_subjects, n_partitions, n_conditions, n_voxels), 0)

        for subjI in np.arange(n_subjects):
            X[subjI, 0:n_partitions_each_subject[subjI]] = data[subjI]

    else:
        [n_subjects, n_partitions, n_conditions, n_voxels] = data.shape
        X = data

    if criterion == 'voxel_wise':
        Y = X.transpose([3, 0, 1, 2])
    elif criterion == 'condition_wise':
        Y = X.transpose([2, 0, 1, 3])
    elif criterion == 'global':
        Y = X.reshape((1, n_subjects, n_partitions, n_conditions * n_voxels))
    else:
        Y = np.empty(1)
        print('invalid criterion')

    [n_reps, _, _, n_features] = Y.shape
    variances = np.zeros((n_reps, 3))

    for repI in np.arange(n_reps):
        Z = np.zeros((n_subjects * n_partitions, n_features))
        indx = 0
        indx_mat_1 = np.zeros((Z.shape[0], Z.shape[0]))     # for across subjects
        indx_mat_2 = np.zeros(indx_mat_1.shape)             # for within subject, across partitions
        indx_mat_3 = np.zeros(indx_mat_1.shape)             # for within-observations

        for subjI in np.arange(n_subjects):
            subj_start_indx = subjI * n_partitions
            subj_end_indx = (subjI + 1) * n_partitions
            indx_mat_2[subj_start_indx: subj_end_indx, subj_start_indx: subj_end_indx] = 1
            other_subj_inds = np.setdiff1d(np.arange(Z.shape[0]), np.arange(subj_start_indx, subj_end_indx))
            indx_mat_1[subj_start_indx: subj_end_indx, other_subj_inds] = 1
            for partI in np.arange(n_partitions):
                Z[indx] = Y[repI, subjI, partI]
                indx_mat_3[indx, indx] = 1
                indx = indx + 1

        variance_mat = np.matmul(Z, Z.T)

        # indices for across subjects (i.e., v_g)
        indx_mat_1 = np.triu(indx_mat_1, 1)
        where_across_subject = np.where(indx_mat_1 > 0)
        var_1 = np.nanmean(variance_mat[where_across_subject])

        # indices for within subject across partitions (i.e., v_g + v_s)
        indx_mat_2 = np.triu(indx_mat_2, 1)
        where_within_subject_across_partitions = np.where(indx_mat_2 > 0)
        var_2 = np.nanmean(variance_mat[where_within_subject_across_partitions])

        # indices for within-observation pairs (i.e., v_g + v_s + v_e)
        indx_mat_3 = np.triu(indx_mat_3, 0)
        where_within_observation = np.where(indx_mat_3 > 0)
        var_3 = np.nanmean(variance_mat[where_within_observation])

        v_e = var_3 - var_2
        v_s = var_2 - var_1
        v_g = var_1
        variances[repI] = np.array([v_g, v_s, v_e]) / (v_g + v_s + v_e)

    return variances


if __name__ == '__main__':
    # n_subjects = 4
    # n_partitions = 5
    # n_conditions = 10
    # n_voxels = 100
    # data = np.random.randn(n_subjects, n_partitions, n_conditions, n_voxels)

    projectPath, mainResultsPath = setProjectPath()
    # base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion'
    surface_helpers_dir = os.path.join(projectPath, 'surface_helpers')

    resultsPath = os.path.join(mainResultsPath, os.path.basename(__file__).replace('.py', ''))
    if not os.path.exists(resultsPath):
        os.makedirs(resultsPath)

    output = dict()
    PKL_output = os.path.join(resultsPath, 'output_whole_cortex.pkl')

    X_individuals, info_individuals, dataset_obj_individuals = ds.get_dataset(base_dir,
                                                                              dataset='MDTB',
                                                                              atlas='fs32k',
                                                                              subj=None,
                                                                              sess='all',
                                                                              type='CondHalf')
    task_conds = list(info_individuals.names)
    half1_inds = np.array([x for x in np.arange(len(task_conds)) if task_conds[x].__contains__('half1')])
    half2_inds = np.setdiff1d(np.arange(len(task_conds)), half1_inds)

    X_individuals_half_1 = X_individuals[:, half1_inds, :]
    X_individuals_half_2 = X_individuals[:, half2_inds, :]

    data = np.array([X_individuals_half_1, X_individuals_half_2])
    data = data.transpose([1, 0, 2, 3])

    # fill nans with 0
    data[np.isnan(data)] = 0

    criterion = 'global'
    variances = decompose_pattern_into_group_indv_noise(data, criterion=criterion)
    output[criterion] = variances

    criterion = 'voxel_wise'
    variances = decompose_pattern_into_group_indv_noise(data, criterion=criterion)
    output[criterion] = variances

    criterion = 'condition_wise'
    variances = decompose_pattern_into_group_indv_noise(data, criterion=criterion)
    output[criterion] = variances

    with open(PKL_output, 'wb') as pk:
        pickle.dump(output, pk)


"""
Reliability module for establishing pattern-reliability on fMRI data

In general, this module provides functions to calculate the reliability of a fMRI dataset across different conditions, sessions, or subjects.
The data is assumed to be in the form of (n_subjects x n_measures x n_voxels) or simple in (n_measures x n_voxels) matrix.
"""

import numpy as np
import Functional_Fusion.matrix as matrix
import Functional_Fusion.util as util
from numpy import eye,sqrt,nansum
from numpy.linalg import inv

def within_subj(data, cond_vec, part_vec,
                separate = 'none',
                subtract_mean=True):
    """ Calculates the within-subject reliability of a data set
    Data (X) is grouped by condition vector, and the
    partition vector indicates the independent measurements
    Does the calculation for each subject if X is a 3d array
    Args:
        data (ndarray): (num_subj x) num_trials x num_voxel tensor of data
        cond_vec (ndarray): num_trials condition vector
        part_vec (ndarray): num_trials partition vector
        separate (str): {'none','voxel_wise','condition_wise'}
        subtract_mean (bool): Remove the mean per voxel in partition across conditions?
    Returns:
        r (ndarray): (n_subjects x) n_separate array of correlations
    """
    partitions = np.unique(part_vec)
    n_part = partitions.shape[0]
    if data.ndim == 2:
        X = data.reshape(1, data.shape[0], data.shape[1])
        single_subj = True
    else:
        single_subj = False
        X = data
    n_subj = X.shape[0]

        # Transform the data into a 4d array
    X = flat2ndarray(X, part_vec, cond_vec)
    [n_subjects, n_part,n_cond, n_voxels] = X.shape
    # remove mean for each voxel for each partition
    if subtract_mean:
        X = X - np.nanmean(X, axis=2, keepdims=True)

    # rearrange the data to be in the form of (n_separate, n_subjects,n_partitions, n_features)
    if separate == 'voxel_wise':
        Y = X.transpose([0, 3,  1, 2])
    elif separate == 'condition_wise':
        Y = X.transpose([0, 2, 1, 3])
    elif separate in ['none']:
        Y = X.reshape((n_subjects, 1, n_part, n_cond * n_voxels))
    else:
        raise(NameError('separate needs to be none, voxel_wise, or condition_wise'))
    [n_split, _, _, n_features] = Y.shape

    # This computes the sums of squares-matrix for each split separately (broadcasting)
    YY = np.matmul(Y,Y.transpose([0,1,3,2]))
    # Get the mean of ondiagonal and offdiagonal elements
    ondiag = np.where(np.eye(n_part))
    offdiag = np.where(1-np.eye(n_part))

    # Average the within-partition and across-partition cross-products for each subject and separate split
    SS_1 = np.nanmean(YY[:,:,ondiag[0],ondiag[1]],axis=2)
    SS_2 = np.nanmean(YY[:,:,offdiag[0],offdiag[1]],axis=2)

    # Compute the variances from the sums of squares
    if subtract_mean and separate == 'voxel_wise':
        n_df = (n_cond - 1)
    elif subtract_mean and separate == 'cond_wise':
        n_df = n_voxels
    elif subtract_mean:
        n_df = (n_cond-1) * n_voxels
    else:
        n_df = n_features
    v_e = (SS_1 - SS_2) / n_df
    v_s = (SS_2) / n_df
    r = v_s / (v_s + v_e)
    if single_subj:
        r = r[0]
    return r

def between_subj(data, cond_vec=None,
                             separate='none',
                             subtract_mean=True):
    """ Calculates the average between-subject reliability
    The data is averaged across multiple measurements
    first - so the reliability is for the mean patterns

    Args:
        datas (ndarray): num_subj x num_trials x num_voxel tensor of data
        cond_vec (ndarray): num_trials condition vector (otherwise assumed to be identity)
        separate (str): {'none','voxel_wise','condition_wise'}
        subtract_mean (bool): Remove the mean per voxel before correlation calc?

    Returns:
        r (ndarray): num_subj vector of correlations
    """
    n_subj,n_trials,n_voxels = data.shape
    if cond_vec is not None:
        Z = matrix.indicator(cond_vec)
    else:
        Z = eye(n_trials)
    n_cond = Z.shape[1]
    X = np.zeros((n_subj,n_cond,n_voxels))
    for s in range(n_subj):
        X[s,:,:] = util.nan_linear_model(Z, data[s,:,:])
    if subtract_mean:
        X=X-np.nanmean(X,axis=1,keepdims=True)

    # rearrange the data to be in the form of (n_separate, n_subjects,n_partitions, n_features)
    if separate == 'voxel_wise':
        Y = X.transpose([2, 0, 1])
    elif separate == 'condition_wise':
        Y = X.transpose([1, 0, 2])
    elif separate in ['none']:
        Y = X.reshape((1,n_subj, n_cond * n_voxels))
    else:
        raise(NameError('separate needs to be none, voxel_wise, or condition_wise'))

    # This computes the sums of squares-matrix for each split separately (broadcasting)
    YY = np.matmul(Y,Y.transpose([0,2,1]))
    # Get the mean of ondiagonal and offdiagonal elements
    ondiag = np.where(np.eye(n_subj))
    offdiag = np.where(1-np.eye(n_subj))

    # Average the within-partition and across-partition cross-products for each subject and separate split
    SS_1 = np.nanmean(YY[:,ondiag[0],ondiag[1]],axis=1)
    SS_2 = np.nanmean(YY[:,offdiag[0],offdiag[1]],axis=1)

    v_e = (SS_1 - SS_2)
    v_s = (SS_2)
    r = v_s / (v_s + v_e)
    return r

def within_subj_loo(data, cond_vec,
                          part_vec, 
                          separate='voxel_wise',
                          subtract_mean=True):
    """ Calculates the within-subject reliability of a data set
    Data (X) is grouped by condition vector, and the
    partition vector indicates the independent measurements
    Does the calculation for each subejct if X is a 3d array
    Args:
        X (ndarray): (num_subj x) num_trials x num_voxel tensor of data
        cond_vec (ndarray): num_trials condition vector
        part_vec (ndarray): num_trials partition vector
        separate (str): {'none','voxel_wise','condition_wise'}
        subtract_mean (bool): Remove the mean per voxel before correlation calc?
    Returns:
        r (ndarray): (num_subj x) num_partition matrix of correlations
    """
    partitions = np.unique(part_vec)
    n_part = partitions.shape[0]
    if len(data.shape) == 2:
        X = data.reshape(1, data.shape[0], data.shape[1])
        single_subj = True
    else:
        single_subj = False
        X= data.copy()
    
    n_subj = X.shape[0]
    if separate=='voxel_wise':
        r = np.zeros((n_subj, n_part, data.shape[2]))
    elif separate=='condition_wise':
        raise(NameError('condition_wise not implemented yet'))
    elif separate=='none':
        r = np.zeros((n_subj, n_part))
    else:
        raise(NameError('separate needs to be none, voxel_wise, or condition_wise'))
    Z = matrix.indicator(cond_vec)
    for s in np.arange(n_subj):
        for pn, part in enumerate(partitions):
            i1 = part_vec == part
            i2 = part_vec != part
            X1 = util.nan_linear_model(Z[i1, :], X[s, i1, :])
            X2 = util.nan_linear_model(Z[i2, :], X[s, i2, :])
            # Check if this partition contains nan row
            if subtract_mean:
                X1 -= np.nanmean(X1, axis=0)
                X2 -= np.nanmean(X2, axis=0)
            if separate=='voxel_wise':
                r[s, pn, :] = nansum(X1 * X2, axis=0) / \
                    sqrt(nansum(X1 * X1, axis=0)
                         * nansum(X2 * X2, axis=0))
            else:
                r[s, pn] = nansum(X1 * X2) / \
                    sqrt(nansum(X1 * X1) * nansum(X2 * X2))
    if single_subj:
        r = r[0, :]    
    return r

def between_subj_loo(data, cond_vec=None,
                            separate='none',
                            subtract_mean=True):
    """ Calculates the correlation of the responses of each of the subjects with the mean of the other subjects. 
    If cond_vec is given, the data is averaged across multiple measurem
    first.

    Args:
        data (ndarray): num_subj x num_trials x num_voxel tensor of data
        cond_vec (ndarray): num_trials condition vector
        separate (str): {'none','voxel_wise','condition_wise'}
        subtract_mean (bool): Remove the mean per voxel before correlation calc?

    Returns:
        r (ndarray): num_subj vector of correlations 
    """
    n_subj = data.shape[0]
    n_trials = data.shape[1]
    if cond_vec is not None:
        Z = matrix.indicator(cond_vec)
    else:
        Z = eye(n_trials)
    subj_vec = np.arange(n_subj)
    if separate == 'voxel_wise':
        r = np.zeros((n_subj, data.shape[2]))
    elif separate=='condition_wise':
        raise(NameError('condition_wise not implemented yet'))
    elif separate=='none':
        r = np.zeros((n_subj,))
    else:
        raise(NameError('separate needs to be none, voxel_wise, or condition_wise'))
    for s, i in enumerate(subj_vec):
        X1 = util.nan_linear_model(Z, data[s, :, :])
        i2 = subj_vec != s
        X2 = util.nan_linear_model(Z, np.nanmean(data[i2, :, :], axis=0))
        if subtract_mean:
            X1 -= np.nanmean(X1, axis=0)
            X2 -= np.nanmean(X2, axis=0)
        if separate=='voxel_wise':
            r[i, :] = nansum(X1 * X2, axis=0) / \
                sqrt(nansum(X1 * X1, axis=0)
                     * nansum(X2 * X2, axis=0))
        else:
            r[i] = nansum(X1 * X2) / sqrt(nansum(X1 * X1) * nansum(X2 * X2))
    return r

def decompose_subj_group(data, cond_vec, part_vec,
                         separate='none',
                         subtract_mean=True):
    """
    this function decompose a collection of (across subjects and partitions) activity patterns (N condition x P voxels)
    into group, individual and noise components, returns the variance estimates of each component.

    Args:
        data (ndarray):
            n_subjects x n_trials x n_voxels array
        cond_vec (ndarray):
            n_trials condition vector
        part_vec (ndarray):
            n_trials partition vector
        separate (str):
            * 'none':         partition variance components for the whole pattern (N x P) -> returns a single row
            * 'voxel_wise':     partition variance components for each voxel separately -> returns as many rows as voxels
            * 'condition_wise':     partition variance components for each condition separately -> returns as many rows as conditions
            * 'subject_wise':   partition variance components for the whole pattern (NxP) -> but return split by Subjects
    Returns:
        variances: (K x 3 ndarray): v_g, v_s, v_e (variance for group, subject, and noise), where K is the number of voxels, conditions, subjects, or 1
    """
    if data.ndim != 3:
        raise(NameError('data needs to be a 3d array'))

    # Transform the data into a 4d array
    X = flat2ndarray(data, part_vec, cond_vec)
    [n_subjects, n_part,n_cond, n_voxels] = X.shape
    # remove mean for each voxel for each partition
    if subtract_mean:
        X = X - np.nanmean(X, axis=2, keepdims=True)

    # rearrange the data to be in the form of (n_separate, n_subjects,n_partitions, n_features)
    if separate == 'voxel_wise':
        Y = X.transpose([3, 0, 1, 2])
    elif separate == 'condition_wise':
        Y = X.transpose([2, 0, 1, 3])
    elif separate in ['global','none','subject_wise']:
        Y = X.reshape((1, n_subjects, n_part, n_cond * n_voxels))
    else:
        raise(NameError('separate needs to be none, voxel_wise, condition_wise, or subject_wise'))
    [n_split, _, _, n_features] = Y.shape

    # reshape the data to be in the form of (n_split, n_subjects*n_partitions, n_features)
    Y = Y.reshape((n_split, n_subjects * n_part, n_features))
    subj_vec = np.kron(np.arange(n_subjects), np.ones(n_part))
    part_vec = np.kron(np.ones(n_subjects), np.arange(n_part))
    N = n_subjects * n_part
    same_subj = np.equal(subj_vec.reshape(N,1),subj_vec.reshape(1,N))
    same_part = np.equal(part_vec.reshape(N,1),part_vec.reshape(1,N))

    # This computes the sums of squares-matrix for each split separately (broadcasting)
    YY = np.matmul(Y,Y.transpose((0,2,1)))

    # Find cross-products of same subject and same partition
    if separate == 'subject_wise':
        SS_1 = np.zeros((n_subjects,))
        SS_2 = np.zeros((n_subjects,))
        SS_3 = np.zeros((n_subjects,))
        for s in np.arange(n_subjects):
            # Get the rows pertaining to the subject
            YYs = YY[:,subj_vec==s,:]
            same_subj_s = same_subj[subj_vec==s,:]
            same_part_s = same_part[subj_vec==s,:]
            SS_1[s] = np.nanmean(YYs[:,~same_subj_s],axis=1)
            SS_2[s] = np.nanmean(YYs[:,same_subj_s & ~same_part_s],axis=1)
            SS_3[s] = np.nanmean(YYs[:,same_subj_s & same_part_s],axis=1)
    else:
        SS_1 = np.nanmean(YY[:,~same_subj],axis=1)
        SS_2 = np.nanmean(YY[:,same_subj & ~same_part],axis=1)
        SS_3 = np.nanmean(YY[:,same_subj & same_part],axis=1)

    # Compute the variances from the sums of squares
    if subtract_mean and separate == 'voxel_wise':
        n_df = (n_cond - 1)
    elif subtract_mean and separate == 'cond_wise':
        n_df = n_voxels
    elif subtract_mean:
        n_df = (n_cond-1) * n_voxels
    else:
        n_df = n_features
    v_e = (SS_3 - SS_2) / n_df
    v_s = (SS_2 - SS_1) / n_df
    v_g = SS_1 / n_df
    variances = np.c_[v_g, v_s, v_e]

    return variances


def flat2ndarray(flat_data, part_vec, cond_vec):
    """
    convert flat data (n_subjects x n_trials x n_voxels) into a 4d ndarray (n_subjects x n_partitions x n_conditions x n_voxels). Also works on
    (n_trials x n_voxels) array that is converted into a
    (n_partitions x n_conditions x n_voxels) array.

    Args:
        flat_data (nd-array): (n_subjects x ) n_trials x n_voxels array
        part_vec: (nd-array): n_trials -vector for partition index
        cond_vec: (nd-array): n_trials -vector for condition index

    Returns:
        data: (nd-array): (n_subjects x ) n_partitions x n_conditions x n_voxels
    """

    [n_subjects, n_trials, n_voxels] = flat_data.shape

    unique_partitions = np.unique(part_vec)
    n_partitions = unique_partitions.size

    unique_conditions = np.unique(cond_vec)
    n_conditions = unique_conditions.size

    data = np.zeros((n_subjects, n_partitions, n_conditions, n_voxels))

    for p,part in enumerate(unique_partitions):
        for c,cond in enumerate(unique_conditions):
            trial_inds = np.where((cond_vec == cond) & (part_vec == part))[0]
            data[:, p, c, :] = np.mean(flat_data[:, trial_inds, :], axis=1)
    return data
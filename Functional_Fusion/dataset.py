#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data fusion project dataset module

The class for converting and mapping raw data from multi-dataset
to a standard data structure that can be used in Diedrichsen lab

"""
from pathlib import Path
import numpy as np
import pandas as pd
import os, time, sys
import os.path as op

import Functional_Fusion.util as util
import Functional_Fusion.matrix as matrix
import Functional_Fusion.atlas_map as am
import scipy.linalg as sl
import nibabel as nb
from numpy import eye, zeros, ones, empty, nansum, sqrt
from numpy.linalg import pinv, solve
import warnings
import glob
import re


def get_dataset_class(base_dir, dataset):
    """ Returns the dataset object without loading the data

    Args:
        base_dir (str): Functional Fusion base directory
        dataset (str): _description_

    Returns:
        my_dataset (Dataset): Dataset object
    """
    T = pd.read_csv(base_dir + '/dataset_description.tsv', sep='\t')
    T.name = [n.casefold() for n in T.name]
    i = np.where(dataset.casefold() == T.name)[0]
    if len(i) == 0:
        raise (NameError(f'Unknown dataset: {dataset}'))
    dsclass = getattr(sys.modules[__name__], T.class_name[int(i)])
    dir_name = T.dir_name[int(i)]
    if dir_name[0] == '/':
        abs_path = dir_name
    elif dir_name[0] == '.':  # Relative path relative to fusion project
        abs_path = Path(base_dir).parent + Path(dir_name[1:])
    else:
        abs_path = base_dir + '/' + T.dir_name[int(i)]
    my_dataset = dsclass(abs_path)
    return my_dataset

def get_dataset(base_dir, dataset, atlas='SUIT3', sess='all', subj=None,
                type=None, smooth=None, info_only=False):
    """get_dataset tensor and data set object

    Args:
        base_dir (str): Basis directory for the Functional Fusion datastructure
        dataset (str): Data set indicator
        atlas (str): Atlas indicator. Defaults to 'SUIT3'.
        sess (str or list): Sessions. Defaults to 'all'.
        subj (ndarray, str, or list):  Subject numbers /names to get [None = all]
        type (str): 'CondHalf','CondRun', etc....
    Returns:
        data (nd.array):nsubj x ncond x nvox data tensor
        info (pd.DataFrame): Dataframe with info about the data
        my_dataset (DataSet): Dataset object
    """

    my_dataset = get_dataset_class(base_dir, dataset)

    # Get defaults sessions from dataset itself
    if sess == 'all':
        sess = my_dataset.sessions
    elif not isinstance(sess, (list, np.ndarray)):
        sess = [sess]
    # Empty default type (Future change: define per file?)
    if type is None:
        type = my_dataset.default_type

    # Load all data and concatenate
    # across sessions
    info_l = []
    data_l = []
    for s in sess:
        dat, inf = my_dataset.get_data(atlas, s, type, subj, smooth=smooth)
        data_l.append(dat)
        inf['sess'] = [s] * inf.shape[0]
        info_l.append(inf)

    info = pd.concat(info_l, ignore_index=True, sort=False)
    data = np.concatenate(data_l, axis=1)
    return data, info, my_dataset


def build_dataset_from_fusionProject(dataset, atlas, base_dir, sess='all',
                                     cond_ind='cond_num_uni', type='CondHalf',
                                     part_ind='half', subj=None, join_sess=False,
                                     join_sess_part=False, smooth=None):
    """ Builds dataset, cond_vec, part_vec, subj_ind from the given
        dataset in Functional fusion project

    Args:
        dataset (str): Names of the dataset to build
        atlas (object): Atlas indicator
        sess (list): list of 'all' or list of sessions
        design_ind (list, optional): _description_. Defaults to None.
        part_ind (list, optional): _description_. Defaults to None.
        subj (list, optional): _description_. Defaults to None.
        join_sess (bool, optional): Model the sessions with a single model.
            Defaults to True.
    Returns:
        data, cond_vec, part_vec, subj_ind
    """
    sub = 0
    data, cond_vec, part_vec, subj_ind = [],[],[],[]
    # Run over datasets get data + design
    dat, info, tds = get_dataset(base_dir, dataset, atlas=atlas.name,
                                 sess=sess, type=type, smooth=smooth)

    # Sub-index the subjects: (REMOVE AND PASS TO GET_DATA)
    if subj is not None:
        dat = dat[subj, :, :]
    n_subj = dat.shape[0]
    # Find correct indices
    if cond_ind is None:
        cond_ind = tds.cond_ind
    if part_ind is None:
        part_ind = tds.part_ind

    # Make different sessions either the same or different
    if join_sess:
        data.append(dat)
        cond_vec.append(info[cond_ind].values.reshape(-1, ))

        # Check if we want to set no partition after join sessions
        if join_sess_part:
            part_vec.append(np.ones(info[part_ind].shape))
        else:
            part_vec.append(info[part_ind].values.reshape(-1, ))
        subj_ind.append(np.arange(sub, sub + n_subj))
    else:
        if sess == 'all':
            sessions = tds.sessions
        else:
            sessions = sess
        # Now build and split across the correct sessions:
        for s in sessions:
            indx = info.sess == s
            data.append(dat[:, indx, :])
            cond_vec.append(info[cond_ind].values[indx].reshape(-1, ))
            part_vec.append(info[part_ind].values[indx].reshape(-1, ))
            subj_ind.append(np.arange(sub, sub + n_subj))

    return data, cond_vec, part_vec, subj_ind


def prewhiten_data(data):
    """ prewhitens a list of data matrices.
    It assumes that the last row of each data matrix is the ResMS-value
    Returns a list of data matrices that is one shorter

    Args:
        data (list of ndarrays): List of data arrays
    """
    # Get the resms and prewhiten the data
    data_n = []
    for i in range(len(data)):
        # Prewhiten the data univariately
        resms = data[i][-1, :]
        data_n.append(data[i][0:-1, :])
        resms[resms <= 0] = np.nan
        data_n[i] = data_n[i] / np.sqrt(resms)
    return data_n


def agg_data(info, by, over, subset=None):
    """ Aggregates data over rows (condition) safely by sorting them by the fields in "by"
    while integrating out "over".
    Adds a n_rep field to count how many instances are of each
    Returns condensed data frame + Contrast matrix.

    Args:
        info (DataFrame): Original DataFrame
        by (list): Fields that define the index of the new data
        over (list): Fields to ignore / integrate over. All other fields
            will be pulled through.
        subset (bool array): If given, ignores certain rows from the
            original data frame
    Return
        data_info (DataFrame): Reduced data frame
        C (ndarray): Indicator matrix defining the mapping from full to reduced
    Example: 
        data,info,mdtb= ds.get_data('MDTB','MNISymDentate1',ses_id='ses-s1',type='CondRun')   
        cinfo,C = ds.agg_data(info,['cond_num_uni'],['run','half','reg_num','names'])
        cdata = np.linalg.pinv(C) @ data
    """
    # Subset original data frame as needed
    info_n = info.copy()
    if subset is None:
        indx = np.arange(info.shape[0])
    else:
        indx = np.nonzero(subset.values)[0]
        info_n = info_n[subset]

    # Generate n_rep field if not present
    if 'n_rep' not in info.columns:
        info_n['n_rep'] = np.ones((info_n.shape[0],))

    # Other contains the fields to remain constant
    other = list(info.columns.values)
    for ov in over + by:
        other.remove(ov)

    # Define operations on data
    operations = {'n_rep': np.sum}
    for o in other:
        operations[o] = max

    # Group the new data frame
    info_gb = info_n.groupby(by)
    data_info = info_gb.agg(operations).reset_index()

    # Build indicator matrix for averaging
    C = np.zeros((info.shape[0], data_info.shape[0]))
    for i, (k, v) in enumerate(info_gb.indices.items()):
        C[indx[v], i] = 1
    return data_info, C


def agg_parcels(data, label_vec, fcn=np.nanmean):
    """ Aggregates data over colums to condense to parcels

    Args:
        data (ndarray): Either 2d or 3d data structure
        labels (ndarray): 1d-array that gives the labels
        fcn (function): Function to use to aggregate over these
    Returns:
        aggdata (ndarray): Aggregated either 2d or 3d data structure
        labels (ndarray): Region number corresponding to each "column"
    """
    # Subset original data frame as needed
    labels = np.unique(label_vec[label_vec > 0])
    n_parcels = len(labels)
    psize = np.array(data.shape)
    psize[-1] = n_parcels
    parcel_data = np.zeros(psize)
    for i, l in enumerate(labels):
        parcel_data[..., i] = fcn(
            data[..., label_vec == l], axis=len(psize) - 1)
    return parcel_data, labels

def combine_parcel_labels(labels_org,labelvec_org,labels_new):
    """ Combines parcel labels from a new atlas to an existing atlas
    Example call: 
    combine_parcel_labels(labels_org,labelvec_org,['A.L','A.R','S..','M3.')
    To get different aggregations of the Nettekoven atlas
    * A.L includes A1-4L 
    * A.R includes A1-4R
    * S.. includes all S-areas
    * M3. includes M3L and M3R 

    Args:
        labels_org (list of str): Original label names (should include '0' for first)
        labelvec_org (ndarray): Original label vector 
        labels_new (list of str): List of regexpressions for new labels
    Returns:
        labelvec_new (ndarray): New label vector

    """
    labelvec_new = np.zeros(labelvec_org.shape)
    for i,l in enumerate(labels_new):
        for j,lo in enumerate(labels_org):
            if re.match(l,lo):
                labelvec_new[labelvec_org== j] = i
    
    return labelvec_new

def optimal_contrast(data, C, X, reg_in=None, baseline=None):
    """Recombines betas from a GLM into an optimal new contrast, taking into account a design matrix

    Args:
        data (list of ndarrays): List of N x P_i arrays of data
        C (ndarray): N x Q array indicating contrasts
        X (ndarray): Optimal design matrix - Defaults to None.
        reg_in (ndarray): Contrast of interest: Logical vector indicating
            which rows of C we will put in the matrix
        baseline (ndarray): Fixed effects contrast removed after estimation
    """
    # Check the sizes
    N, Q = C.shape
    T, Np = X.shape
    # infer number of regressors of no interest that are not in the data structure
    num_nointerest = Np - N
    # Add to the contrast matrix
    Cn = sl.block_diag(C, np.eye(num_nointerest))
    # Make new design matrix
    Xn = X @ Cn
    # Loop over the data:
    data_new = []
    for i in range(len(data)):
        # Append the regressors of no interest regressors
        dat = np.concatenate([data[i],
                              np.zeros((num_nointerest, data[i].shape[1]))])
        # Do the averaging / reweighting:
        d = solve(Xn.T @ Xn, Xn.T @ X @ dat)
        # Subset to the contrast of interest
        if reg_in is not None:
            d = d[reg_in, :]
        # Now subtract baseline
        d = remove_baseline(d,baseline)
        # Put the data in a list:
        data_new.append(d)
    return data_new

def remove_baseline(data, baseline):
    """ Removes a baseline from the data"""
    if baseline is None:
        return data 
    Q = data.shape[0]
    R = eye(Q) - baseline @ pinv(baseline)
    return R @ data

def reliability_within_subj(X, part_vec, cond_vec,
                            voxel_wise=False,
                            subtract_mean=True):
    """ Calculates the within-subject reliability of a data set
    Data (X) is grouped by condition vector, and the
    partition vector indicates the independent measurements

    Args:
        X (ndarray): num_subj x num_trials x num_voxel tensor of data
        part_vec (ndarray): num_trials partition vector
        cond_vec (ndarray): num_trials condition vector
        voxel_wise (bool): Return the results as map or overall?
        subtract_mean (bool): Remove the mean per voxel before correlation calc?
    Returns:
        r (ndarray)L: num_subj x num_partition matrix of correlations
    """
    partitions = np.unique(part_vec)
    n_part = partitions.shape[0]
    n_subj = X.shape[0]
    if voxel_wise:
        r = np.zeros((n_subj, n_part, X.shape[2]))
    else:
        r = np.zeros((n_subj, n_part))
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
            if voxel_wise:
                r[s, pn, :] = nansum(X1 * X2, axis=0) / \
                    sqrt(nansum(X1 * X1, axis=0)
                         * nansum(X2 * X2, axis=0))
            else:
                r[s, pn] = nansum(X1 * X2) / \
                    sqrt(nansum(X1 * X1) * nansum(X2 * X2))
    return r


def reliability_between_subj(X, cond_vec=None,
                             voxel_wise=False,
                             subtract_mean=True):
    """ Calculates the correlation of the responses of each of the subjects with the mean of the other subjects.
    If cond_vec is given, the data is averaged across multiple measurem
    first.

    Args:
        X (ndarray): num_subj x num_trials x num_voxel tensor of data
        part_vec (ndarray): num_trials partition vector
        voxel_wise (bool): Return the results as map or overall?
        subtract_mean (bool): Remove the mean per voxel before correlation calc?

    Returns:
        r (ndarray): num_subj vector of correlations
    """
    n_subj = X.shape[0]
    n_trials = X.shape[1]
    if cond_vec is not None:
        Z = matrix.indicator(cond_vec)
    else:
        Z = eye(n_trials)
    subj_vec = np.arange(n_subj)
    if voxel_wise:
        r = np.zeros((n_subj, X.shape[2]))
    else:
        r = np.zeros((n_subj,))
    for s, i in enumerate(subj_vec):
        X1 = util.nan_linear_model(Z, X[s, :, :])
        i2 = subj_vec != s
        X2 = util.nan_linear_model(Z, np.nanmean(X[i2, :, :], axis=0))
        if subtract_mean:
            X1 -= np.nanmean(X1, axis=0)
            X2 -= np.nanmean(X2, axis=0)
        if voxel_wise:
            r[i, :] = nansum(X1 * X2, axis=0) / \
                sqrt(nansum(X1 * X1, axis=0)
                     * nansum(X2 * X2, axis=0))
        else:
            r[i] = nansum(X1 * X2) / sqrt(nansum(X1 * X1) * nansum(X2 * X2))
    return r

def reliability_maps(base_dir, dataset_name, atlas='MNISymC3', type='CondHalf',
                     subtract_mean=True, voxel_wise=True, subject_wise=False):
    """    Calculates the average within subject reliability maps across sessions for a single dataset

    Args:
        base_dir (str / path): Base directory
        dataset_name (str): Name of data set
        atlas (str): _description_. Defaults to 'MNISymC3'.
        subtract_mean (bool): Remove the mean per voxel before correlation calc?

    Returns:
        _type_: _description_
    """
    data, info, dataset = get_dataset(base_dir, dataset_name, atlas=atlas, type = type)
    n_sess = len(dataset.sessions)
    n_vox = data.shape[2]
    Rel = np.zeros((n_sess, n_vox))
    if subject_wise:
        Rel = np.zeros((n_sess, data.shape[0], n_vox))
    for i, s in enumerate(dataset.sessions):
        indx = info.sess == s
        r = reliability_within_subj(data[:, indx, :],
                                    part_vec=info[dataset.part_ind][indx],
                                    cond_vec=info[dataset.cond_ind][indx],
                                    voxel_wise=voxel_wise,
                                    subtract_mean=subtract_mean)
        if subject_wise:
            Rel[i, :, :] = np.nanmean(r, axis=1)
        else:
            Rel[i, :] = np.nanmean(np.nanmean(r, axis=0), axis=0)
    return Rel, dataset.sessions

def decompose_pattern_into_group_indiv_noise(data, criterion='global'):
    """
    this function decompose a collection of (across subjects and partitions) activity patterns (N condition x P voxels)
    into group, individual and noise components, returns the variance estimates of each component.

    Args:
        data (list,ndarray):
            * either a list of numpy ndarrays [sub-01: (n_partitions_01 x n_conditions x n_voxels), sub-02: (n_partitions_02 x n_conditions x n_voxels), ...]
            * or an ndarray of shape n_subjects x n_partitions x n_conditions x n_voxels, i.e., S x R x N x P
        criterion (str):
            * 'global':         partition variance components for the whole pattern (N x P) -> returns a single row
            * 'voxel_wise':     partition variance components for each voxel separately -> returns as many rows as voxels
            * 'condition_wise':     partition variance components for each condition separately -> returns as many rows as conditions
            * 'subject_wise':   partition variance components for the whole pattern (NxP) -> but return split by Subjects
    Returns:
        variances: (K x 3 ndarray): v_g, v_s, v_e (variance for group, subject, and noise), where K is the number of voxels, conditions, subjects, or 1
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

    # rearrange the data to be in the form of (n_split, n_subjects,n_partitions, n_features)
    if criterion == 'voxel_wise':
        Y = X.transpose([3, 0, 1, 2])
    elif criterion == 'condition_wise':
        Y = X.transpose([2, 0, 1, 3])
    elif criterion in ['global','subject_wise']:
        Y = X.reshape((1, n_subjects, n_partitions, n_conditions * n_voxels))
    else:
        Y = np.empty(1)
        print('invalid criterion')
    [n_split, _, _, n_features] = Y.shape

    # reshape the data to be in the form of (n_split, n_subjects*n_partitions, n_features)
    Y = Y.reshape((n_split, n_subjects * n_partitions, n_features))
    subj_vec = np.kron(np.arange(n_subjects), np.ones(n_partitions))
    part_vec = np.kron(np.ones(n_subjects), np.arange(n_partitions))
    N = n_subjects * n_partitions
    same_subj = np.equal(subj_vec.reshape(N,1),subj_vec.reshape(1,N))
    same_part = np.equal(part_vec.reshape(N,1),part_vec.reshape(1,N))

    # This computes the sums of squares-matrix for each split separately (broadcasting)
    YY = np.matmul(Y,Y.transpose((0,2,1)))

    # Find cross-products of same subject and same partition
    if criterion == 'subject_wise':
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
    v_e = (SS_3 - SS_2) / n_features
    v_s = (SS_2 - SS_1) / n_features
    v_g = SS_1 / n_features
    variances = np.c_[v_g, v_s, v_e]

    return variances


def flat2ndarray(flat_data, part_vec, cond_vec):
    """
    convert flat data (n_subjects x n_trials x n_voxels) into a 4d ndarray (n_subjects x n_partitions x n_conditions x n_voxels)

    Args:
        flat_data (nd-array): n_subjects x n_obs x n_voxels array
        part_vec: (nd-array): n_obs -vector for partition index
        cond_vec: (nd-array): n_obs -vector for condition index

    Returns:
        data: (nd-array): n_subjects x n_partitions x n_conditions x n_voxels
    """

    [n_subjects, n_obs, n_voxels] = flat_data.shape

    unique_partitions = np.unique(part_vec)
    n_partitions = unique_partitions.size

    unique_conditions = np.unique(cond_vec)
    n_conditions = unique_conditions.size

    data = np.zeros((n_subjects, n_partitions, n_conditions, n_voxels))

    for partI in np.arange(n_partitions):
        for condI in np.arange(n_conditions):
            trial_inds = np.where(cond_vec == unique_conditions[condI] and part_vec == unique_partitions[partI])
            data[:, partI, condI, :] = np.mean(flat_data[:, trial_inds, :], axis=1)

class DataSet:
    def __init__(self, base_dir):
        """DataSet class:
        Implements the interface for each of the data set
        Note that the actual preprocessing and glm estimate
        do not have to be performed with functionality provided by
        this class. The class is just a instrument to present the user with
        a uniform interface of how to get subject info

        Args:
            basedir (str): basis directory
        """
        self.base_dir = base_dir
        self.surface_dir = base_dir + '/derivatives/{0}/anat'
        self.anatomical_dir = base_dir + '/derivatives/{0}/anat'
        self.estimates_dir = base_dir + '/derivatives/{0}/estimates'
        self.func_dir = base_dir + '/derivatives/{0}/func'
        self.suit_dir = base_dir + '/derivatives/{0}/suit'
        self.data_dir = base_dir + '/derivatives/{0}/data'
        # assume that the common atlas directory is on the level before
        self.atlas_dir = os.path.join(os.path.dirname(base_dir), 'Atlases')
        # Some information that a standard data set should have
        self.sessions = [None]
        self.default_type = None
        self.cond_ind = None  # Condition Indicator (field in tsv file )
        self.part_ind = None  # Partition Indicator (field in tsv file )
        self.cond_name = None  # Condition Names (field in tsv file )

    def get_participants(self):
        """ returns a data frame with all participants
        available in the study. The fields in the data frame correspond to the
        standard columns in participant.tsv.
        https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html

        Returns:
            Pinfo (pandas data frame): participant information in standard bids format
        """
        self.part_info = pd.read_csv(
            self.base_dir + '/participants.tsv', delimiter='\t')
        return self.part_info

    def get_data_fnames(self, participant_id, session_id=None, type='Cond'):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
            type (str): Type of data. Defaults to 'Cond' for task-based data. For rest data use 'Tseries'.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'

        if type[:4] == 'Cond' or type[:4] == 'Task':
            T = pd.read_csv(
                dirw + f'/{participant_id}_{session_id}_reginfo.tsv', sep='\t')
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
            fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        elif type == 'Tseries' or type == 'FixTseries':
            # Find all run files of the structure f'{dirw}/{participant_id}_{session_id}_run-??.nii'
            fnames = glob.glob(f'{dirw}/{participant_id}_{session_id}_run-??.nii')
            runs = [int(fname.split('run-')[-1].split('_')[0].split('.')[0]) for fname in fnames]
            runs = np.unique(runs)
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{r:02}.nii' for r in runs]
            if type == 'FixTseries':
                # Make sure to load fix-cleaned timeseries
                fnames = [f'{dirw}/{participant_id}_{session_id}_run-{r:02}_fix.nii' for r in runs]
            try:
                # Load timeseries info file if it exists
                T = pd.read_csv(
                    dirw + f'/{participant_id}_{session_id}_tinfo.tsv', sep='\t')
            except:
                # Create timeseries info yourself if it doesn't exist
                timepoints = [nb.load(fname).shape[-1] for fname in fnames]
                runs = np.repeat(runs, timepoints)
                # Timepoints start counting at 1, not 0!
                timepoints = np.arange(sum(timepoints))+1
                timepoints_string = [f'T{timepoint:04d}' for timepoint in timepoints]
                T = pd.DataFrame({'run': runs,
                                    'timepoint': timepoints_string,
                                    'time_id':timepoints}, index=None)
                
        return fnames, T

    def get_info(self, ses_id='ses-s1', type=None, subj=None, fields=None):
        """Loads all the CIFTI files in the data directory of a certain space / type and returns they content as a Numpy array

        Args:
            space (str): Atlas space (Defaults to 'SUIT3').
            ses_id (str): Session ID (Defaults to 'ses-s1').
            type (str): Type of data (Defaults to 'CondHalf').
            subj (ndarray): Subject numbers to get - by default all
            fields (list): Column names of info stucture that are returned
                these are also be tested to be equivalent across subjects
        Returns:
            Data (ndarray): (n_subj, n_contrast, n_voxel) array of data
            info (DataFramw): Data frame with common descriptor
        """
        T = self.get_participants()

        # only get data from subjects that have rest, if specified in dataset description
        if type == 'Tseries' and 'ses-rest' in T.columns:
                subj = T[T['ses-rest'] == 1].participant_id.tolist()

        # Deal with subset of subject option
        if subj is None:
            subj = np.arange(T.shape[0])
        else:
            subj = [T.participant_id.tolist().index(i) for i in subj]

        if type is None:
            type = self.default_type

        max = 0
        # Loop over the different subjects to find the most complete info
        for s in T.participant_id.iloc[subj]:
            # Get an check the information
            info_raw = pd.read_csv(self.data_dir.format(s)
                                   + f'/{s}_{ses_id}_{type}.tsv', sep='\t')
            # Reduce tsv file when fields are given
            if fields is not None:
                info = info_raw[fields]
            else:
                info = info_raw

            # Keep the most complete info
            if info.shape[0] > max:
                info_com = info
                max = info.shape[0]
        return info_com

    def get_atlasmaps(self, atlas, sub, ses_id, smooth=None):
        """This function generates atlas map for the data of a specific subject into a specific atlas space. The general DataSet.get_atlasmaps defines atlas maps for different spaces
            - SUIT: Using individual normalization from source space. 
            - MNI152NLin2009cSymC: Via indivual SUIT normalization + group
            - MNI152NLin6AsymC: Via indivual SUIT normalization + group
            - MNI152Lin2009cSym: Via individual MNI normalization
            - MNI152NLin6Asym: Via individual MNI normalization
        fs32k: Via individual pial and white surfaces (need to be in source space)
        Other dataset classes will overwrite and extend this function.

        Args:
            atlas (FunctionFusion.Atlas):
                Functional Fusion atlas object
            sub (str):
                Subject_id for the individual subject
            ses_id (str):
                Session_id for the individual subject if atlasmap is session dependent. (defaults to none)
            smooth (float):
                Width of smoothing kernel for extraction. Defaults to None.
        Returns:
            AtlasMap:
                Built AtlasMap object
        """
        atlas_maps = []
        adir = self.anatomical_dir.format(sub)
        edir = self.estimates_dir.format(sub)
        if atlas.space == 'SUIT':
            deform = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.space in ['MNI152NLin2009cSymC','MNI152NLin6AsymC']:
            # This is nornmalization over SUIT->MNI (cerebellum only)
            deform1  = am.get_deform(atlas.space, 'SUIT')
            deform2 = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            deform = [deform1, deform2]
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.space in ['MNI152NLin2009cSym']:
            # This is direct MNI normalization 
            deform = adir + f'/{sub}_space-{atlas.space}_xfm.nii'
            mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.space == 'fs32k':
            for i, struc in enumerate(atlas.structure):
                if struc=='cortex_left':
                    hem = 'L'
                elif struc=='cortex_right':
                    hem = 'R'
                else:
                    raise ValueError('Structure for fs32k needs to be cortex_left or cortex_right.')
                pial = adir + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
                mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(atlas.vertex[i],
                                                  white, pial, mask))
                atlas_maps[i].build()
        else:
            raise ValueError(f'Atlas space {atlas.space} not supported for extraction')
        return atlas_maps

    def extract_all(self,
                    ses_id='ses-s1',
                    type='CondHalf',
                    atlas='SUIT3',
                    smooth=2.0,
                    subj='all'):
        """Extracts data in Volumetric space from a dataset in which the data is stored in Native space. Saves the results as CIFTI files in the data directory.

        Args:
            ses_id (str):
                Session. Defaults to 'ses-s1'.
            type (str):
                Type for condense_data. Defaults to 'CondHalf'.
            atlas (str):
                Short atlas string. Defaults to 'SUIT3'.
            smooth (float):
                Smoothing kernel. Defaults to 2.0.
            subj (list / str):
                List of Subject numbers to get use. Default = 'all'
        """
        myatlas, _ = am.get_atlas(atlas)
        # create and calculate the atlas map for each participant
        T = self.get_participants()
        if subj != 'all':
            T = T.iloc[subj]
        for s in T.participant_id:
            print(f'Atlasmap {s}')
            atlas_maps = self.get_atlasmaps(myatlas, s, ses_id,
                                                  smooth=smooth)
            print(f'Extract {s}')
            fnames, info = self.get_data_fnames(s, ses_id, type=type)
            data = am.get_data_nifti(fnames, atlas_maps)
            data, info = self.condense_data(data, info, type,
                                            participant_id=s, ses_id=ses_id)
            # Write out data as CIFTI file
            C = myatlas.data_to_cifti(data, info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_{type}.tsv', sep='\t', index=False)


    def get_data(self, space='SUIT3', ses_id='ses-s1', type=None,
                 subj=None, fields=None, smooth=None, verbose=False):
        """Loads all the CIFTI files in the data directory of a certain space / type and returns they content as a Numpy array

        Args:
            space (str): Atlas space (Defaults to 'SUIT3').
            ses_id (str): Session ID (Defaults to 'ses-s1').
            type (str): Type of data (Defaults to 'CondHalf').
            subj (ndarray, str, or list):  Subject numbers /names to get [None = all]
            fields (list): Column names of info stucture that are returned
                these are also be tested to be equivalent across subjects
        Returns:
            Data (ndarray): (n_subj, n_contrast, n_voxel) array of data
            info (DataFramw): Data frame with common descriptor
        """
        T = self.get_participants()
        # Assemble the data
        Data = None
        # Deal with subset of subject option
        if subj is None:
            subj = T.participant_id
            # only get data from subjects that have rest, if specified in dataset description
            if type == 'Tseries' and 'ses-rest' in T.columns:
                subj = T[T['ses-rest'] == 1].participant_id.tolist()
        elif isinstance(subj, str):
            subj = [subj]
        elif isinstance(subj, (int,np.integer)):
            subj = [T.participant_id.iloc[subj]]
        elif isinstance(subj, (list, np.ndarray)):
            if isinstance(subj[0], (int,np.integer)):
                subj = T.participant_id.iloc[subj]
            elif isinstance(subj[0], str):
                subj = subj
            else:
                raise (NameError('subj must be a list of strings or integers'))
        else:
            raise (NameError('subj must be a str, int, list or ndarray'))
        if type is None:
            type = self.default_type

        info_com = self.get_info(
            subj=subj, ses_id=ses_id, type=type, fields=fields)

        # Loop again to assemble the data
        Data_list = []
        for i, s in enumerate(subj):
            # If you add verbose printout, make it so
            # that by default it is suppressed by a verbose=False option
            if verbose:
                print(f'- Getting data for {s} in {space}')
            # Load the data
            if smooth is not None:
                C = nb.load(self.data_dir.format(s)
                            + f'/{s}_space-{space}_{ses_id}_{type}_desc-sm{int(smooth)}.dscalar.nii')
            else:
                C = nb.load(self.data_dir.format(s)
                            + f'/{s}_space-{space}_{ses_id}_{type}.dscalar.nii')
            this_data = C.get_fdata()

            # Check if this subject data in incomplete
            if this_data.shape[0] != info_com.shape[0]:
                this_info = pd.read_csv(self.data_dir.format(s)
                                        + f'/{s}_{ses_id}_{type}.tsv', sep='\t')
                base = np.asarray(info_com['names'])
                incomplete = np.asarray(this_info['names'])
                for j in range(base.shape[0]):
                    if base[j] != incomplete[j]:
                        warnings.warn(
                            f'{s}, {ses_id}, {type} - missing data {base[j]}')
                        incomplete = np.insert(
                            np.asarray(incomplete), j, np.nan)
                        this_data = np.insert(this_data, j, np.nan, axis=0)
            Data_list.append(this_data[np.newaxis, ...])
        # concatenate along the first dimension (subjects)
        Data = np.concatenate(Data_list, axis=0)
        # Ensure that infinite values (from div / 0) show up as NaNs
        Data[np.isinf(Data)] = np.nan
        return Data, info_com

    def group_average_data(self, ses_id=None,
                               type=None,
                               atlas='SUIT3',
                               subj=None):
            """Loads group data in SUIT space from a standard experiment structure
            averaged across all subjects. Saves the results as CIFTI files in the data/group directory.

            Args:
                ses_id (str, optional): Session ID. If not provided, the first session ID in the dataset will be used.
                type (str, optional): Type of data. If not provided, the default type will be used.
                atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
                subj (list or None, optional): Subset of subjects to include in the group average. If None, all subjects will be included.

            """
            if ses_id is None:
                ses_id = self.sessions[0]
            if type is None:
                type = self.default_type

            data, info = self.get_data(space=atlas, ses_id=ses_id,
                                       type=type, subj=subj)
            # average across participants
            X = np.nanmean(data, axis=0)
            # make output cifti
            s = self.get_participants().participant_id[0]
            C = nb.load(self.data_dir.format(s) +
                        f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            C = nb.Cifti2Image(dataobj=X, header=C.header)
            # save output
            dest_dir = op.join(self.data_dir.format('group'))
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/group_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            if 'sn' in info.columns:
                info = info.drop(columns=['sn'])

            info.to_csv(dest_dir +
                        f'/group_{ses_id}_{type}.tsv', sep='\t', index=False)

class DataSetNative(DataSet):
    """Data set with estimates data stored as
    nifti-files in Native space.
    """

    def get_atlasmaps(self, atlas, sub, ses_id, smooth=None):
        """This function generates atlas map for the data of a specific subject into a specific atlas space.
        For Native space, we are using indivdual maps for SUIT and surface space.
        Addtiionally, we defines deformations MNI space via the individual normalization into MNI152NLin6Asym (FSL, SPM Segement).
        Other MNI space (symmetric etc) are not implemented yet.
        Args:
            atlas (FunctionFusion.Atlas):
                Functional Fusion atlas object
            sub (str):
                Subject_id for the individual subject
            ses_id (str):
                Session_id for the individual subject if atlasmap is session dependent. (defaults to none)
            smooth (float):
                Width of smoothing kernel for extraction. Defaults to None.
        Returns:
            AtlasMap: Built AtlasMap object
        """
        atlas_maps = []

        if atlas.space == ['MNI152NLin2009cSym','MNI152NLin6Asym']:
            # This is for MNI standard space)
            deform = self.anatomical_dir.format(sub) + f'/{sub}_space-{atlas.space}_xfm.nii'
            if not os.path.exists(deform):
                warnings.warn(f'No individual deformation found for {atlas.space} in {sub} - resortinh to MNI_xfm.nii')
                deform = self.anatomical_dir.format(sub) + f'/{sub}_space-MNI_xfm.nii'
            edir = self.estimates_dir.format(sub)
            mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        else:
            atlas_maps = super().get_atlasmaps(atlas,sub,ses_id,smooth=smooth)
        return atlas_maps

class DataSetMNIVol(DataSet):
    def __init__(self, base_dir,space='MNI152NLin6Asym'):
        """Data set with estimates data stored as nifti-files in a standard group space. The exact MNI template should be indicated in the space-argument ('MNI152NLin6Asym','MNI152N2009cAsym','MNI152N2009cSym'). The small deformations between the different MNI spaces are implemented when extracting the data.

        Args:
            basedir (str): basis directory
            space (str): Group Space in which data is stored (Defaults to 'MNI152NLin6Asym').
        """
        super().__init__(base_dir)
        self.space = space

    def get_atlasmaps(self, atlas, sub, ses_id, smooth=None):
        """This function generates atlas map for the data stored in MNI space.
        For SUIT and surface space, it goes over deformations estimated on the individual anatomy. If atlas.space matches dataset.space, it uses no deformation, but a direct readout. For mismatching MNI space it tries to find the correct transformation file.
        Args:
            atlas (FunctionFusion.Atlas):
                Functional Fusion atlas object
            sub (str):
                Subject_id for the individual subject
            ses_id (str):
                Session_id for the individual subject if atlasmap is session dependent. (defaults to none)
            smooth (float):
                Width of smoothing kernel for extraction. Defaults to None.
        Returns:
            AtlasMap: Built AtlasMap object
        """
        atlas_spaces = ['MNI152NLin6Asym','MNI152NLin2009cAsym','MNI152NLin2009cSym']
        atlas_maps = []
        edir = self.estimates_dir.format(sub)
        mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
        # Matching MNI space
        if atlas.space == self.space:
            atlas_maps.append(am.AtlasMapDeform(atlas.world, None, mask))
            atlas_maps[0].build(smooth=smooth)
        # Mis-matching MNI space
        elif atlas.space in atlas_spaces:
            deform = am.get_deform(atlas.space, self.space)
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        # Any other space (SUIT + fs32k)
        else:
            atlas_maps = super().get_atlasmaps(atlas,sub,ses_id,smooth=smooth)
        return atlas_maps

class DataSetCifti(DataSet):
    """Data set that comes in HCP-format in already pre-extracted cifti files.
    """

    def get_data_fnames(self, participant_id, session_id=None):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T = pd.read_csv(
            dirw + f'/{participant_id}_{session_id}_reginfo.tsv', sep='\t')
        fnames = [f'{dirw}/{participant_id}_{session_id}_beta.dscalar.nii']
        fnames.append(
            f'{dirw}/{participant_id}_{session_id}_resms.dscalar.nii')
        return fnames, T

    def extract_all(self, ses_id='ses-s1', type='CondHalf', atlas='SUIT3'):
        """Extracts cerebellar data. Saves the results as CIFTI files in the data directory.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        myatlas, _ = am.get_atlas(atlas)
        # Get the correct map into CIFTI-format
        if isinstance(myatlas, am.AtlasVolumetric):
            deform = am.get_deform(myatlas.space,'MNI152NLin6Asym')
            atlas_map = am.AtlasMapDeform(myatlas.world,
                                          deform, None)
            atlas_map.build(smooth=2.0)
        elif isinstance(myatlas, am.AtlasSurface):
            atlas_map = myatlas
        # Extract the data for each participant
        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}')
            fnames, info = self.get_data_fnames(s, ses_id)
            data = am.get_data_cifti(fnames, [atlas_map])

            data, info = self.condense_data(data, info, type,
                                            participant_id=s, ses_id=ses_id)
            C = myatlas.data_to_cifti(data[0], info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_{type}.tsv', sep='\t', index=False)


class DataSetMDTB(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-s1', 'ses-s2']
        self.default_type = 'CondHalf'
        self.cond_ind = 'cond_num_uni'
        self.part_ind = 'half'
        self.cond_name = 'cond_name'

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Condense the data in a certain way optimally
        'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
        'CondRun': Conditions with seperate estimates per run. Defaults to 'CondHalf'.

        Args:
            data (list): List of extracted datasets
            info (DataFrame): Data Frame with description of data - row-wise
            type (str): Type of extraction:
            participant_id (str): ID of participant
            ses_id (str): Name of session

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provided
        """

        # Depending on the type, make a new contrast
        info['half'] = 2 - (info.run < 9)
        if type == 'Tseries' or type == 'FixTseries':
            info['names'] = info['timepoint']
            data_new, data_info = data, info

        else:
            if type == 'CondHalf':
                data_info, C = agg_data(info,
                                        ['half', 'cond_num'],
                                        ['run', 'reg_num'],
                                        subset=(info.instruction == 0))
                data_info['names'] = [
                    f'{d.cond_name.strip()}-half{d.half}' for i, d in data_info.iterrows()]
                # Baseline substraction
                B = matrix.indicator(data_info.half, positive=True)

            elif type == 'CondRun':
                data_info, C = agg_data(info,
                                        ['run', 'cond_num'],
                                        [],
                                        subset=(info.instruction == 0))
                data_info['names'] = [
                    f'{d.cond_name}-run{d.run:02d}' for i, d in data_info.iterrows()]

                # Baseline substraction
                B = matrix.indicator(data_info.run, positive=True)
            elif type == 'CondAll':

                data_info, C = agg_data(info,
                                        ['cond_num'],
                                        ['run', 'half', 'reg_num'],
                                        subset=(info.instruction == 0))
                data_info['names'] = [
                    f'{d.cond_name}' for i, d in data_info.iterrows()]

                # Baseline substraction
                B = np.ones((data_info.shape[0],1))

            # Prewhiten the data
            data_n = prewhiten_data(data)

            # Load the designmatrix and perform optimal contrast
            dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
            X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
            reg_in = np.arange(C.shape[1], dtype=int)
            #  contrast for all instructions
            CI = matrix.indicator(info.run * info.instruction, positive=True)
            C = np.c_[C, CI]
            data_new = optimal_contrast(data_n, C, X, reg_in, baseline=B)

        return data_new, data_info


class DataSetHcpResting(DataSetCifti):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-rest1', 'ses-rest2']
        self.hem_name = ['cortex_left', 'cortex_right']
        self.default_type = 'Net67Run'
        self.cond_ind = 'net_id'
        self.cond_name = 'names'
        self.part_ind = 'half'

    def get_data_fnames(self, participant_id, ses_id):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
        Returns:
            fnames (list): List of fnames
        """

        dirw = self.func_dir.format(participant_id)
        fnames = []
        if ses_id == "ses-rest1":
            runs = np.arange(0, 2)
        elif ses_id == "ses-rest2":
            runs = np.arange(2, 4)
        # idx = self.sessions.index(ses_id)
        T = pd.read_csv(
            dirw + f'/{participant_id}_{ses_id}_reginfo.tsv', sep='\t')
        for r in runs:
            fnames.append(
                f'{dirw}/sub-{participant_id}_run-{r}_space-MSMSulc.dtseries.nii')
        return fnames, T

    def condense_data(self, data, info, type, participant_id=None, ses_id=None):
        if type == 'Tseries' or type == 'FixTseries':
            info['names'] = info['timepoint']
        return data, info

    def regress_networks(self, X, Y):
        """Regresses a spatial map (X) into data (Y).
        Returns the network timecourses.

        Args:
            X (np.arry): 4D Network data of the signal components
                (default input networks are in fs32k Space: 59518 vertices x nComponents )
            Y (<nibabel CIFTI image object>): fMRI timeseries in volume
                Has to be in the same space as networks (59518 vertices x nTimepoints )
        Returns:
            network_timecourse (np.ndarray):
                A numpy array (nTimepoints x nNetworks) with the fMRI timecourse for
                each resting-state network
        """
        X = X.T
        Y = Y.T.squeeze()
        d_excluded = np.where(np.isnan(Y))[0].shape[0]
        v_excluded = np.unique(np.where(np.isnan(Y))[0]).shape[0]
        print(
            f'Setting nan datapoints ({v_excluded} unique vertices) to zero. Entire timeseries: {d_excluded/v_excluded == Y.shape[1]}')
        Y[np.isnan(Y)] = 0
        network_timecourse = np.matmul(np.linalg.pinv(X), Y)

        return network_timecourse

    def average_within_Icos(self, label_file, data, atlas="fs32k"):
        """Average the raw time course for voxels within a parcel

        Args:
            label_file (str): cortical parcellation label file
            Y (np.ndarray): fMRI timeseries in volume
                Has to be in the same space as networks (59518 vertices x nTimepoints)
        Returns:
            A numpy array (nNetworks x nTimepoints) with the fMRI timecourse for
            each resting-state network
        """

        # create an instance of atlas to get the label vector
        atlas, ainfo = am.get_atlas(atlas)

        # create label_vector by passing on the label file
        # Set unite_struct to true if you want to integrate over left and right hemi
        atlas.get_parcel(label_file, unite_struct=False)

        # use agg_parcel to aggregate data over parcels and get the list of unique parcels
        parcel_data, parcels = agg_parcels(
            data, atlas.label_vector, fcn=np.nanmean)

        # fill nan value in Y to zero
        print("Setting nan datapoints (%d unique vertices) to zero"
              % np.unique(np.where(np.isnan(parcel_data))[1]).shape[0])
        # Y = np.nan_to_num(np.transpose(Y))
        parcel_data = np.nan_to_num(parcel_data)

        # return np.matmul(np.linalg.pinv(indicator_mat.T), Y)
        return parcel_data, parcels

    def correlate(self, X, Y):
        """ Correlate X and Y numpy arrays after standardizing them"""
        X = util.zstandarize_ts(X)
        Y = util.zstandarize_ts(Y)
        return Y.T @ X / X.shape[0]

    def connectivity_fingerprint(self, source, target, info, type):
        """ Calculate the connectivity fingerprint of a target region

        Args:
            source (np.ndarray): Source data
            target (np.nzdarray): Target timecourse
            info (pandas.DataFrame): Information about the source data
            type (str): Type of fingerprint to calculate ('Run' or 'All').
                        Estimates fingerprint from each run seperately or from all concatenated runs.

        Returns:
            coef (np.ndarray): Connectivity fingerprint
        """
        coefs = []
        if type == 'Run':
            for run in info.run.unique():
                data_run = source[info.run == run]
                net_run = target.T[info.run == run]
                coef = self.correlate(data_run, net_run)
                coefs.append(coef)

        elif type == 'All':
            coef = self.correlate(source, target)
            coefs.append(coef)

        return np.vstack(coefs)


class DataSetPontine(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-01']
        self.default_type = 'TaskHalf'
        self.cond_ind = 'task_num'
        self.cond_name = 'task_name'
        self.part_ind = 'half'

    def condense_data(self, data, info,
                      type='TaskHalf',
                      participant_id=None,
                      ses_id=None):
        """ Condense the data from the pontine project after extraction

        Args:
            data (list of ndarray)
            info (dataframe)
            type (str): Type of extraction:
                'TaskHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'TaskRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.
            participant_id (str): ID of participant
            ses_id (str): Name of session

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """

        # Depending on the type, make a new contrast
        info['half'] = 2 - (info.run < 9)
        n_cond = np.max(info.reg_id)

        if type == 'TaskHalf':
            data_info, C = agg_data(info,
                                    ['half', 'reg_id'],
                                    ['run', 'reg_num'],
                                    subset=(info.reg_id > 0))
            data_info['names'] = [
                f'{d.task_name.strip()}-half{d.half}' for i, d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.half, positive=True)

        elif type == 'TaskRun':

            data_info, C = agg_data(info,
                                    ['run', 'reg_id'],
                                    ['reg_num'],
                                    subset=(info.reg_id > 0))
            data_info['names'] = [
                f'{d.task_name.strip()}-run{d.run}' for i, d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.half, positive=True)

        elif type == 'TaskAll':
            data_info, C = agg_data(info,
                                    ['reg_id'],
                                    ['run', 'half', 'reg_num'],
                                    subset=(info.reg_id > 0))
            data_info['names'] = [
                f'{d.task_name.strip()}' for i, d in data_info.iterrows()]
            # Baseline substraction
            B = np.ones((data_info.shape[0],1))

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(self.estimates_dir.format(participant_id) + f'/{ses_id}/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1], dtype=int)
        CI = matrix.indicator(info.run * info.instruction, positive=True)
        C = np.c_[C, CI]

        data_new = optimal_contrast(data_n, C, X, reg_in)

        return data_new, data_info


class DataSetNishi(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-01', 'ses-02']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'task_name'
        self.part_ind = 'half'

    def condense_data(self, data, info,
                      type='TaskHalf',
                      participant_id=None,
                      ses_id=None):
        """ Condense the data from the pontine project after extraction

        Args:
            data (list of ndarray)
            info (dataframe)
            type (str): Type of extraction:
                'TaskHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'TaskRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.
            participant_id (str): ID of participant
            ses_id (str): Name of session

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        # Depending on the type, make a new contrast
        info['half'] = 2 - (info.run < (len(np.unique(info.run)) / 2 + 1))
        n_cond = np.max(info.reg_id)

        if type == 'CondHalf':
            data_info, C = agg_data(info,
                                    ['half', 'reg_id'],
                                    ['run', 'reg_num'])
            data_info['names'] = [
                f'{d.task_name.strip()}-half{d.half}' for i, d in data_info.iterrows()]

            # Baseline substraction
            B = matrix.indicator(data_info.half, positive=True)

        elif type == 'CondRun':
            data_info, C = agg_data(info,
                                    ['run', 'reg_id'],
                                    ['reg_num'])

            data_info['names'] = [
                f'{d.task_name.strip()}-run{d.run:02d}' for i, d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.run, positive=True)
        elif type == 'CondAll':
            data_info, C = agg_data(info,
                                    ['reg_id'],
                                    ['run', 'half'])
            # Baseline substraction
            B = np.ones((data_info.shape[0],))

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(self.estimates_dir.format(participant_id) + f'/{ses_id}/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1], dtype=int)
        data_new = optimal_contrast(data_n, C, X, reg_in, baseline=B)

        return data_new, data_info


class DataSetIBC(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-archi',
                         'ses-clips4',
                         'ses-enumeration',
                         'ses-hcp1', 'ses-hcp2',
                         'ses-lyon1', 'ses-lyon2',
                         'ses-mathlang',
                         'ses-mtt1', 'ses-mtt2',
                         'ses-preference',
                         'ses-rsvplanguage',
                         'ses-spatialnavigation',
                         'ses-tom']
        self.default_type = 'CondHalf'
        self.cond_ind = 'cond_num_uni'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def get_participants(self):
        """ returns a data frame with all participants complete participants
        Returns:
            Pinfo (pandas data frame): participant information in standard bids format
        """
        self.part_info = pd.read_csv(
            self.base_dir + '/participants.tsv', delimiter='\t')
        return self.part_info[self.part_info.complete == 1]

    def get_data_fnames(self, participant_id, session_id=None, type = "CondHalf"):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T = pd.read_csv(
            dirw + f'/{participant_id}_{session_id}_reginfo.tsv', sep='\t')
        fnames = [
            f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

    def get_atlasmaps(self, atlas, sub, ses_id, smooth=None):
        """This function generates atlas map for the data of a specific subject into a specific atlas space.
        Uses the general ones, but overwrites the choice of masks
        Args:
            atlas (FunctionFusion.Atlas):
                Functional Fusion atlas object
            sub (str):
                Subject_id for the individual subject
            ses_id (str):
                Session_id for the individual subject if atlasmap is session dependent. (defaults to none)
            smooth (float):
                Width of smoothing kernel for extraction. Defaults to None.
        Returns:
            AtlasMap: Built AtlasMap object
        """
        atlas_maps = []
        if atlas.space == 'SUIT':
            deform = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
            add_mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'
            atlas_map = am.AtlasMapDeform(atlas.world, deform, mask)
            atlas_map.build(smooth=2.0, additional_mask=add_mask)
        elif atlas.space in ['MNI152NLin2009cSymC','MNI152NLin6AsymC']:
            # This is nornmalization over SUIT->MNI (cerebellum only)
            deform1, m = am.get_deform(atlas.space, 'SUIT')
            deform2 = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            deform = [deform1, deform2]
            mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
            add_mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth,additional_mask=add_mask)
        else:
            atlas_maps=super().get_atlasmaps(atlas, sub, ses_id, smooth=None)
        return atlas_maps

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Condense the data in a certain way optimally
        Args:
            data (list): List of extracted datasets
            info (DataFrame): Data Frame with description of data - row-wise
            type (str): Type of extraction:
                'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.
            participant_id (str): ID of participant
            ses_id (str): Name of session

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        n_cond = np.max(info.reg_id)
        info['n_rep'] = np.ones((info.shape[0],))
        if type == 'CondHalf':
            data_info, C = agg_data(info,
                                    ['half', 'cond_num_uni'],
                                    ['run', 'reg_num'])
            data_info['names'] = [
                f'{d.cond_name.strip()}-half{d.half}' for i, d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)

        for i in range(len(data_n)):
            data_n[i] = pinv(C) @ data_n[i]
        return data_n, data_info

class DataSetDemand(DataSetCifti):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-01']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Extract data in a specific atlas space
        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        # Depending on the type, make a new contrast
        info['half'] = (info.run % 2) + 1
        n_cond = np.max(info.reg_id)
        if type == 'CondHalf':
            data_info, C = agg_data(info, ['half', 'reg_id'], ['run'])
            data_info['names'] = [f'{d.cond_name.strip()}-half{d.half}'
                                  for i, d in data_info.iterrows()]
        elif type == 'CondAll':
            data_info, C = agg_data(info, ['reg_id'], ['half', 'run'])
            data_info['names'] = [f'{d.cond_name.strip()}-half{d.half}'
                                  for i, d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Combine with contrast
        for i in range(len(data_n)):
            data_n[i] = pinv(C) @ data_n[i]
        return data_n, data_info


class DataSetWMFS(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-01', 'ses-02']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Condense the data in a certain way optimally
        Args:
            data (list): List of extracted datasets
            info (DataFrame): Data Frame with description of data - row-wise
            type (str): Type of extraction:
                'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.
            participant_id (str): ID of participant
            ses_id (str): Name of session

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """

        # Depending on the type, make a new contrast
        info['half'] = 2 - (info.run < 3)
        n_cond = np.max(info.loc[info.error == 0].reg_id)

        if type == 'CondHalf':
            # Make new data frame for the information of the new regressors

            data_info, C = agg_data(info,
                                    ['half', 'cond_num'],
                                    ['run', 'reg_num'],
                                    subset=(info.error == 0))
            data_info['names'] = [
                f'{d.cond_name.strip()}-half{d.half}' for i, d in data_info.iterrows()]
        elif type == 'CondRun':

            # Subset of info sutructure
            data_info, C = agg_data(info,
                                    ['run', 'cond_num'],
                                    ['reg_num'],
                                    subset=(info.error == 0))
            data_info['names'] = [
                f'{d.cond_name.strip()}-run{d.run:02d}' for i, d in data_info.iterrows()]
        elif type == 'CondAll':
            data_info, C = agg_data(info,
                                    ['cond_num'],
                                    ['run', 'half', 'reg_num'],
                                    subset=(info.error == 0))
            # data_info = info_gb.agg({'n_rep':np.sum}).reset_index(drop=True)
            data_info['names'] = [
                f'{d.cond_name.strip()}' for i, d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
        X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1], dtype=int)
        data_new = optimal_contrast(data_n, C, X, reg_in, baseline=None)

        return data_new, data_info


class DataSetSomatotopic(DataSetMNIVol):
    def __init__(self, dir):
        super().__init__(dir)
        self.space = 'MNI152NLin6Asym'
        self.sessions = ['ses-motor']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def get_atlasmaps(self, atlas, sub=None, ses_id = None, smooth=None):
        """ Gets group atlasmaps.
        Assumes that all scans are in the same space (self.space)

        Args:
            sub (str; optional): Subject
        Returns:
            atlas_maps (list): List of atlasmaps

        """
        # if you have group xfm file, then you get it from the atlas directory
        # if you have xfm files per subject, then you can get it from the anat dir under individual subject
        atlas_maps = []
        if atlas.structure == 'cerebellum':
            deform = self.atlas_dir + \
                f'/tpl-{self.space}/tpl-{self.space}_from-SUIT_xfm.nii'
            if atlas.name[0:4] != 'SUIT':
                deform1 = am.get_deform(atlas.space, 'SUIT')
                deform = [deform1, deform]
            mask = self.atlas_dir + \
                f'/{self.space}/{self.space}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.space == 'fs32k':
            for i, hem in enumerate(['L', 'R']):
                adir = self.anatomical_dir.format(sub)
                pial = adir + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
                mask = self.atlas_dir + \
                    f'/{self.space}/{self.space}_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(atlas.vertex[i],
                                                  white, pial, mask))
                atlas_maps[i].build()
        return atlas_maps

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Extract data in a specific atlas space
        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row

        N.B.: Because some runs are missing for session 1-3, CondRun can only be run for session 04 (which has all runs for all subjects).
        Missing runs are: S3_sess03_MOTOR6, S3_sess01_MOTOR3, S3_sess01_MOTOR4, S3_sess01_MOTOR5, S4_sess01_MOTOR6, S4_sess02_MOTOR6 & S6_sess02_MOTOR2
        """
        # Depending on the type, make a new contrast
        info['half'] = (info.run % 2) + 1
        n_cond = np.max(info.reg_id)
        if type == 'CondHalf':
            data_info, C = agg_data(info, ['half', 'reg_id'], ['run'])
            data_info['names'] = [f'{d.cond_name.strip()}-half{d.half}'
                                  for i, d in data_info.iterrows()]
        elif type == 'CondAll':
            data_info, C = agg_data(info, ['reg_id'], ['half', 'run'])
            data_info['names'] = [
                f'{d.cond_name}' for i, d in data_info.iterrows()]
        elif type == 'CondRun':
            data_info, C = agg_data(info, ['run', 'half', 'reg_id'], [])
            data_info['names'] = [f'{d.cond_name.strip()}-run{d.run}'
                                  for i, d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Combine with contrast
        for i in range(len(data_n)):
            data_n[i] = pinv(C) @ data_n[i]
        return data_n, data_info


class DataSetDmcc(DataSetMNIVol):
    def __init__(self, dir):
        super().__init__(dir)
        self.space = 'MNI152NLin2009cAsym'
        self.sessions = ['ses-axcpt-bas-mixed', 'ses-cuedts-bas-mixed', 'ses-stern-bas-mixed', 'ses-stroop-bas-mixed']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'cond_name'
        self.part_ind = 'knot_num'

    def get_data_fnames(self, participant_id, session_id=None, type='Cond'):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
            type (str): Type of data. Defaults to 'Cond' for task-based data. For rest data use 'Tseries'.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        # handle subjects with missing pro or rea sessions
        T = pd.read_csv(
            dirw + f'/{participant_id}_{session_id}_reginfo.tsv', sep='\t')
        
        if type == 'Contrast': # if you wanna look/work with at contrasts
            T = pd.read_csv(
                dirw + f'/{participant_id}_{session_id}_coninfo.tsv', sep='\t')
        
        
        if type[:4] == 'Cond' or type[:4] == 'Task' or type[:4] == 'Blck':
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
            # fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        elif type == 'Contrast':
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.con_id:02}_con.nii' for i, t in T.iterrows()]
        elif type == 'Tseries' or type == 'FixTseries':
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{r:02}.nii' for r in T.run.unique().tolist()]
        return fnames, T

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Extract data in a specific atlas space
        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row

        N.B.: Because some runs are missing for session 1-3, CondRun can only be run for session 04 (which has all runs for all subjects).
        Missing runs are: S3_sess03_MOTOR6, S3_sess01_MOTOR3, S3_sess01_MOTOR4, S3_sess01_MOTOR5, S4_sess01_MOTOR6, S4_sess02_MOTOR6 & S6_sess02_MOTOR2
        """
        # Depending on the type, make a new contrast
        info['half'] = (info.run % 2) + 1
        # n_cond = np.max(info.reg_id)
        
        if type == 'CondAll':
            data_info, C = agg_data(info, ['cond_num', 'cond_name'], ['knot_num', 'run'])
            data_info['names'] = [
                f'{d.cond_name}' for i, d in data_info.iterrows()]
        elif type == 'Contrast':
            data_info, C = agg_data(info, ['contrast_num', 'contrast_name'], ['knot_num', 'run'])
            data_info['names'] = [
                f'{d.contrast_name}' for i, d in data_info.iterrows()]
            
        

        # Prewhiten the data
        # data_n = prewhiten_data(data)
        # NOTE: I am currently using betas estimated using AFNI TentZero
        # It does not output ResMS and based on the documentation, it prewhitens the data
        # so the betas produced are already prewhitened.
        # data_n = prewhiten_data(data)
        data_n = data

        # Load the designmatrix and perform optimal contrast
        if type != 'Contrast':
            dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
            X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
            reg_in = np.arange(C.shape[1], dtype=int)
            data_new = optimal_contrast(data_n, C, X, reg_in, baseline=None)
        else:
            data_new = data_n
            for i in range(len(data_n)):
                data_new[i] = pinv(C) @ data_n[i]

        
        return data_new, data_info


class DataSetLanguage(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-01','ses-02']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'taskName'
        self.part_ind = 'half'

    def condense_data(self, data, info,
                      type='TaskHalf',
                      participant_id=None,
                      ses_id=None):
        """ Condense the data from the language localizer project after extraction

        Args:
            data (list of ndarray)
            info (dataframe)
            type (str): Type of extraction:
                'TaskHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'TaskRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.
            participant_id (str): ID of participant
            ses_id (str): Name of session

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """

        # Depending on the type, make a new contrast
        info['half'] = 2 - (info.run < 5)
        if type == 'Tseries' or type == 'FixTseries':
            info['names'] = info['timepoint']
            data_new, data_info = data, info

        else:
            if type == 'CondHalf':
                data_info, C = agg_data(info,
                                        ['half', 'reg_id'],
                                        ['run', 'reg_num'],
                                        subset=(info.reg_id >0))
                data_info['names'] = [
                    f'{d.taskName.strip()}-half{d.half}' for i, d in data_info.iterrows()]
                # Baseline substraction

            elif type == 'CondRun':

                data_info, C = agg_data(info,
                                        ['run', 'reg_id'],
                                        ['reg_num'],
                                        subset=(info.reg_id > 0))
                data_info['names'] = [
                    f'{d.taskName.strip()}-run{d.run}' for i, d in data_info.iterrows()]
                # Baseline substraction

            elif type == 'CondAll':
                data_info, C = agg_data(info,
                                        ['reg_id'],
                                        ['run', 'half', 'reg_num'],
                                        subset=(info.reg_id > 0))
                data_info['names'] = [
                    f'{d.taskName.strip()}' for i, d in data_info.iterrows()]
                # Baseline substraction

            # Prewhiten the data
            data_n = prewhiten_data(data)

            dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'

            # Load the designmatrix and perform optimal contrast
            X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
            reg_in = np.arange(C.shape[1], dtype=int)
            CI = matrix.indicator(info.run * info.inst, positive=True)
            C = np.c_[C, CI]

            data_new = optimal_contrast(data_n, C, X, reg_in)

        return data_new, data_info

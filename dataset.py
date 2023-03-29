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
import os
import os.path as op
import sys
import sys

sys.path.append("..")
import Functional_Fusion.util as util
import Functional_Fusion.matrix as matrix
import Functional_Fusion.atlas_map as am
import scipy.linalg as sl
import nibabel as nb
import nitools as nt
from numpy import eye, zeros, ones, empty, nansum, sqrt
from numpy.linalg import pinv, solve
import warnings
import SUITPy as suit
import matplotlib.pyplot as plt
import re


def get_dataset_class(base_dir, dataset):
    T = pd.read_csv(base_dir + '/dataset_description.tsv', sep='\t')
    T.name = [n.casefold() for n in T.name]
    i = np.where(dataset.casefold() == T.name)[0]
    if len(i) == 0:
        raise (NameError(f'Unknown dataset: {dataset}'))
    dsclass = getattr(sys.modules[__name__], T.class_name[int(i)])
    dir_name = base_dir + '/' + T.dir_name[int(i)]
    my_dataset = dsclass(dir_name)
    return my_dataset


def get_dataset(base_dir, dataset, atlas='SUIT3', sess='all', subj=None,
                type=None, info_only=False):
    """get_dataset
    Args:
        base_dir (str): Basis directory for the Functional Fusion repro
        dataset (str): Data set indicator
        atlas (str): Atlas indicator. Defaults to 'SUIT3'.
        sess (str or list): Sessions. Defaults to 'all'.
        type (str): 'CondHalf','CondRun', etc....
    Returns:
        _type_: _description_
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

        dat, inf = my_dataset.get_data(atlas, s, type, subj)
        data_l.append(dat)
        inf['sess'] = [s] * inf.shape[0]
        info_l.append(inf)

    info = pd.concat(info_l, ignore_index=True, sort=False)
    data = np.concatenate(data_l, axis=1)
    return data, info, my_dataset


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
        C (ndarray): Contrast matrix defining the mapping from full to reduced
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

    # Build contrast matrix for averaging
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
        if baseline is not None:
            Q = d.shape[0]
            R = eye(Q) - baseline @ pinv(baseline)
            d = R @ d
        # Put the data in a list:
        data_new.append(d)
    return data_new


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
    """ Calculates the between-subject reliability of a data set
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
        r = np.zeroes((n_subj, X.shape[2]))
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


def reliability_maps(base_dir, dataset_name, atlas='MNISymC3',
                     subtract_mean=True, voxel_wise=True):
    """    Calculates the average within subject reliability maps across sessions for a single data

    Args:
        base_dir (str / path): Base directory
        dataset_name (str): Name of data set
        atlas (str): _description_. Defaults to 'MNISymC3'.
        subtract_mean (bool): Remove the mean per voxel before correlation calc?

    Returns:
        _type_: _description_
    """
    data, info, dataset = get_dataset(base_dir, dataset_name, atlas=atlas)
    n_sess = len(dataset.sessions)
    n_vox = data.shape[2]
    Rel = np.zeros((n_sess, n_vox))
    for i, s in enumerate(dataset.sessions):
        indx = info.sess == s
        r = reliability_within_subj(data[:, indx, :],
                                    part_vec=info[dataset.part_ind][indx],
                                    cond_vec=info[dataset.cond_ind][indx],
                                    voxel_wise=voxel_wise,
                                    subtract_mean=subtract_mean)
        Rel[i, :] = np.nanmean(np.nanmean(r, axis=0), axis=0)
    return Rel, dataset.sessions


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

    def get_data(self, space='SUIT3', ses_id='ses-s1',
                 type=None, subj=None, fields=None, verbose=False):
        """Loads all the CIFTI files in the data directory of a certain space / type and returns they content as a Numpy array

        Args:
            space (str): Atlas space (Defaults to 'SUIT3').
            ses_id (str): Session ID (Defaults to 'ses-s1').
            type (str): Type of data (Defaults to 'CondHalf').
            subj (ndarray, str, or list):  Subject numbers /names to get
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
        elif isinstance(subj,str):
            subj = [subj]
        elif isinstance(subj,(list,np.ndarray)):
            if isinstance(subj[0],int):
                subj = T.participant_id.iloc[subj]
        if type is None:
            type = self.default_type

        max = 0

        # Loop over the different subjects to find the most complete info
        for i, s in enumerate(subj):

            # Get an check the information
            info_raw = pd.read_csv(self.data_dir.format(s)
                                   + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t')
            # Reduce tsv file when fields are given
            if fields is not None:
                info = info_raw[fields]
            else:
                info = info_raw

            # Keep the most complete info
            if info.shape[0] > max:
                info_com = info
                max = info.shape[0]

        # Loop again to assemble the data
        Data_list = []
        for i, s in enumerate(subj):
            # If you add verbose printout, make it so
            # that by default it is suppressed by a verbose=False option
            if verbose:
                print(f'- Getting data for {s} in {space}')
            # Load the data
            C = nb.load(self.data_dir.format(s)
                        + f'/{s}_space-{space}_{ses_id}_{type}.dscalar.nii')
            this_data = C.get_fdata()

            # Check if this subject data in incomplete
            if this_data.shape[0] != info_com.shape[0]:
                this_info = pd.read_csv(self.data_dir.format(s)
                                        + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t')
                base = np.asarray(info_com['names'])
                incomplete = np.asarray(this_info['names'])
                for j in range(base.shape[0]):
                    if base[j] != incomplete[j]:
                        warnings.warn(
                            f'{s}, {ses_id}, {type} - missing data {base[j]}')
                        incomplete = np.insert(
                            np.asarray(incomplete), j, np.nan)
                        this_data = np.insert(this_data, j, np.nan, axis=0)

            # if Data is None:
            #     Data = np.zeros((len(subj), C.shape[0], C.shape[1]))
            # Data[i, :, :] = this_data
            Data_list.append(this_data[np.newaxis, ...])

        # concatenate along the first dimension (subjects)
        Data = np.concatenate(Data_list, axis=0)
        # Ensure that infinite values (from div / 0) show up as NaNs
        Data[np.isinf(Data)] = np.nan
        return Data, info_com

    def get_info(self, ses_id='ses-s1',
                 type=None, subj=None, fields=None):
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

        # Deal with subset of subject option
        if subj is None:
            subj = np.arange(T.shape[0])

        if type is None:
            type = self.default_type

        max = 0

        # Loop over the different subjects to find the most complete info
        for s in T.participant_id.iloc:
            # Get an check the information
            info_raw = pd.read_csv(self.data_dir.format(s)
                                   + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t')
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

    def group_average_data(self, ses_id=None,
                           type=None,
                           atlas='SUIT3'):
        """Loads group data in SUIT space from a standard experiment structure
        averaged across all subjects. Saves the results as CIFTI files in the data/group directory.
        Args:
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        if ses_id is None:
            ses_id = self.sessions[0]
        if type is None:
            type = self.default_type

        data, info = self.get_data(space=atlas, ses_id=ses_id,
                                   type=type)
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
                    f'/group_{ses_id}_info-{type}.tsv', sep='\t', index=False)

    def plot_cerebellum(self, subject='group', sessions=None, atlas='SUIT3', type=None, cond='all', savefig=False, cmap='hot', colorbar=False):
        """Loads group data in SUIT3 space from a standard experiment structure
        averaged across all subjects and projects to SUIT flatmap. Saves the results as .png figures in the data/group/figures directory.
        Args:
            sub (str, optional): Subject string. Defaults to group to plot data averaged across all subjects.
            session (str, optional): Session string. Defaults to first session of in session list of dataset.
            atlas (str, optional): Atlas string. Defaults to 'SUIT3'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            cond (str or list): List of condition indices (e.g. [0,1,2] for the first three conditions) or 'all'. Defaults to 'all'.
            savefig (str, optional): Boolean indicating whether figure should be saved. Defaults to 'False'.
            cmap (str, optional): Matplotlib colour map. Defaults to 'hot'.
            colorbar (str, optional): Boolean indicating whether colourbar should be plotted in figure. Defaults to 'False'.


        """
        if sessions is None:
            sessions = self.sessions
        if type is None:
            type = self.default_type
        if subject == 'all':
            subjects = self.get_participants().participant_id.tolist()
        else:
            subjects = [subject]

        atlasmap, atlasinfo = am.get_atlas(atlas, self.atlas_dir)

        for sub in subjects:
            print(f'Plotting {sub}')
            for session in sessions:
                info = self.data_dir.split(
                    '/{0}')[0] + f'/{sub}/data/{sub}_{session}_info-{type}.tsv'
                data = self.data_dir.split(
                    '/{0}')[0] + f'/{sub}/data/{sub}_space-{atlas}_{session}_{type}.dscalar.nii'

                # Load average
                C = nb.load(data)
                D = pd.read_csv(info, sep='\t')
                X = C.get_fdata()
                # limes = [X[np.where(~np.isnan(X))].min(), X[np.where(~np.isnan(X))].max()] # cannot use nanmax or nanmin because memmap does not have this attribute
                limes = [np.percentile(X[np.where(~np.isnan(X))], 5), np.percentile(
                    X[np.where(~np.isnan(X))], 95)]

                if cond == 'all':
                    conditions = D[self.cond_name]
                else:
                    conditions = D[self.cond_name][cond]

                # -- Plot each condition in seperate figures --
                dest_dir = self.data_dir.split('/{0}')[0] + f'/{sub}/figures/'
                Path(dest_dir).mkdir(parents=True, exist_ok=True)
                for c in conditions:
                    condition_name = c.strip()
                    D[D[self.cond_name] == c].index
                    Nifti = atlasmap.data_to_nifti(
                        X[D[D[self.cond_name] == c].index, :].mean(axis=0))
                    surf_data = suit.flatmap.vol_to_surf(
                        Nifti, space=atlasinfo['normspace'])
                    fig = suit.flatmap.plot(
                        surf_data, render='matplotlib', new_figure=True, cscale=limes, cmap=cmap, colorbar=colorbar)
                    fig.set_title(condition_name)

                    # save figure
                    if savefig:
                        plt.savefig(
                            dest_dir + f'{sub}_{session}_{condition_name}.png')
                    plt.clf()

    
    def condense_data(self,data, info, type, participant_id=None, ses_id=None):
        """Empty function to be overwritten by child classes."""
        return data, info

class DataSetNative(DataSet):
    """Data set with estimates data stored as
    nifti-files in Native space.
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
        fnames = [
            f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

    def get_indiv_atlasmaps(self, atlas, sub, ses_id, smooth=None):
        atlas_maps = []
        if atlas.structure == 'cerebellum':
            deform = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            if atlas.name[0:4] != 'SUIT':
                deform1, m = am.get_deform(
                    self.atlas_dir, atlas.name, source='SUIT2')
                deform = [deform1, deform]
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.name == 'fs32k':
            for i, hem in enumerate(['L', 'R']):
                adir = self.anatomical_dir.format(sub)
                edir = self.estimates_dir.format(sub)
                pial = adir + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
                mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(atlas.vertex[i],
                                                  white, pial, mask))
                atlas_maps[i].build()

        elif atlas.name == 'fs32k_Asym':
            for i, hem in enumerate(['L', 'R']):
                adir = self.anatomical_dir.format(sub)
                edir = self.estimates_dir.format(sub)
                pial = adir + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
                mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(atlas.vertex_mask[i],
                                                  white, pial, mask))
                atlas_maps[i].build()
        return atlas_maps

    def extract_all(self,
                    ses_id='ses-s1',
                    type='CondHalf',
                    atlas='SUIT3',
                    smooth=2.0):
        """Extracts data in Volumetric space from a dataset in which the data is stored in Native space. Saves the results as CIFTI files in the data directory.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type for condense_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        myatlas, _ = am.get_atlas(atlas, self.atlas_dir)
        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:

            print(f'Atlasmap {s}')
            atlas_maps = self.get_indiv_atlasmaps(myatlas, s, ses_id,
                                                  smooth=smooth)
            print(f'Extract {s}')
            fnames, info = self.get_data_fnames(s, ses_id)

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
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)


class DataSetMNIVol(DataSet):
    """Data set with estimates data stored as
    nifti-files in MNI space.
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
        fnames = [
            f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

    def get_group_atlasmaps(self, atlas, sub=None, smooth=None):
        """ Gets group atlasmaps.
        Assumes that all scans are in the same space (self.group_space)

        Args:
            sub (str; optional): Subject
        Returns:
            atlas_maps (list): List of atlasmaps

        """
        atlas_maps = []
        if atlas.structure == 'cerebellum':
            deform = self.atlas_dir + \
                f'/{self.group_space}/{self.group_space}_space-SUIT_xfm.nii'
            if atlas.name[0:4] != 'SUIT':
                deform1, m = am.get_deform(self.atlas_dir, atlas.name, 'SUIT2')
                deform = [deform1, deform]
            mask = self.atlas_dir + \
                f'/{self.group_space}/{self.group_space}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world, deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.name == 'fs32k':
            for i, hem in enumerate(['L', 'R']):
                adir = self.anatomical_dir.format(sub)
                pial = adir + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
                mask = self.atlas_dir + \
                    f'/{self.group_space}/{self.group_space}_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(atlas.vertex[i],
                                                  white, pial, mask))
                atlas_maps[i].build()
        return atlas_maps

    def extract_all(self,
                    ses_id='ses-01',
                    type='CondHalf',
                    atlas='SUIT3',
                    smooth=None):
        """Extracts data in Volumetric space from a dataset in which the data is stored in Native space. Saves the results as CIFTI files in the data directory.
        Args:
            type (str, optional): Type for condense_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        myatlas, _ = am.get_atlas(atlas, self.atlas_dir)

        # extract and save data for each participant
        T = self.get_participants()
        for idx, s in enumerate(T.participant_id):
            if idx == 0 and myatlas.structure == 'cerebellum':
                print(f'Atlasmap group')
                atlas_maps = self.get_group_atlasmaps(myatlas,
                                                      smooth=smooth)
            elif myatlas.name == 'fs32k':
                atlas_maps = self.get_group_atlasmaps(myatlas, sub=s,
                                                      smooth=smooth)
            print(f'Extract {s}')
            fnames, info = self.get_data_fnames(s, ses_id)
            data = am.get_data_nifti(fnames, atlas_maps)
            data, info = self.condense_data(data, info, type)
            # Write out data as CIFTI file
            C = myatlas.data_to_cifti(data, info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)


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
        myatlas, _ = am.get_atlas(atlas, self.atlas_dir)
        # Get the correct map into CIFTI-format
        if isinstance(myatlas, am.AtlasVolumetric):
            deform, mask = am.get_deform(self.atlas_dir,
                                         target=myatlas,
                                         source='MNIAsymC2')
            atlas_map = am.AtlasMapDeform(myatlas.world,
                                          deform, mask)
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
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)


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
        info['half'] = 2 - (info.run < 9)
        n_cond = np.max(info.cond_num)
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
            B = np.ones((data_info.shape[0],))

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
        self.default_type = 'Net69Run'
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

    def average_within_Icos(self, label_file, data, atlas = "fs32k"):
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
        atlas, ainfo = am.get_atlas(atlas,self.atlas_dir)

        # create label_vector by passing on the label file
        # Set unite_struct to true if you want to integrate over left and right hemi 
        atlas.get_parcel(label_file, unite_struct = False)

        # use agg_parcel to aggregate data over parcels and get the list of unique parcels
        parcel_data, parcels =  agg_parcels(data, atlas.label_vector, fcn=np.nanmean)

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
            B = matrix.indicator(data_info.half, positive=True)

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
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
                                    ['half', 'cond_num'],
                                    ['run', 'reg_num'])
            data_info['names'] = [
                f'{d.cond_name.strip()}-half{d.half}' for i, d in data_info.iterrows()]

            # Baseline substraction
            B = matrix.indicator(data_info.half, positive=True)

        elif type == 'CondRun':
            data_info, C = agg_data(info,
                                    ['run', 'cond_num'],
                                    [])

            data_info['names'] = [
                f'{d.cond_name.strip()}-run{d.run:02d}' for i, d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.run, positive=True)
        elif type == 'CondAll':
            data_info, C = agg_data(info,
                                    ['cond_num'],
                                    ['run', 'half', 'reg_num'])
            # Baseline substraction
            B = np.ones((data_info.shape[0],))

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
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
        fnames = [
            f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

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

    def extract_all_suit(self, ses_id='ses-archi', type='CondHalf', atlas='SUIT3'):
        """Extracts data in SUIT space - we need to overload this from the standard,
        as the voxel-orientation (and therefore the atlasmap) is different from session to session in IBC.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        suit_atlas, _ = am.get_atlas(atlas, self.atlas_dir)
        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:
            print(f'Atlasmap {s}')
            deform = self.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
            if atlas[0:7] == 'MNISymC':
                xfm_name = self.atlas_dir + \
                    '/tpl-SUIT/tpl-SUIT_space-MNI152NLin2009cSymC_xfm.nii'
                deform = [xfm_name, deform]

            mask = self.estimates_dir.format(
                s) + f'/{ses_id}/{s}_{ses_id}_mask.nii'
            add_mask = self.suit_dir.format(s) + f'/{s}_desc-cereb_mask.nii'
            atlas_map = am.AtlasMapDeform(suit_atlas.world, deform, mask)
            atlas_map.build(smooth=2.0, additional_mask=add_mask)

            print(f'Extract {s}')
            fnames, info = self.get_data_fnames(s, ses_id)
            data = am.get_data_nifti(fnames, [atlas_map])
            data, info = self.condense_data(data, info, type,
                                            participant_id=s, ses_id=ses_id)
            C = suit_atlas.data_to_cifti(data[0], info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)


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
        self.group_space = 'tpl-MNI152NLin6Asym'
        self.sessions = ['ses-motor']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def condense_data(self, data, info,
                      type='CondHalf'):
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

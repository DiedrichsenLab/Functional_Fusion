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
import glob
import Functional_Fusion.util as util
import Functional_Fusion.matrix as matrix
import Functional_Fusion.atlas_map as am
import scipy.linalg as sl
import nibabel as nb
import nitools as nt
from numpy import eye,zeros,ones,empty,nansum, sqrt
from numpy.linalg import pinv,solve
import warnings

def get_dataset(base_dir,dataset,atlas='SUIT3',sess='all',type=None):
    """get_dataset

    Args:
        base_dir (str): Basis directory for the Functional Fusion repro
        dataset (str): Data set indicator
        atlas (str): Atlas indicator. Defaults to 'SUIT3'.
        sess (str or list): Sessiions. Defaults to 'all'.
        type (str): 'CondHalf','CondRun', etc....

    Returns:
        _type_: _description_
    """
    fiel = None
    if dataset.casefold() == 'MDTB'.casefold():
        my_dataset = DataSetMDTB(base_dir + '/MDTB')
        fiel = ['study','half','common','cond_name','cond_num','cond_num_uni','common']
        if type=='CondRun':
            fiel = fiel+['run']
        # Extract all sessions
    elif dataset.casefold() == 'Pontine'.casefold():
        my_dataset = DataSetPontine(base_dir + '/Pontine')
        fiel = ['task_name','task_num','half']
    elif dataset.casefold() == 'Nishimoto'.casefold():
        my_dataset = DataSetNishi(base_dir + '/Nishimoto')
        fiel = ['task_name','reg_id','half']
    elif dataset.casefold() == 'IBC'.casefold():
        my_dataset = DataSetIBC(base_dir + '/IBC')
    elif dataset.casefold() == 'HCP'.casefold():
        my_dataset = DataSetHcpResting(base_dir + '/HCP')
    elif dataset.casefold() == 'Demand'.casefold():
        my_dataset = DataSetDemand(base_dir + '/Demand')
    else:
        raise(NameError('Unknown data set'))

    # Get defaults sessions from dataset itself
    if sess=='all':
        sess=my_dataset.sessions
    if type is None:
        type = my_dataset.default_type

    # Load all data and concatenate
    # across sessions
    info_l = []
    data_l = []
    for s in sess:
        dat,inf = my_dataset.get_data(atlas,s,type,fields=fiel)
        data_l.append(dat)
        inf['sess']=[s]*inf.shape[0]
        info_l.append(inf)

    info = pd.concat(info_l,ignore_index=True,sort=False)
    data = np.concatenate(data_l,axis=1)
    return data,info,my_dataset

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
        resms = data[i][-1,:]
        data_n.append(data[i][0:-1,:])
        resms[resms<=0]=np.nan
        data_n[i] = data_n[i] / np.sqrt(resms)
    return data_n

def agg_data(info,by,over,subset=None):
    """ Aggregates data safely by sorting them by the fields in "by"
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
        info_n['n_rep']=np.ones((info_n.shape[0],))

    # Other contains the fields to remain constant
    other = list(info.columns.values)
    for ov in over+by:
        other.remove(ov)

    # Define operations on data
    operations = {'n_rep':np.sum}
    for o in other:
        operations[o]=max

    # Group the new data frame
    info_gb = info_n.groupby(by)
    data_info = info_gb.agg(operations).reset_index()

    # Build contrast matrix for averaging
    C = np.zeros((info.shape[0],data_info.shape[0]))
    for i,(k,v) in enumerate(info_gb.indices.items()):
        C[indx[v],i]=1
    return data_info,C

def optimal_contrast(data,C,X,reg_in=None,baseline=None):
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
    Cn = sl.block_diag(C,np.eye(num_nointerest))
    # Make new design matrix
    Xn = X @ Cn
    # Loop over the data:
    data_new = []
    for i in range(len(data)):
        # Append the regressors of no interest regressors
        dat = np.concatenate([data[i],
                    np.zeros((num_nointerest,data[i].shape[1]))])
        # Do the averaging / reweighting:
        d = solve(Xn.T @ Xn, Xn.T @ X @ dat)
        # Subset to the contrast of interest
        if reg_in is not None:
            d=d[reg_in,:]
        # Now subtract baseline
        if baseline is not None:
            Q = d.shape[0]
            R = eye(Q)-baseline @ pinv(baseline)
            d = R @ d
        # Put the data in a list:
        data_new.append(d)
    return data_new

def reliability_within_subj(X,part_vec,cond_vec,
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
        r = np.zeros((n_subj,n_part,X.shape[2]))
    else:
        r = np.zeros((n_subj,n_part))
    Z = matrix.indicator(cond_vec)
    for s in np.arange(n_subj):
        for pn,part in enumerate(partitions):
            i1 = part_vec==part
            X1= pinv(Z[i1,:]) @ X[s,i1,:]
            i2 = part_vec!=part
            X2 = pinv(Z[i2,:]) @ X[s,i2,:]
            if subtract_mean:
                X1 -= X1.mean(axis=0)
                X2 -= X2.mean(axis=0)
            if voxel_wise:
                r[s,pn,:] = nansum(X1*X2,axis=0)/ \
                    sqrt(nansum(X1*X1,axis=0)
                        *nansum(X2*X2,axis=0))
            else:
                r[s,pn] = nansum(X1*X2)/sqrt(nansum(X1*X1)*nansum(X2*X2))
    return r

def reliability_between_subj(X,cond_vec=None,
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
        r = np.zeroes((n_subj,X.shape[2]))
    else:
        r = np.zeros((n_subj,))
    for s,i in enumerate(subj_vec):
        X1= pinv(Z) @ X[s,:,:]
        i2 = subj_vec!=s
        X2 = pinv(Z) @ X[i2,:,:].mean(axis=0)
        if subtract_mean:
            X1 -= X1.mean(axis=0)
            X2 -= X2.mean(axis=0)
        if voxel_wise:
            r[i,:] = nansum(X1*X2,axis=0)/ \
                    sqrt(nansum(X1*X1,axis=0)
                        *nansum(X2*X2,axis=0))
        else:
            r[i] = nansum(X1*X2)/sqrt(nansum(X1*X1)*nansum(X2*X2))
    return r

def reliability_maps(base_dir,dataset_name,
                    atlas = 'MNISymC3',
                    subtract_mean=True):
    """    Calculates the average within subject reliability maps across sessions for a single data

    Args:
        base_dir (str / path): Base directory
        dataset_name (str): Name of data set
        atlas (str): _description_. Defaults to 'MNISymC3'.
        subtract_mean (bool): Remove the mean per voxel before correlation calc?

    Returns:
        _type_: _description_
    """
    data,info,dataset = get_dataset(base_dir,dataset_name,atlas=atlas)
    n_sess = len(dataset.sessions)
    n_vox = data.shape[2]
    Rel = np.zeros((n_sess,n_vox))
    for i,s in enumerate(dataset.sessions):
        indx = info.sess==s
        r = reliability_within_subj(data[:,indx,:],
                    part_vec=info[dataset.part_ind][indx],
                    cond_vec=info[dataset.cond_ind][indx],
                    voxel_wise=True,
                    subtract_mean=subtract_mean)
        Rel[i,:] = np.nanmean(np.nanmean(r,axis=0),axis=0)
    return Rel,dataset.sessions

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
        self.base_dir  = base_dir
        self.surface_dir = base_dir + '/derivatives/{0}/anat'
        self.anatomical_dir = base_dir + '/derivatives/{0}/anat'
        self.estimates_dir = base_dir + '/derivatives/{0}/estimates'
        self.func_dir = base_dir + '/derivatives/{0}/func'
        self.suit_dir = base_dir + '/derivatives/{0}/suit'
        self.data_dir = base_dir + '/derivatives/{0}/data'
        # assume that the common atlas directory is on the level before
        self.atlas_dir = os.path.join(os.path.dirname(base_dir),'Atlases')
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
        self.part_info = pd.read_csv(self.base_dir + '/participants.tsv',delimiter='\t')
        return self.part_info

    def get_data(self, space='SUIT3', ses_id='ses-s1',
                 type='CondHalf', subj=None, fields=None):
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
        # Assemble the data
        Data = None
        # Deal with subset of subject option
        if subj is None:
            subj = np.arange(T.shape[0])

        max = 0

        # Loop over the different subjects to find the most complete info
        for i, s in enumerate(T.participant_id.iloc):
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
        for i, s in enumerate(T.participant_id.iloc):
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
                        warnings.warn(f'{s}, {ses_id}, {type} - missing data {base[j]}')
                        incomplete = np.insert(np.asarray(incomplete), j, np.nan)
                        this_data = np.insert(this_data, j, np.nan, axis=0)

            if Data is None:
                Data = np.zeros((len(subj), C.shape[0], C.shape[1]))
            Data[i, :, :] = this_data

        # Ensure that infinite values (from div / 0) show up as NaNs
        Data[np.isinf(Data)] = np.nan
        return Data, info_com

    def group_average_data(self, ses_id='ses-s1',
                                 type='CondHalf',
                                 atlas='SUIT3'):
        """Loads group data in SUIT space from a standard experiment structure
        averaged across all subjects. Saves the results as CIFTI files in the data/group directory.
        Args:
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
            info_column (str, optional): Column of info tsv file for which each average should be calculated. Defaults to 'task_name'
        """

        data, info = self.get_data(space=atlas, ses_id=ses_id,
                                   type=type)
        # average across participants
        X = np.nanmean(data,axis=0)
        # make output cifti
        s = self.get_participants().participant_id[0]
        C = nb.load(self.data_dir.format(s) +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
        C = nb.Cifti2Image(dataobj=X, header=C.header)
        # save output
        dest_dir = op.join(self.data_dir.format('group'))
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        nb.save(C, dest_dir +
                f'/group_{ses_id}_space-{atlas}_{type}.dscalar.nii')
        info.drop(columns=['sn']).to_csv(dest_dir +
            f'/group_{ses_id}_info-{type}.tsv', sep='\t',index=False)

class DataSetNative(DataSet):
    """Data set with estimates data stored as
    nifti-files in Native space.
    """
    def get_data_fnames(self,participant_id,session_id=None):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T=pd.read_csv(dirw+f'/{participant_id}_{session_id}_reginfo.tsv',sep='\t')
        fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i,t in T.iterrows()]
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

    def get_indiv_atlasmaps(self,atlas,sub,ses_id,smooth=None):
        atlas_maps=[]
        if atlas.structure=='cerebellum':
            deform = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            if atlas.name[0:4]!='SUIT':
                deform1,m = am.get_deform(self.atlas_dir,atlas.name,'SUIT2') 
                deform = [deform1, deform]
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'
            atlas_maps.append(am.AtlasMapDeform(atlas.world,deform, mask))
            atlas_maps[0].build(smooth=smooth)
        elif atlas.name=='fs32k':
            for i,hem in enumerate(['L','R']):
                adir = self.anatomical_dir.format(sub)
                edir = self.estimates_dir.format(sub)
                pial = adir + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
                mask = edir + f'/{ses_id}/{sub}_{ses_id}_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(atlas.vertex[i],
                            white,pial, mask))
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
        myatlas = am.get_atlas(atlas,self.atlas_dir)
        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:

            print(f'Atlasmap {s}')
            atlas_maps = self.get_indiv_atlasmaps(myatlas,s,ses_id,
                            smooth=smooth)
            print(f'Extract {s}')
            fnames,info = self.get_data_fnames(s,ses_id)
            data = am.get_data_nifti(fnames,atlas_maps)
            data,info = self.condense_data(data,info,type,
                    participant_id=s,ses_id=ses_id)
            # Write out data as CIFTI file
            C=myatlas.data_to_cifti(data,info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir + f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv',sep='\t', index = False)

class DataSetCifti(DataSet):
    """Data set that comes in HCP-format in already pre-extracted cifti files.
    """
    def get_data_fnames(self,participant_id,session_id=None):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T=pd.read_csv(dirw+f'/{participant_id}_{session_id}_reginfo.tsv',sep='\t')
        fnames = [f'{dirw}/{participant_id}_{session_id}_beta.dscalar.nii']
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.dscalar.nii')
        return fnames, T

    def extract_all(self,ses_id='ses-s1',type='CondHalf',atlas='SUIT3'):
        """Extracts cerebellar data. Saves the results as CIFTI files in the data directory.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        myatlas = am.get_atlas(atlas,self.atlas_dir)
        # Get the correct map into CIFTI-format 
        if isinstance(myatlas,am.AtlasVolumetric):
            deform,mask = am.get_deform(self.atlas_dir,
                target=myatlas,
                source='MNIAsymC2')
            atlas_map = am.AtlasMapDeform(myatlas.world,
                deform,mask)
            atlas_map.build(smooth=2.0)
        elif isinstance(myatlas,am.AtlasSurface):
            atlas_map = myatlas
        # Extract the data for each participant 
        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}')
            fnames,info = self.get_data_fnames(s,ses_id)
            data = am.get_data_cifti(fnames,[atlas_map])

            data,info = self.condense_data(data,info,type,
                    participant_id=s,ses_id=ses_id)
            C=myatlas.data_to_cifti(data[0],info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir + f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv',sep='\t', index = False)


class DataSetMDTB(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions=['ses-s1','ses-s2']
        self.default_type = 'CondHalf'
        self.cond_ind = 'cond_num_uni'
        self.part_ind = 'half'
        self.cond_name = 'cond_name'

    def condense_data(self,data,info,
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
        info['half']=2-(info.run<9)
        n_cond = np.max(info.cond_num)
        if type == 'CondHalf':
            data_info,C = agg_data(info,
                        ['half','cond_num'],
                        ['run','reg_num'],
                        subset=(info.instruction==0))
            data_info['names']=[f'{d.cond_name.strip()}-half{d.half}' for i,d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'CondRun':
            data_info,C = agg_data(info,
                        ['run','cond_num'],
                        [],
                        subset=(info.instruction==0))
            data_info['names']=[f'{d.cond_name}-run{d.run:02d}' for i,d in data_info.iterrows()]

            # Baseline substraction
            B = matrix.indicator(data_info.run,positive=True)
        elif type == 'CondAll':

            data_info,C = agg_data(info,
                        ['cond_num'],
                        ['run','half','reg_num'],
                        subset=(info.instruction==0))
            data_info['names']=[f'{d.cond_name}' for i,d in data_info.iterrows()]

            # Baseline substraction
            B = np.ones((data_info.shape[0],))

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
        X = np.load(dir+f'/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1],dtype=int)
        #  contrast for all instructions
        CI = matrix.indicator(info.run*info.instruction,positive=True)
        C = np.c_[C,CI]
        data_new = optimal_contrast(data_n,C,X,reg_in,baseline=B)

        return data_new, data_info

class DataSetHcpResting(DataSetCifti):
    def __init__(self, dir):
        super(DataSetHcpResting, self).__init__(base_dir=dir)
        # self.func_dir = self.base_dir + '/{0}/estimates'
        self.derivative_dir = self.base_dir + '/derivatives'
        self.sessions=['ses-s1','ses-s2']
        self.hem_name = ['cortex_left', 'cortex_right']
        self.default_type = 'NetRun'
        self.cond_ind = 'reg_id'
        self.cond_name = 'region_name'
        self.part_ind = 'half'

    def get_data_fnames(self, participant_id):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
        Returns:
            fnames (list): List of fnames
        """
        dirw = self.derivative_dir + f"/{participant_id}" + "/func"
        fnames = []
        for r in range(4):
            fnames.append(f'{dirw}/sub-{participant_id}_run-{r}_space-MSMSulc.dtseries.nii')
        return fnames

    def extract_all_suit(self, ses_id='ses-s1', type='IcoAll', atlas='SUIT3', res=162):
        """ MDTB extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'IcoAll': Single estimate per session, Correlation with cortical icosahedron parcels
                'IcoRun': Seperate estimates per run, Correlation with cortical icosahedron parcels
                'NetAll': Single estimate per session, Correlation with cortico-cerebellar resting-state networks
                'NetRun': Seperate estimates per run, Correlation with cortico-cerebellar resting-state networks
                    Defaults to 'IcoAll'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        suit_atlas = am.get_atlas(atlas, self.atlas_dir)

        # Get the deformation map from MNI to SUIT
        mni_atlas = self.atlas_dir + '/tpl-MNI152NLin6AsymC'
        deform = mni_atlas + '/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
        mask = mni_atlas + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
        atlas_map = am.AtlasMapDeform(suit_atlas.world, deform, mask)
        atlas_map.build(smooth=2.0)

        if type[0:3] == 'Ico':
            networks=None
            # Get the parcelation
            surf_parcel = []
            for i, h in enumerate(['L', 'R']):
                dir = self.atlas_dir + '/tpl-fs32k'
                gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
                surf_parcel.append(am.AtlasSurfaceParcel(self.hem_name[i], gifti))
            bpa = surf_parcel[0].get_parcel_axis() + surf_parcel[1].get_parcel_axis()
            seed_names=list(bpa.name)
        elif type[0:3] == 'Net':
            surf_parcel = None
            # Get the networks
            networkdir = self.base_dir + '/group_ica/dim_25/'
            networkimg = nb.load(networkdir +
                'melodic_IC.nii.gz')
            networks = networkimg.get_fdata()
            net_selected = pd.read_csv(
                networkdir + 'classified_components.txt', sep=', ', skiprows=[0], skipfooter=1, engine='python', header=None, names=['Network', 'Classification', 'IsNoise'], dtype="category")
            networks = networks[:, :, :,
                                net_selected.Classification == 'Signal']
            seed_names = [
                f'Network-{n+1:02}' for n in np.arange(networks.shape[-1])]

        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}')
            if ses_id == 'ses-s1':
                runs = [0, 1]
            elif ses_id == 'ses-s2':
                runs = [2, 3]
            else:
                raise ValueError('Unknown session id.')

            coef = self.get_cereb_connectivity(
                s, atlas_map, runs=runs, type=type, cortical_atlas_parcels=surf_parcel, networks=networks)

            if type[3:7] == 'All':  # Average across runs
                coef = np.nanmean(coef, axis=0)

                # Make info structure
                reg_ids = np.arange(len(seed_names)) + 1
                info = pd.DataFrame({'sn': [s] * coef.shape[0],
                                    'sess': [ses_id] * coef.shape[0],
                                     'half': [1] * coef.shape[0],
                                     'reg_id': reg_ids,
                                     'region_name': seed_names,
                                     'names': seed_names})

            elif type[3:7] == 'Run':  # Concatenate over runs
                coef = np.concatenate(coef, axis=0)

                # Make info structure
                run_ids = np.repeat(runs, int(coef.shape[0] / len(runs)))
                reg_ids = np.tile(np.arange(len(seed_names)), 2) + 1
                names = ["{}_run-{}".format(reg_name, run_id)
                         for reg_name, run_id in zip(list(seed_names) * 2, run_ids)]
                info = pd.DataFrame({'sn': [s] * coef.shape[0],
                                    'sess': [ses_id] * coef.shape[0],
                                    'run': run_ids,
                                     'half': 2 - (run_ids < run_ids[-1]),
                                    'reg_id': reg_ids,
                                     'region_name': list(seed_names) * 2,
                                    'names': names})

                # update brain parcel axis (repeat names)
                # bpa = bpa + bpa


            # --- Save cerebellar data as dscalar CIFTI-file and write info to tsv ---
            C = am.data_to_cifti(coef, [atlas_map], info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)

            # --- Build a connectivity CIFTI-file and save ---
            # bmc = suit_atlas.get_brain_model_axis()
            # header = nb.Cifti2Header.from_axes((bpa, bmc))
            # cifti_img = nb.Cifti2Image(dataobj=coef, header=header)
            # nb.save(cifti_img, dest_dir +
            #         f'/{s}_space-{atlas}_{ses_id}_{type}_{res}.dpconn.nii')

    def extract_all_fs32k(self, ses_id='ses-s1', type='IcoAll', res=162):
        """ MDTB extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'IcoAll': Single estimate per session,
                          Correlation with cortical icosahedron parcels
                'IcoRun': Seperate estimates per run,
                          Correlation with cortical icosahedron parcels
                'NetAll': Single estimate per session,
                          Correlation with cortico-cerebellar resting-state networks
                'NetRun': Seperate estimates per run,
                          Correlation with cortico-cerebellar resting-state networks
                Defaults to 'IcoAll'.
            res: the resolution of underlying icosahedron. Default 162
        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        # Make the atlas object
        mask_L = self.atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii'
        mask_R = self.atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii'
        atlas = am.AtlasSurface('fs32k', mask_gii=[mask_L, mask_R],
                                   structure=['CORTEX_LEFT', 'CORTEX_RIGHT'])
        cortex_mask = [mask_L, mask_R]
        bmc = atlas.get_brain_model_axis()
        seed_names=[]

        if type[0:3] == 'Ico':
            networks=None
            # Get the parcelation
            surf_parcel = []
            for i, h in enumerate(['L', 'R']):
                dir = self.atlas_dir + '/tpl-fs32k'
                gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
                mask = dir + f'/tpl-fs32k_hemi-{h}_mask.label.gii'
                surf_parcel.append(am.AtlasSurfaceParcel(self.hem_name[i], gifti, mask_gii=mask))
            bpa = surf_parcel[0].get_parcel_axis() + surf_parcel[1].get_parcel_axis()
            seed_names=list(bpa.name)
        elif type[0:3] == 'Net':
            surf_parcel = None
            # Get the networks
            networkdir = self.base_dir + '/group_ica/dim_25/'
            networkimg = nb.load(networkdir + 'melodic_IC.nii.gz')
            networks = networkimg.get_fdata()
            net_selected = pd.read_csv(networkdir + 'classified_components.txt',
                                       sep=', ', skiprows=[0], skipfooter=1,
                                       engine='python', header=None,
                                       names=['Network', 'Classification', 'IsNoise'],
                                       dtype="category")
            networks = networks[:, :, :, net_selected.Classification == 'Signal']
            seed_names = [f'Network-{n+1:02}' for n in np.arange(networks.shape[-1])]
        else:
            raise NameError("type must start with either 'Ico' or 'Net'!")
        # Making cifti2 axis for these network name
        bpa = nb.cifti2.ScalarAxis(seed_names)

        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}, type {type}')
            if ses_id == 'ses-s1':
                runs = [0, 1]
            elif ses_id == 'ses-s2':
                runs = [2, 3]
            else:
                raise ValueError('Unknown session id.')

            coef = self.get_cortical_connectivity(s, cortex_mask=cortex_mask, runs=runs, type=type,
                                                  cortical_atlas_parcels=surf_parcel,
                                                  networks=networks)

            if type[3:7] == 'All':  # Average across runs
                coef = [np.nanmean(c, axis=0) for c in coef]

                # Make info structure
                info = []
                for i, d in enumerate(coef):
                    reg_ids = np.arange(len(seed_names)) + 1
                    this_info = pd.DataFrame({'sn': [s] * d.shape[0],
                                              'hemis': i + 1,
                                              'sess': [ses_id] * d.shape[0],
                                              'half': [1] * d.shape[0],
                                              'reg_id': reg_ids,
                                              'region_name': seed_names,
                                              'names': seed_names})
                    info.append(this_info)

            elif type[3:7] == 'Run':  # Concatenate over runs
                coef = [np.concatenate(c, axis=0) for c in coef]

                # Make info structure
                info = []
                for i, d in enumerate(coef):
                    run_ids = np.repeat(runs, int(d.shape[0] / len(runs)))
                    reg_ids = np.tile(np.arange(len(seed_names)), 2) + 1
                    names = ["{}_run-{}".format(reg_name, run_id)
                             for reg_name, run_id in zip(list(seed_names) * 2, run_ids)]
                    this_info = pd.DataFrame({'sn': [s] * d.shape[0],
                                              'hemis': i + 1,
                                              'sess': [ses_id] * d.shape[0],
                                              'run': run_ids,
                                              'half': 2 - (run_ids < run_ids[-1]),
                                              'reg_id': reg_ids,
                                              'region_name': list(seed_names) * 2,
                                              'names': names})
                    info.append(this_info)
                # update brain parcel axis (repeat names)
                bpa = bpa + bpa

            info = pd.concat([info[0], info[1]])

            # --- Build a connectivity CIFTI-file and save ---
            print(f'Writing {s}, type {type} ...')
            header = nb.Cifti2Header.from_axes((bpa, bmc))
            cifti_img = nb.Cifti2Image(dataobj=np.c_[coef[0], coef[1]], header=header)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(cifti_img, dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_'
                                          f'{res}.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_space-fs32k_{ses_id}_info-{type}_{res}.tsv',
                           sep='\t', index=False)

    def extract_ts_volume(self,
                participant_id,
                atlas_map,
                ses_id=[0,1,2,3],
                type='Run'):
        """ Returns the time series data for an atlas map
                runs=[0,1,2,3]):

        Args:
            participant_id (_type_): _description_
            atlas_map (_type_): _description_
        """
        # get the file name for the cifti time series
        fnames, info = self.get_data_fnames(participant_id)
        ts_volume = []
        for r in runs:
            # load the cifti
            ts_cifti = nb.load(fnames[r])

            # get the ts in volume for subcorticals
            ts_vol = nt.volume_from_cifti(ts_cifti)
            # transfer to suit space applying the deformation
            ts_vol = am.get_data4D(ts_vol, [atlas_map])
            ts_volume.append(ts_vol[0])
        return ts_volume

    def extract_ts_surface(self,
                    participant_id,
                    atlas_parcels,
                    runs=[0,1,2,3]):
        """Returns the information from the CIFTI file
        in the 32K surface for left and right hemisphere.

        Args:
            participant_id (_type_): _description_
            atlas_parcel (_type_): _description_
        """
        hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT', 'CIFTI_STRUCTURE_CORTEX_RIGHT']
        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef = None
        ts_cortex=[]
        for r in runs:
            # load the cifti
            ts_cifti = nb.load(fnames[r])

            # get the ts in surface for corticals
            ts_32k = util.surf_from_cifti(ts_cifti,hem_name)
            ts_parcel = []
            for hem in range(2):

                # get the average within parcels
                ts_parcel.append(
                    atlas_parcels[hem].agg_data(ts_32k[hem])
                    )

            # concatenate them into a single array for correlation calculation
            ts_cortex.append(ts_parcel)
        return ts_cortex  # shape (n_tessl,P)

    def get_network_timecourse(self, networks, ts):
        """Regresses the group spatial map into the fMRI run.
        Returns the run-specific network timecourse.

        Args:
            networks (np.arry): 4D Network data of the signal components
                (default input networks are in MNI Space: 91 x 109 x 91 x nComponents )
            ts_vol (<nibabel CIFTI image object>): fMRI timeseries in volume
                Has to be in the same space as networks (91 x 109 x 91 x nTimepoints )
        Returns:
            ts_networks (np.ndarray):
                A numpy array (nTimepoints x nNetworks) with the fMRI timecourse for
                each resting-state network
        """
        X = networks.reshape(-1,networks.shape[3])
        Y = ts.reshape(-1, ts.shape[3])
        ts_networks = np.matmul(np.linalg.pinv(X), Y)

        return ts_networks  # shape (n_tessl, time_course)

    def get_cereb_connectivity(self, participant_id, cereb_atlas_map,runs=[0, 1, 2, 3],
                               type='IcoAll', cortical_atlas_parcels=None, networks=None):
        """Uses the original CIFTI files to produce cerebellar connectivity
           file

        """
        if type[0:3] == 'Ico':
            hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT', 'CIFTI_STRUCTURE_CORTEX_RIGHT']

        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef = None
        for r,run in enumerate(runs):
            # load the cifti
            ts_cifti = nb.load(fnames[run])

            # get the ts in volume for subcorticals
            ts_vol = util.volume_from_cifti(ts_cifti)
            # transfer to suit space applying the deformation
            ts_cerebellum = am.get_data4D(ts_vol, [cereb_atlas_map])
            ts_cerebellum = ts_cerebellum[0]

            if type[0:3] == 'Ico':
                # get the ts in surface for corticals
                ts_32k = util.surf_from_cifti(ts_cifti,hem_name)
                ts_parcel = []
                for hem in range(2):

                    # get the average within parcels
                    ts_parcel.append(
                        cortical_atlas_parcels[hem].agg_data(ts_32k[hem])
                        )

                # concatenate them into a single array for correlation calculation
                ts_seed = np.concatenate(ts_parcel, axis = 1)
            elif type[0:3] == 'Net':
                # Regress network spatial map into the run's wholebrain data
                # (returns run-specific timecourse for each network)
                ts = ts_vol.get_fdata()
                ts_seed = self.get_network_timecourse(networks, ts)
                ts_seed = ts_seed.T


            # Standardize the time series for easier calculation
            ts_cerebellum = util.zstandarize_ts(ts_cerebellum)
            ts_seed = util.zstandarize_ts(ts_seed)

            # Correlation calculation
            if coef is None:
                coef=np.empty((len(runs),ts_seed.shape[1],
                                ts_cerebellum.shape[1]))
            N = ts_cerebellum.shape[0]
            coef[r,:,:] = ts_seed.T @ ts_cerebellum / N

        return coef  # shape (n_tessl,P)

    def get_cortical_connectivity(self, participant_id, cortex_mask=None,
                                  runs=[0,1,2,3], type='IcoAll',
                                  cortical_atlas_parcels=None,
                                  networks=None):
        """Uses the original CIFTI files to produce cortical connectivity
           file
        Args:
            participant_id (int): participant id
            cortex_mask (list): the list of L and R hemis cortex mask
            runs: index of runs
            type: the underlying cortical parcellation and type of extraction
                'IcoAll': Single estimate per session,
                          Correlation with cortical icosahedron parcels
                'IcoRun': Seperate estimates per run,
                          Correlation with cortical icosahedron parcels
                'NetAll': Single estimate per session,
                          Correlation with cortico-cerebellar resting-state networks
                'NetRun': Seperate estimates per run,
                          Correlation with cortico-cerebellar resting-state networks
                Defaults to 'IcoAll'.
            cortical_atlas_parcels: cortical random tessellation parcel object
            networks: group ICA networks
        Returns:
            [coef_1,coef_2]: List of cortical functional connectivity.
                             [left hemisphere, right hemisphere]
        """
        hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT', 'CIFTI_STRUCTURE_CORTEX_RIGHT']

        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef_1, coef_2 = None, None
        for r,run in enumerate(runs):
            # load the cifti
            ts_cifti = nb.load(fnames[run])

            # get the ts in surface for corticals
            ts_32k = util.surf_from_cifti(ts_cifti, hem_name, mask_gii=cortex_mask)

            if type[0:3] == 'Ico':
                assert cortical_atlas_parcels is not None, \
                    "cortical_atlas_parcels must be given if extraction type is `Ico`!"
                ts_parcel = []
                for hem in range(2):
                    # get the average within parcels
                    ts_parcel.append(cortical_atlas_parcels[hem].agg_data(ts_32k[hem]))

                # concatenate them into a single array for correlation calculation
                ts_parcel = np.concatenate(ts_parcel, axis=1)
            elif type[0:3] == 'Net':
                assert networks is not None, \
                    "networks must be given if extraction type is `Net`!"
                # Regress network spatial map into the run's wholebrain data
                # (returns run-specific timecourse for each network)
                ts_vol = util.volume_from_cifti(ts_cifti)
                ts = ts_vol.get_fdata()
                ts_parcel = self.get_network_timecourse(networks, ts)
                ts_parcel = ts_parcel.T

            # Standardize the time series for easier calculation
            ts_cortex = [util.zstandarize_ts(timeseries) for timeseries in ts_32k]
            ts_parcel = util.zstandarize_ts(ts_parcel)

            # Correlation calculation
            if coef_1 is None:
                coef_1 = np.empty((len(runs),ts_parcel.shape[1],
                                   ts_cortex[0].shape[1])) # (runs, parcels, vertices)
            if coef_2 is None:
                coef_2 = np.empty((len(runs),ts_parcel.shape[1],
                                   ts_cortex[1].shape[1])) # (runs, parcels, vertices)

            N1 = ts_cortex[0].shape[0]
            N2 = ts_cortex[1].shape[0]
            coef_1[r,:,:] = ts_parcel.T @ ts_cortex[0] / N1
            coef_2[r,:,:] = ts_parcel.T @ ts_cortex[1] / N2

        return [coef_1,coef_2]  # shape (n_tessl,P)

class DataSetLanguage(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)

    def extract_all_suit(self, ses_id='ses-01', type='TaskAll', atlas='SUIT3'):
        """Extracts data in SUIT space from a standard experiment structure
        across all subjects. Saves the results as CIFTI files in the data directory.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        # Make the atlas object
        suit_atlas = am.get_atlas(atlas,self.atlas_dir)


        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}')
            data, info = self.extract_data(s, [atlas_map],
                                           ses_id=ses_id,
                                           type=type)
            C = am.data_to_cifti(data, [atlas_map], info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)

    def condense_data(self, participant_id,
                     atlas_maps,
                     ses_id,
                     type='TaskAll'):
        """ Language extraction of atlasmap locations
        from nii files

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'TaskAll': Contrasts with one estimate for all the whole experiment

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
        fnames, info = self.get_data_fnames(participant_id, ses_id)
        data = am.get_data3D(fnames, atlas_maps)
        # For debugging: data = [np.random.normal(0,1,(len(fnames),atlas_maps[0].P))]

        # Depending on the type, make a new contrast
        info['half'] = 2 - (info.run < 9)
        n_cond = np.max(info.reg_id)
        if type == 'TaskAll':

            # Make new data frame for the information of the new regressors
            ii = (info.run == 1) & (info.reg_id > 0)
            data_info = info[ii].copy().reset_index(drop=True)
            data_info['names'] = [
                f'{d.task_name.strip()}' for i, d in data_info.iterrows()]

           # Contrast for the regressors of interest
            reg = info.cond_num.copy()
            reg[info.instruction == 1] = 0
            C = matrix.indicator(reg, positive=True)  # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.instruction, positive=True)
            C = np.c_[C, CI]
            reg_in = np.arange(n_cond, dtype=int)

            # Baseline substraction
            B = matrix.indicator(data_info.run, positive=True)

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
        data_new = optimal_contrast(data_n, C, X, reg_in)

        return data_new, data_info

class DataSetPontine(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions=['ses-01']
        self.default_type = 'TaskHalf'
        self.cond_ind = 'task_num'
        self.cond_name = 'task_name'
        self.part_ind = 'half'

    def condense_data(self, data,info,
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
            data_info,C = agg_data(info,
                        ['half','reg_id'],
                        ['run','reg_num'],
                        subset=(info.reg_id>0))
            data_info['names']=[f'{d.task_name.strip()}-half{d.half}' for i,d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'TaskRun':

            data_info,C = agg_data(info,
                        ['run','reg_id'],
                        ['reg_num'],
                        subset=(info.reg_id>0))
            data_info['names']=[f'{d.task_name.strip()}-run{d.run}' for i,d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'TaskAll':
            data_info,C = agg_data(info,
                        ['reg_id'],
                        ['run','half','reg_num'],
                        subset=(info.reg_id>0))
            data_info['names']=[f'{d.task_name.strip()}' for i,d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.half,positive=True)

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1],dtype=int)
        CI = matrix.indicator(info.run*info.instruction,positive=True)
        C = np.c_[C,CI]

        data_new = optimal_contrast(data_n, C, X, reg_in)

        return data_new, data_info

class DataSetNishi(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions=['ses-01','ses-02']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'task_name'
        self.part_ind = 'half'

    def condense_data(self, data,info,
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
        info['half']=2-(info.run< (len(np.unique(info.run))/2+1))
        n_cond = np.max(info.reg_id)

        if type == 'CondHalf':
            data_info,C = agg_data(info,
                        ['half','cond_num'],
                        ['run','reg_num'])
            data_info['names']=[f'{d.cond_name.strip()}-half{d.half}' for i,d in data_info.iterrows()]


            # Baseline substraction
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'CondRun':
            data_info,C = agg_data(info,
                        ['run','cond_num'],
                        [])

            data_info['names']=[f'{d.cond_name.strip()}-run{d.run:02d}' for i,d in data_info.iterrows()]
            # Baseline substraction
            B = matrix.indicator(data_info.run,positive=True)
        elif type == 'CondAll':
            data_info,C = agg_data(info,
                        ['cond_num'],
                        ['run','half','reg_num'])
            # Baseline substraction
            B = np.ones((data_info.shape[0],))

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir+f'/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1],dtype=int)
        data_new = optimal_contrast(data_n,C,X,reg_in,baseline=B)

        return data_new, data_info

class DataSetIBC(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-archi',
                         'ses-clips4',
                         'ses-enumeration',
                         'ses-hcp1','ses-hcp2',
                         'ses-lyon1','ses-lyon2',
                         'ses-mathlang',
                         'ses-mtt1','ses-mtt2',
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
        self.part_info = pd.read_csv(self.base_dir + '/participants.tsv',delimiter='\t')
        return self.part_info[self.part_info.complete==1]

    def get_data_fnames(self,participant_id,session_id=None):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T=pd.read_csv(dirw+f'/{participant_id}_{session_id}_reginfo.tsv',sep='\t')
        fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i,t in T.iterrows()]
        fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

    def condense_data(self,data,info,
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
        info['n_rep']=np.ones((info.shape[0],))
        if type == 'CondHalf':
            data_info,C = agg_data(info,
                        ['half','cond_num'],
                        ['run','reg_num'])
            data_info['names']=[f'{d.cond_name.strip()}-half{d.half}' for i,d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)
        for i in range(len(data_n)):
            data_n[i]=pinv(C) @ data_n[i]
        return data_n, data_info

    def extract_all_suit(self,ses_id='ses-archi',type='CondHalf',atlas='SUIT3'):
        """Extracts data in SUIT space - we need to overload this from the standard,
        as the voxel-orientation (and therefore the atlasmap) is different from session to session in IBC.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        suit_atlas = am.get_atlas(atlas,self.atlas_dir)
        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:
            print(f'Atlasmap {s}')
            deform = self.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
            if atlas[0:7]=='MNISymC':
                xfm_name = self.atlas_dir + \
                           '/tpl-MNI152NLIn2009cSymC/tpl-SUIT_space-MNI152NLin2009cSymC_xfm.nii'
                deform = [xfm_name,deform]
            mask = self.estimates_dir.format(s) + f'/{ses_id}/{s}_{ses_id}_mask.nii'
            add_mask = self.suit_dir.format(s) + f'/{s}_desc-cereb_mask.nii'
            atlas_map = am.AtlasMapDeform(suit_atlas.world, deform, mask)
            atlas_map.build(smooth=2.0,additional_mask = add_mask)
            print(f'Extract {s}')
            data,info = self.extract_data(s,[atlas_map],
                                          ses_id=ses_id,
                                          type=type)
            C=suit_atlas.data_to_cifti(data[0],info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir + f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv',sep='\t', index = False)


class DataSetDemand(DataSetCifti):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions = ['ses-01']
        self.default_type = 'CondHalf'
        self.cond_ind = 'cond_num'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def condense_data(self,data,info,
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
        info['half']=(info.run%2)+1
        n_cond = np.max(info.reg_id)
        if type == 'CondHalf':
            data_info,C=agg_data(info,['half','reg_id'],['run'])
            data_info['names']=[f'{d.cond_name.strip()}-half{d.half}'
                    for i,d in data_info.iterrows()]
        elif type == 'CondAll':
            data_info,C=agg_data(info,['reg_id'],['half','run'])
            data_info['names']=[f'{d.cond_name.strip()}-half{d.half}'
                    for i,d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Combine with contrast
        for i in range(len(data_n)):
            data_n[i]=pinv(C) @ data_n[i]
        return data_n, data_info

class DataSetWMFS(DataSetNative):
    def __init__(self, dir):
        super().__init__(dir)
        self.sessions=['ses-01','ses-02']
        self.default_type = 'CondHalf'
        self.cond_ind = 'reg_id'
        self.cond_name = 'cond_name'
        self.part_ind = 'half'

    def condense_data(self,data,info,
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
        info['half']=2-(info.run<3)
        n_cond = np.max(info.loc[info.error == 0].reg_id)

        if type == 'CondHalf':
            # Make new data frame for the information of the new regressors

            data_info,C = agg_data(info,
                        ['half','cond_num'],
                        ['run','reg_num'],
                        subset=(info.error==0))
            data_info['names']=[f'{d.cond_name.strip()}-half{d.half}' for i,d in data_info.iterrows()]
        elif type == 'CondRun':

            # Subset of info sutructure
            data_info,C = agg_data(info,
                        ['run','cond_num'],
                        ['reg_num'],
                        subset=(info.error==0))
            data_info['names']=[f'{d.cond_name.strip()}-run{d.run:02d}' for i,d in data_info.iterrows()]
        elif type == 'CondAll':
            data_info,C = agg_data(info,
                        ['cond_num'],
                        ['run','half','reg_num'],
                        subset=(info.error==0))
            # data_info = info_gb.agg({'n_rep':np.sum}).reset_index(drop=True)
            data_info['names']=[f'{d.cond_name.strip()}' for i,d in data_info.iterrows()]

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir+f'/{participant_id}_{ses_id}_designmatrix.npy')
        reg_in = np.arange(C.shape[1],dtype=int)
        data_new = optimal_contrast(data_n,C,X,reg_in,baseline=None)

        return data_new, data_info


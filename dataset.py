#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data fusion project dataset module

The class for converting and mapping raw data from multi-dataset
to a standard data structure that can be used in Diedrichsen lab

"""
import numpy as np
import nibabel as nib
import pandas as pd
import os
import util
import matrix
import atlas_map as am
import scipy.linalg as sl
import nibabel as nb
import nitools as nt
from pathlib import Path
from numpy import eye
from numpy.linalg import pinv,solve

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
        data_n[i] = data_n[i] / np.sqrt(np.abs(resms))
    return data_n

def optimal_contrast(data,C,X,reg_in=None,baseline=None):
    """Recombines betas from a GLM into an optimal new contrast, taking into account a design matrix

    Args:
        data (list of ndarrays): List of N x P_i arrays of data
        C (ndarray): N x Q array indicating contrasts
        X (ndarray): Optimal design matrix - Defaults to None.
        reg_in (ndarray): Contrast of interest: Logical vector indicating which rows of C we will put in the matrix
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

    def get_data_fnames(self,participant_id,session_id=None):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dir = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T=pd.read_csv(dir+f'/{participant_id}_{session_id}_reginfo.tsv',sep='\t')
        fnames = [f'{dir}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i,t in T.iterrows()]
        fnames.append(f'{dir}/{participant_id}_{session_id}_resms.nii')
        return fnames, T


    def get_data(self, participant_id, atlas_maps):
        """the main function to output the processed data
        Args:
            participant_id: standard participant_id
            atlas_maps: List of atlasmaps to find the voxels

        Returns:
            Y (np.ndarray):
                A numatlasses list with N x P numpy array of data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
        """
        pass

    def get_all_suit(self,ses_id='ses-s1',type='CondSes',atlas='SUIT3'):
        """Extracts data in SUIT space from a standard experiment structure
        across all subjects. Saves the results as CIFTI files in the data directory.
        Args:
            ses_id (str, optional): Session. Defaults to 'ses-s1'.
            type (str, optional): Type - defined in ger_data. Defaults to 'CondSes'.
            atlas (str, optional): Short atlas string. Defaults to 'SUIT3'.
        """
        # Make the atlas object
        if (atlas=='SUIT3'):
            mask = self.atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
        if (atlas=='SUIT2'):
            mask = self.atlas_dir + '/tpl-SUIT/tpl-SUIT_res-2_gmcmask.nii'
        if (atlas=='MNISymC3'):
            mask = self.atlas_dir + '/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
        if (atlas=='MNISymC2'):
            mask = self.atlas_dir + '/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-2_gmcmask.nii'
        suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)

        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:
            print(f'Atlasmap {s}')
            deform = self.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
            if atlas[0:7]=='MNISymC':
                xfm_name = self.atlas_dir + '/tpl-MNI152NLIn2000cSymC/tpl-SUIT_space-MNI152NLin2009cSymC_xfm.nii'
                deform = [xfm_name,deform]
            mask = self.suit_dir.format(s) + f'/{s}_desc-cereb_mask.nii'
            atlas_map = am.AtlasMapDeform(self, suit_atlas, s,deform, mask)
            atlas_map.build(smooth=2.0)
            print(f'Extract {s}')
            data,info = self.get_data(s,[atlas_map],
                                                    ses_id=ses_id,
                                                    type=type)
            C=am.data_to_cifti(data,[atlas_map],info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir + f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv',sep='\t')

    def get_all_fs32k(self,ses_id='ses-s1',type='CondSes'):
        """Extracts data in fs32K space from a standard experiment structure
        across all subjects. Saves the results as CIFTI files in the data directory.

        Args:
            ses_id (str, optional): _description_. Defaults to 'ses-s1'.
            type (str, optional): _description_. Defaults to 'CondSes'.
        """
        # Make the atlas object
        atlas =[]
        bm_name = ['cortex_left','cortex_right']
        for i,hem in enumerate(['L','R']):
            mask = self.atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-{hem}_mask.label.gii'
            atlas.append(am.AtlasSurface(bm_name[i],mask_gii=mask))

        # create and calculate the atlas map for each participant
        T = self.get_participants()
        for s in T.participant_id:
            atlas_maps = []
            data = []
            for i,hem in enumerate(['L','R']):
                adir = self.anatomical_dir.format(s)
                edir = self.estimates_dir.format(s)
                pial = adir + f'/{s}_space-32k_hemi-{hem}_pial.surf.gii'
                white = adir + f'/{s}_space-32k_hemi-{hem}_white.surf.gii'
                mask = edir + f'/ses-s1/{s}_ses-s1_mask.nii'
                atlas_maps.append(am.AtlasMapSurf(self, atlas[i],
                            s,white,pial, mask))
                atlas_maps[i].build()
            print(f'Extract {s}')
            data,info = self.get_data(s,atlas_maps,
                                                ses_id=ses_id,
                                                type=type)
            C=am.data_to_cifti(data,atlas_maps,info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
            pass

class DataSetMDTB(DataSet):
    def __init__(self, dir):
        super().__init__(dir)


    def get_data(self,participant_id,
                    atlas_maps,
                    ses_id,
                    type='CondSes'):
        """ MDTB extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondSes': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondSes'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
        fnames,info = self.get_data_fnames(participant_id,ses_id)
        data = am.get_data3D(fnames,atlas_maps)
        # For debugging: data = [np.random.normal(0,1,(len(fnames),atlas_maps[0].P))]

        # Depending on the type, make a new contrast
        info['half']=2-(info.run<9)
        n_cond = np.max(info.cond_num)
        if type == 'CondSes':

            # Make new data frame for the information of the new regressors
            ii = ((info.run == 1) | (info.run == 9)) & (info.cond_num>0)
            data_info = info[ii].copy().reset_index()
            data_info['names']=[f'{d.cond_name.strip()}-sess{d.half}' for i,d in data_info.iterrows()]

            # Contrast for the regressors of interest
            reg = (info.half-1)*n_cond + info.cond_num
            reg[info.instruction==1] = 0
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.half*info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond*2,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'CondRun':

            # Subset of info sutructure
            ii = (info.cond_num>0)
            data_info = info[ii].copy().reset_index()
            data_info['names']=[f'{d.cond_name}-run{d.run:02d}' for i,d in data_info.iterrows()]

            reg = (info.run-1)*n_cond + info.cond_num
            reg[info.instruction==1] = 0
            # Contrast for the regressors of interst
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.run*info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond*16,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.run,positive=True)
        elif type == 'CondAll':

            # Make new data frame for the information of the new regressors
            ii = (info.run == 1)  & (info.cond_num>0)
            data_info = info[ii].copy().reset_index()
            data_info['names']=[f'{d.cond_name}' for i,d in data_info.iterrows()]

            # Contrast for the regressors of interest
            reg = info.cond_num.copy()
            reg[info.instruction==1] = 0
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.run,positive=True)

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir+f'/{participant_id}_{ses_id}_designmatrix.npy')
        data_new = optimal_contrast(data_n,C,X,reg_in,baseline=B)
        
        return data_new, data_info

class DataSetHcpResting(DataSet):
    def __init__(self, dir):
        super(DataSetHcpResting, self).__init__(base_dir=dir)
        # self.func_dir = self.base_dir + '/{0}/estimates'
        self.all_sub = self.get_participants()
        self.derivative_dir = self.base_dir + '/derivatives'

    def get_data_fnames(self, participant_id):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
        Returns:
            fnames (list): List of fnames
        """
        dir = self.derivative_dir + f"/{participant_id}" + "/func"
        fnames = []
        for r in range(4):
            fnames.append(f'{dir}/sub-{participant_id}_run-{r}_space-MSMSulc.dtseries.nii')
        return fnames

    def get_ts_volume(self,
                participant_id,
                atlas_map,
                runs=[0,1,2,3]):
        """ Returns the time series data for an atlas map
            sample from voxels
        Args:
            participant_id (_type_): _description_
            atlas_map (_type_): _description_
        """
        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        ts_volume = []
        for r in runs:
            # load the cifti
            ts_cifti = nb.load(fnames[r])

            # get the ts in volume for subcorticals
            ts_vol = util.volume_from_cifti(ts_cifti)
            # transfer to suit space applying the deformation
            ts_vol = am.get_data4D(ts_vol, [atlas_map])
            ts_volume.append(ts_vol[0])

        return ts_volume

    def get_ts_surface(self,
                    participant_id,
                    atlas_parcel,
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


    def get_cereb_connectivity(self,
                 participant_id,
                 cereb_atlas_map,
                 cortical_atlas_parcels,
                 runs=[0,1,2,3]):
        """
        Uses the original CIFTI files to produce cerebellar connectivity
        file
        """
        hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT', 'CIFTI_STRUCTURE_CORTEX_RIGHT']
        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef = None
        for r in runs:
            # load the cifti
            ts_cifti = nb.load(fnames[r])

            # get the ts in volume for subcorticals
            ts_vol = util.volume_from_cifti(ts_cifti)
            # transfer to suit space applying the deformation
            ts_cerebellum = am.get_data4D(ts_vol, [cereb_atlas_map])
            ts_cerebellum = ts_cerebellum[0]

            # get the ts in surface for corticals
            ts_32k = util.surf_from_cifti(ts_cifti,hem_name)
            ts_parcel = []
            for hem in range(2):

                # get the average within parcels
                ts_parcel.append(
                    cortical_atlas_parcels[hem].agg_data(ts_32k[hem])
                    )

            # concatenate them into a single array for correlation calculation
            ts_cortex = np.concatenate(ts_parcel, axis = 1)

            # Standardize the time series for easier calculation
            ts_cerebellum = util.zstandarize_ts(ts_cerebellum)
            ts_cortex = util.zstandarize_ts(ts_cortex)

            # Correlation calculation
            if coef is None:
                coef=np.empty((len(runs),ts_cortex.shape[1],
                                ts_cerebellum.shape[1]))
            N = ts_cerebellum.shape[0]
            coef[r,:,:] = ts_cortex.T @ ts_cerebellum / N

        return coef  # shape (n_tessl,P)

class DataSetPontine(DataSet):
    def __init__(self, dir):
        super().__init__(dir)

    def get_data(self, participant_id,
                 atlas_maps,
                 ses_id,
                 type='CondSes'):
        """ Pontine extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondSes': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondSes'.

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
        if type == 'TaskSes':

            # Make new data frame for the information of the new regressors
            ii = ((info.run == 1) | (info.run == 9)) & (info.reg_id > 0)
            data_info = info[ii].copy().reset_index()
            data_info['names'] = [
                f'{d.task_name.strip()}-sess{d.half}' for i, d in data_info.iterrows()]

            # Contrast for the regressors of interest
            reg = (info.half - 1) * n_cond + info.reg_id
            reg[info.instruction == 1] = 0
            C = matrix.indicator(reg, positive=True)  # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.half*info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond*2,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'TaskRun':

            # Subset of info sutructure
            ii = (info.reg_id > 0)
            data_info = info[ii].copy().reset_index()
            data_info['names'] = [
                f'{d.task_name.strip()}-run{d.run:02d}' for i, d in data_info.iterrows()]

            reg = (info.run - 1) * n_cond + info.reg_id
            reg[info.instruction == 1] = 0
            reg = (info.run-1)*n_cond + info.cond_num
            reg[info.instruction==1] = 0
            # Contrast for the regressors of interst
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.run*info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond*16,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.run,positive=True)

        elif type == 'TaskAll':

            # Make new data frame for the information of the new regressors
            ii = (info.run == 1) & (info.reg_id > 0)
            data_info = info[ii].copy().reset_index()
            data_info['names'] = [
                f'{d.task_name.strip()}' for i, d in data_info.iterrows()]

           # Contrast for the regressors of interest
            reg = info.cond_num.copy()
            reg[info.instruction==1] = 0
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            CI = matrix.indicator(info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.run,positive=True)


        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
        data_new = optimal_contrast(data_n, C, X, reg_in)

        return data_new, data_info


class DataSetNishi(DataSet):
    def __init__(self, dir):
        super().__init__(dir)


    def get_data(self,participant_id,
                    atlas_maps,
                    ses_id,
                    type='CondSes'):
        """ MDTB extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondSes': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondSes'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
        fnames,info = self.get_data_fnames(participant_id,ses_id)
        data = am.get_data3D(fnames,atlas_maps)
        # For debugging: data = [np.random.normal(0,1,(len(fnames),atlas_maps[0].P))]

        # Depending on the type, make a new contrast
        info['half']=2-(info.run< (len(np.unique(info.run))/2+1))
        n_cond = np.max(info.reg_id)
        if type == 'CondSes':
            # Make new data frame for the information of the new regressors            
            ii = ((info.run == min(info.run)) | (info.run == len(np.unique(info.run))/2+1)) 
            data_info = info[ii].copy().reset_index()
            data_info['names']=[f'{d.task_name}-sess{d.half}' for i,d in data_info.iterrows()]

            # Contrast for the regressors of interest
            reg = (info.half-1)*n_cond + info.reg_id
            # reg[info.instruction==1] = 0
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            # CI = matrix.indicator(info.half,positive=True)
            # C = np.c_[C,CI]
            reg_in = np.arange(n_cond*2,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.half,positive=True)

        elif type == 'CondRun':

            # Subset of info sutructure
            # ii = (info.cond_num>0)
            data_info = info.copy().reset_index()
            data_info['names']=[f'{d.task_name}-run{d.run:02d}' for i,d in data_info.iterrows()]

            reg = (info.run-1)*n_cond + info.reg_id
            # reg[info.instruction==1] = 0
            # Contrast for the regressors of interst
            C = matrix.indicator(reg,positive=True) # Drop the instructions


            # contrast for all instructions
            # CI = matrix.indicator(info.run*info.instruction,positive=True)
            # C = np.c_[C,CI]
            reg_in = np.arange(n_cond*len(np.unique(info.run)),dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.run,positive=True)
        elif type == 'CondAll':

            # Make new data frame for the information of the new regressors
            ii = (info.run == min(info.run))
            data_info = info[ii].copy().reset_index()
            data_info['names']=[f'{d.task_name}' for i,d in data_info.iterrows()]

            # Contrast for the regressors of interest
            reg = info.creg_id
            # reg[info.instruction==1] = 0
            C = matrix.indicator(reg,positive=True) # Drop the instructions

            # contrast for all instructions
            # CI = matrix.indicator(info.instruction,positive=True)
            # C = np.c_[C,CI]
            # reg_in = np.arange(n_cond-1,dtype=int)

            # Baseline substraction 
            B = matrix.indicator(data_info.run,positive=True)

        # Prewhiten the data
        data_n = prewhiten_data(data)

        # Load the designmatrix and perform optimal contrast
        X = np.load(dir+f'/{participant_id}_{ses_id}_designmatrix.npy')
        data_new = optimal_contrast(data_n,C,X,reg_in,baseline=B)
        
        return data_new, data_info

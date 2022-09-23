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

class DataSetMDTB(DataSet):
    def __init__(self, dir):
        super().__init__(dir)


    def get_data_fnames(self,participant_id,session_id=None):
        """ Gets all raw data files

        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dir = self.estimates_dir.format(participant_id) + f'/{session_id}'
        T=pd.read_csv(dir+f'/{participant_id}_{session_id}_reginfo.tsv',sep='\t')
        fnames = [f'{dir}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i,t in T.iterrows()]
        fnames.append(f'{dir}/{participant_id}_{session_id}_resms.nii')
        return fnames, T

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
        # data = np.random.normal(0,1,(737,atlas_maps[0].P))
        data = am.get_data3D(fnames,atlas_maps)
        # Load design matrix for optimal reweighting
        X = np.load(dir+f'/{participant_id}_{ses_id}_designmatrix.npy')

        # determine the different halfs
        info['half']=2-(info.run<9)

        if type == 'CondSes':
            n_cond = np.max(info.cond_num)
            reg = (info.half-1)*n_cond + info.cond_num
            reg[info.instruction==1] = 0
            # Contrast for the regressors of interst
            C = matrix.indicator(reg,positive=True) # Drop the instructions
            # contrast for all instructions
            CI = matrix.indicator(info.half*info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond*2,dtype=int)
            # Subset of info sutructire
            ii = ((info.run == 1) | (info.run == 9)) & (info.cond_num>0)
            data_info = info[ii].copy().reset_index()
            names=[f'{d.cond_name}-sess{d.half}' for i,d in data_info.iterrows()]
        elif type == 'CondRun':
            n_cond = np.max(info.cond_num)
            reg = (info.run-1)*n_cond + info.cond_num
            reg[info.instruction==1] = 0
            # Contrast for the regressors of interst
            C = matrix.indicator(reg,positive=True) # Drop the instructions
            # contrast for all instructions
            CI = matrix.indicator(info.run*info.instruction,positive=True)
            C = np.c_[C,CI]
            reg_in = np.arange(n_cond*16,dtype=int)
            # Subset of info sutructire
            ii = (info.cond_num>0)
            data_info = info[ii].copy().reset_index()
            names=[f'{d.cond_name}-run{d.run:02d}' for i,d in data_info.iterrows()]
        elif type == 'CondAll':
            pass
        pass
        # Add the block regressors
        C = sl.block_diag(C,np.eye(16))
        Xn = X @ C

        # Get the resms and prewhiten the data
        data_n = []
        for i in range(len(data)):
            # Prewhiten the data univariately
            resms = data[i][-1,:]
            data[i] = data[i][0:-1,:]
            data[i] = data[i] / np.sqrt(np.abs(resms))
            # Append the intercept regressors
            data[i] = np.concatenate([data[i],np.zeros((16,data[i].shape[1]))])
            # Do the averaging / reweighting:
            d = np.linalg.solve(Xn.T @ Xn, Xn.T @ X @ data[i])
            # Put the data in the list
            data_n.append(d[reg_in,:])
        return data_n, data_info, names

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


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
                    sess_id,
                    type='CondSes'):
        dir = self.estimates_dir.format(participant_id) + f'/{sess_id}'
        fnames,info = self.get_data_fnames(participant_id,sess_id)
        # data = np.random.normal(0,1,(737,atlas_maps[0].P))
        data = am.get_data3D(fnames,atlas_maps)
        # Load design matrix for optimal reweighting
        X = np.load(dir+f'/{participant_id}_{sess_id}_designmatrix.npy')

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
        elif type == 'TaskSes':
            pass
        elif type == 'AllSes':
            pass
        pass
        # Add the block regressors
        C = sl.block_diag(C,np.eye(16))
        Xn = X @ C

        # get the resms and prewhiten the data
        data_n = []
        for i in range(len(data)):
            resms = data[i][-1,:]
            data[i] = data[i][0:-1,:]
            data[i] = data[i] / np.sqrt(np.abs(resms))
            # Append the intercept regressors
            data[i] = np.concatenate([data[i],np.zeros((16,data[i].shape[1]))])

            d = np.linalg.solve(Xn.T @ Xn, Xn.T @ X @ data[i])
            data_n.append(d[reg_in,:])
        return data_n, data_info, names


class DataSetHcpResting(DataSet):
    def __init__(self, dir):
        super(DataSetHcpResting, self).__init__(base_dir=dir)
        # self.func_dir = self.base_dir + '/{0}/estimates'
        self.all_sub = self.get_participants()
        self.derivative_dir = self.base_dir + '/derivatives'
        self.func_dir = self.base_dir + '/func'

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
                atlas_map): 
        """ Returns the time series data for an atlas map 
            sample from voxels
        Args:
            participant_id (_type_): _description_
            atlas_map (_type_): _description_
        """
        pass
    # def get_data_vol(self, participant_id,
    #              atlas_maps,
    #              sess_id=None,
    #              ):
    #     # print(f"I'm in get_data")
    #     # get the file name for the cifti time series
    #     fnames = self.get_data_fnames(participant_id,sess_id)

    #     # get the volumetric data for subcorticals
    #     vol_ts = []
    #     for file_name in fnames:
    #         # for subcorticals
    #         vol_4D = self._volume_from_cifti(file_name)



    #         # get cerebellar time series in suit space
    #         vol_ts.append(am.get_data4D(vol_4D,atlas_maps))

    #     return vol_ts

    def get_ts_surface(self,
                    participant_id,
                    atlas_parcel):
        """Returns the information from the CIFTI file
        in the 32K surface for left and right hemisphere. 

        Args:
            participant_id (_type_): _description_
            atlas_parcel (_type_): _description_
        """
        pass 
    # def get_data_surf(self, participant_id,
    #              atlas_maps,
    #              sess_id=None,
    #              ):
    #     # get the file name for the cifti time series
    #     fnames = self.get_data_fnames(participant_id,sess_id)

    #     # get the volumetric data for subcorticals
    #     surf_ts = []
    #     for file_name in fnames:

    #         # for cortical hemispheres
    #         surf_ts.append(self._surf_from_cifti(file_name))

    #     return surf_ts



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
            for hem in range(2):

                # get the average within parcels 
                ts_parcel.append(am.get_average_data(data=ts_32k,
                                                     labels=label_objs[hemi],
                                                     mask=mask_obj[hemi]))

            # concatenate them into a single array for correlation calculation
            ts_cortex = np.concatenate(ts_parcel, axis = 1)

            # calculate functional connectivity matrix
            ## no need to define new variables, but still ...
            ts_cerebellum = util.zstandarize_ts(ts_cerebellum)
            ts_cortex = util.zstandarize_ts(ts_cortex)

        return np.concatenate((coef[0], coef[1]), axis=1)  # shape (P, 2 * n_tessl)


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
            data_n: Univariately prewhitened data
            data_info: Pandas data frame with information
            names: Names for CIFTI-file per roww
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
            reg_in = np.arange(n_cond*2,dtype=int)
            # Subset of info sutructire
            ii = (info.cond_num>0)
            data_info = info[ii].copy().reset_index()
            names=[f'{d.cond_name}-run{d.half}' for i,d in data_info.iterrows()]
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

    def _volume_from_cifti(self, ts_cifti, save=False):
        """
        Gets the 4D nifti object containing the time series
        for all the subcortical structures
        Args:
            ts_cifti (cifti obj ) - cifti object of the time series
            save (Boolean) - set to True if you want to save the 4D nifti (NOT IMPLEMENTED YET)
        Returns:
            nii_vol(nifti vol object) - nifti object containing the time series of subcorticals
        """
        # get brain axis models
        bmf = ts_cifti.header.get_axis(1)
        # get the data array with all the time points, all the structures
        ts_array = ts_cifti.get_fdata(dtype=np.float32)

        # initialize a matrix representing 4D data (x, y, z, time_point)
        subcorticals_vol = np.zeros(bmf.volume_shape + (ts_array.shape[0],))
        for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
            # get the values corresponding to the brain model
            bm_vals = ts_array[:, slc]

            # get the voxels/vertices corresponding to the current brain model
            ijk = bm.voxel
            # fill in data
            if (idx != 0) & (idx != 1): # indices 0 and 1 are cortical hemispheres
                # print(str(nam))
                subcorticals_vol[ijk[:, 0], ijk[:, 1], ijk[:, 2], :] = bm_vals.T 

        # save as nii 
        nii_vol_4d = nb.Nifti1Image(subcorticals_vol,bmf.affine)
        # if save:
        #     ts_nifti = dir+'/sub-100307_ses-01_task-rest_space-subcortex_run-01_bold.nii'
        #     nb.save(nii_vol,ts_nifti)
        return nii_vol_4d

    def _surf_from_cifti(self, ts_cifti, save = False):
        """
        Gets the time series of cortical surface vertices (Left and Right)
        Args:
            ts_cifti (cifti obj) - cifti object of time series
            save (Boolean) - set to True if you want to save the 4D nifti (NOT IMPLEMENTED YET)
        Returns:
            cii (cifti object) - contains the time series for the cortex
        """
        # get brain axis models
        bmf = ts_cifti.header.get_axis(1)
        # print(dir(bmf))
        # get the data array with all the time points, all the structures
        ts_array = ts_cifti.get_fdata()
        ts_list = []
        for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
            # just get the cortical surfaces
            if (idx == 0) | (idx == 1): # just for cortex left (idx = 0) and cortex right (idx = 1)
                # get the values corresponding to the brain model
                bm_vals = ts_array[:, slc]
                ts_list.append(bm_vals)
            else:
                break
        return ts_list

    def get_data_fnames(self, participant_id, session_id=None):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
        Returns:
            fnames (list): List of fnames
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dir = self.func_dir.format(f"{participant_id}") + f'/{session_id}'
        fnames = [f'{dir}/sub-{participant_id}_{session_id}_space-fsLR32k_run-{r:02}.dtseries.nii' for r in [1, 2]]
        return fnames

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

    def get_data(self, participant_id, 
                 atlas_maps, 
                 label_objs, 
                 mask_obj,
                 sess_id = None):
        """
        atlas maps only needed for the cerebellum!
        """

        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id,sess_id)
        coef = []
        for f in fnames:
            # load the cifti 
            ts_cifti = nb.load(f)

            # get the ts in volume for subcorticals
            ts_vol = self._volume_from_cifti(ts_cifti)
            # transfer to suit space applying the deformation
            ts_cerebellum = am.get_data4D(ts_vol, atlas_maps)

            # get the ts in surface for corticals
            ts_32k = self._surf_from_cifti(ts_cifti, save=False)

            # get the average within parcels for left and right hemi
            ts_parcel = []
            for hemi in [0, 1]:
                ts_parcel.append(am.get_average_data(data=ts_32k[hemi],
                                                     labels=label_objs[hemi],
                                                     mask=mask_obj[hemi]))
            
            # concatenate them into a single array for correlation calculation
            ts_cortex = np.concatenate(ts_parcel, axis = 1)

            # calculate functional connectivity matrix
            ## no need to define new variables, but still ... 
            Y = np.asarray(ts_cerebellum.copy()[0])
            X = ts_cortex.copy()

            # subtract mean
            X = X - X.mean(axis = 0, keepdims = True)
            Y = Y - Y.mean(axis = 0, keepdims = True)

            # X is now shape of (1200, n_corticaltess) ndarray
            # Y is now shape of (1200, n_suit_voxel) ndarray
            a = np.nansum((X - np.nanmean(X, axis=0)) ** 2, axis=0)
            b = np.nansum((Y - np.nanmean(Y, axis=0)) ** 2, axis=0)
            var = np.sqrt(b.reshape(-1, 1) @ a.reshape(1, -1))
            cov = Y.T @ X
            coef.append(cov/var)

        return np.concatenate((coef[0], coef[1]), axis=1)  # shape (P, 2 * n_tessl)

if __name__ == '__main__':
    A = DataSetHcpResting('Y:\data\FunctionalFusion\HCP')
    part = A.get_participants()
    data = A.get_data(part['participant_id'][0], sess_id="01")

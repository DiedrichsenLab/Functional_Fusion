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
        data = am.get_data(fnames,atlas_maps)
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

    def _volume_from_cifti(self, data, axis, vol_name):
        assert isinstance(axis, nib.cifti2.BrainModelAxis)
        vol_data = np.full(axis.volume_shape + (data.shape[0],), np.nan)
        for name, data_indices, model in axis.iter_structures():
            if name in vol_name:
                data = data.T[data_indices]  # Assume brainmodels axis is last, move it to front
                vox_indices = tuple(model.voxel.T)
                vol_data[vox_indices] = data  # "Fancy indexing"

        N = nib.Nifti1Image(vol_data, axis.affine)
        return N, vol_data  # Add affine for spatial interpretation

    def _surf_data_from_cifti(self, data, axis, surf_name):
        assert isinstance(axis, nib.cifti2.BrainModelAxis)
        for name, data_indices, model in axis.iter_structures():
            if name == surf_name:
                data = data.T[data_indices]  # Assume brainmodels axis is last, move it to front
                vtx_indices = model.vertex  # Generally 1-N, except medial wall vertices
                surf_data = np.zeros((vtx_indices.max() + 1,) + data.shape[1:], dtype=data.dtype)
                surf_data[vtx_indices] = data
                return surf_data
        raise ValueError(f"No structure named {surf_name}")

    def get_data_fnames(self, participant_id, session_id=None):
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

    def get_data(self, participant_id, atlas_map=None, sess_id=None,
                 run=None, roi=None):
        if run is None:
            run = ['01']
        if sess_id is None:
            sess_id = ['01']
        if roi is None:
            roi = ['CEREBELLUM_LEFT']

        dir = self.estimates_dir.format(participant_id) + f'/{sess_id}'
        fnames, info = self.get_data_fnames(participant_id, sess_id)

        func_dir = os.path.join(self.derivative_dir, f'{participant_id}/estimates')
        for s in sess_id:
            for r in run:
                file = 'sub-%d_ses-%s_task-rest_space-fsLR32k_run-%s_bold.nii' % \
                       (participant_id, s, r)
                G = nib.load(os.path.join(func_dir, file))
                axes = [G.header.get_axis(i) for i in range(G.ndim)]
                cifti_data = G.get_fdata(dtype=np.float32)
                # data = self._surf_data_from_cifti(cifti_data, axes[1], 'CIFTI_STRUCTURE_'+roi)
                this_roi = ['CIFTI_STRUCTURE_'+x for x in roi]
                _, data = self._volume_from_cifti(cifti_data, axes[1], this_roi)


if __name__ == '__main__':
    A = DataSetHcpResting('Y:\data\FunctionalFusion\HCP')
    part = A.get_participants()
    data = A.get_data(part['participant_id'][0], sess_id="ses-s1")
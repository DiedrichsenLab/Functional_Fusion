#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data fusion project dataset module

The class for converting and mapping raw data from multi-dataset
to a standard data structure that can be used in Diedrichsen lab

Created on 3/30/2022 at 12:21 PM
Author: dzhi
"""
import numpy as np
import torch as pt
import nibabel as nib
import pandas as pd
import os


class DataSet:
    def __init__(self, base_dir):
        """DataSet class:
        Implements the interface for each of the data set
        Note that the actual preprocessing and glm estimate
        do not have to be performed with functionality provided by
        this class. The class is just a instrument to present the user with
        a uniform interface of how to get subject info

        Args:
            basedir (str): _description_
        """
        self.base_dir  = base_dir
        self.surface_dir = base_dir + '/{0}/surface'
        self.anatomical_dir = base_dir + '/{0}/anatomical'
        self.contrast_dir = base_dir + '/{0}/contrast'
        self.suit_dir = base_dir + '/{0}/suit'

    def get_participants(self):
        """ returns a data frame with all participants
        available in the study. The fields in the data frame correspond to the
        standard columns in participant.tsv.
        https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
        Returns:
            Pinfo (pandas data frame): participant information in standard bids format
        """
        Pinfo = pd.read_csv(self.base_dir + '/participants.tsv')
        return Pinfo

    def get_data(self, participant_id, atlas_map):
        """the main function to output the processed data
        Args:
            participant_id: standard participant_id
            atlas_map: AtlasMAP to find the voxels

        Returns:
            Y (np.ndarray):
                A N x P numpy array of aggregated data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
        """
        pass


class DataSet_HCP_resting(DataSet):
    def __init__(self, dir='Y:\data\FunctionalFusion\HCP'):
        super(DataSet_HCP_resting, self).__init__(base_dir=dir)
        # self.func_dir = self.base_dir + '/{0}/estimates'
        self.all_sub = self.get_participants()
        self.derivative_dir = self.base_dir + '/derivatives'

    def _volume_from_cifti(self, data, axis, vol_name):
        assert isinstance(axis, nib.cifti2.BrainModelAxis)
        for name, data_indices, model in axis.iter_structures():
            if name == vol_name:
                data = data.T[data_indices]  # Assume brainmodels axis is last, move it to front
                vox_indices = tuple(model.voxel.T)
                vol_data = np.full(axis.volume_shape + data.shape[1:], np.nan)
                vol_data[vox_indices] = data  # "Fancy indexing"
                return vol_data  # Add affine for spatial interpretation

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

    def get_data(self, participant_id, atlas_map=None, session=None,
                 run=None, roi='CEREBELLUM_LEFT'):
        if run is None:
            run = ['01']
        if session is None:
            session = ['01']

        for sub in participant_id:
            func_dir = os.path.join(self.derivative_dir, f'{sub}/estimates')
            for s in session:
                for r in run:
                    file = 'sub-%d_ses-%s_task-rest_space-fsLR32k_run-%s_dtseries.nii' % (sub, s, r)
                    G = nib.load(os.path.join(func_dir, file))
                    axes = [G.header.get_axis(i) for i in range(G.ndim)]
                    cifti_data = G.get_fdata(dtype=np.float32)
                    # data = self._surf_data_from_cifti(cifti_data, axes[1], 'CIFTI_STRUCTURE_'+roi)
                    data = self._volume_from_cifti(cifti_data, axes[1], 'CIFTI_STRUCTURE_'+roi)


if __name__ == '__main__':
    A = DataSet_HCP_resting()
    part = A.get_participants()
    data = A.get_data(part['participant_id'][0:10])


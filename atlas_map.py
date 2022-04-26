#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The functions of atlas definition and atlas mapping

Created on 3/30/2022 at 3:00 PM
Author: dzhi, jdiedrichsen
"""
import numpy as np
import nibabel as nb
import os


class Atlas():
    def __init__(self):
        """The Atlas class implements the general atlas functions
        for mapping from the P brain locations back to nii or gifti files
        Each Atlas is associated with a set of atlas maps
        """
        self.id = 'unknown'
        self.P = np.nan # Number of locations in this atlas

    def map_data(self,data):
        """Maps data back into some atlas form.
        Args:
            data (numpy.ndarray): P or N x P array
        """
        pass

class AtlasVolumetric(Atlas):
    def __init__(self,id,mask_img):
        self.id = id
        self.mask_img = nb.load(mask_img)
        Xmask = self.mask_img.get_data()
        Xmask = (Xmask>0)
        i,j,k = np.where(X>0)
        # self.x,self.y,self.z = affine_transform 



class AtlasMap:
    def __int__(self, dataset, atlas, participant_id):
        """AtlasMap stores the mapping rules from a specific data set (and participant)
        to the desired atlas space in form of a voxel list
        Args:
            dataset_id (string): name of
            participant_id (string): Participant name
        """
        self.participant_id = participant_id
        self.dataset = dataset # Reference to corresponding data set
        self.atlas = atlas     # Reference to corresponding altas
        self.P = atlas.P       #  Number of brain locations

    def build(self):
        """
        Using the dataset, build creates a list of voxel indices of
        For each of the locations, it
        """
        self.voxel_list = [np.empty((1,),dtype=int)]*self.P

    def save(self, file_name):
        """serializes a atlas map to a file
        Args:
            file_name (string): full filepath and name
        """
        pass

    def load(self, file_name):
        """loads a build atlas map from file

        Args:
            file_name (string): full filepath and name
        """
        pass


class AtlasMapDeform(AtlasMap):
    def __int__(self, dataset, atlas, participant_id):
        """AtlasMap stores the mapping rules from a specific data set (and participant)
        to the desired atlas space in form of a voxel list
        Args:
            dataset_id (string): name of
            participant_id (string): Participant name
        """
        super().__init__(participant_id,dataset,atlas)
        self.id = atlas.id # Name of atlas space 


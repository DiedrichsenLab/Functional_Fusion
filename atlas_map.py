#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The functions of atlas definition and atlas mapping

Created on 3/30/2022 at 3:00 PM
Author: dzhi, jdiedrichsen
"""
import numpy as np
from numpy.linalg import inv
import nibabel as nb
import os
import SUITPy as suit
from SUITPy.flatmap import affine_transform
import surfAnalysisPy as surf

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
        """Atlas Volumetric class constructor 

        Args:
            id (str): Name of the altas (e.g. SUIT2)
            mask_img (str): file name of mask image defining atlas location
        """
        self.id = id
        self.mask_img = nb.load(mask_img)
        Xmask = self.mask_img.get_data()
        Xmask = (Xmask>0)
        self.i,self.j,self.k = np.where(Xmask>0)
        self.x,self.y,self.z = affine_transform(self.i,self.j,self.k,self.mask_img.affine)
        self.P = self.x.shape[0]

    def map_data(self,data):
        """Maps data back into a full nifti

        Args:
            data (ndarray): 1-d Numpy array of the size (P,)

        Returns:
            mapped_image (Nifti1Image): Image containing mapped results 
        """
        X=np.zeros(self.mask_img.shape)
        X[self.i,self.j,self.k]=data
        mapped = nb.Nifti1Image(X,self.mask_img.affine)
        return mapped

class AtlasMap():
    def __init__(self, dataset, atlas, participant_id):
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
    def __init__(self, dataset, atlas, participant_id, deform_img,mask_img):
        """AtlasMap stores the mapping rules from a specific data set (and participant)
        to the desired atlas space in form of a voxel list
        Args:
            dataset_id (str): name of
            participant_id (str): Participant name
            deform_img (str): Name for deformation map image
            mask_img (str): Name of masking image that defines the functional data space. 
        """
        super().__init__(dataset,atlas,participant_id)
        self.id = atlas.id
        self.deform_img = nb.load(deform_img)
        self.mask_img = nb.load(mask_img)
    
    def build(self,smooth = None):
        """
        Using the dataset, build creates a list of voxel indices of
        For each of the locations, it
        """
        self.voxel_list = np.empty((self.P,),dtype=object)
        i,j,k=affine_transform(self.atlas.x,self.atlas.y,self.atlas.z,
                            inv(self.deform_img.affine))
        i=i.astype(int)
        j=j.astype(int)
        k=k.astype(int)
        X=self.deform_img.get_fdata()
        
        pass
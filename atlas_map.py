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
import surfAnalysisPy as surf
from util import affine_transform, affine_transform_mat, sample_img_nn,coords_to_linvidxs, sq_eucl_distances

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
        i,j,k = np.where(Xmask>0)
        self.vox = np.vstack((i,j,k))
        self.world = affine_transform_mat(self.vox,self.mask_img.affine)
        self.P = self.world.shape[0]

    def map_data(self,data):
        """Maps data back into a full nifti

        Args:
            data (ndarray): 1-d Numpy array of the size (P,)

        Returns:
            mapped_image (Nifti1Image): Image containing mapped results 
        """
        X=np.zeros(self.mask_img.shape)
        X[self.vox[0],self.vox[1],self.vox[2]]=data
        mapped = nb.Nifti1Image(X,self.mask_img.affine)
        return mapped

class AtlasSurface(Atlas):
    def __init__(self,id,mask_gii):
        """Atlas Surface class constructor 

        Args:
            id (str): Name of the altas (e.g. SUIT2)
            mask_gii (str): gifti file name of mask image defining atlas locations
        """
        self.id = id
        self.mask_gii = nb.load(mask_gii)
        Xmask = self.mask_gii.agg_data()
        Xmask = (Xmask>0)
        self.vertex = np.nonzero(Xmask>0)[0]
        self.P = self.vertex.shape[0]

    def map_data(self,data):
        """Maps data back into a full nifti

        Args:
            data (ndarray): 1-d Numpy array of the size (P,)

        Returns:
            mapped_image (Nifti1Image): Image containing mapped results 
        """
        X=np.zeros(self.mask_img.shape)
        X[self.vox[0],self.vox[1],self.vox[2]]=data
        mapped = nb.Nifti1Image(X,self.mask_img.affine)
        return mapped


class AtlasMap():
    def __init__(self, dataset, atlas, participant_id):
        """AtlasMap stores the mapping rules from a specific data set (and participant) to the desired atlas space in form of a voxel list
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
        pass 

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
        """AtlasMapDeform stores the mapping rules for a non-linear deformation
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
        # Caluculate locations of atlas in individual (deformed) coordinates 
        atlas_ind = sample_img_nn(self.deform_img,self.atlas.world).squeeze().T
        N = atlas_ind.shape[1] # Number of locations in atlas
        if smooth is None: # Use nearest neighbor interpolation  
            self.vox_list,self.vox_weight = coords_to_linvidxs(atlas_ind,self.mask_img,mask=True)
        else:              # Use smoothing kernel of specific size 
            # Get world coordinates and linear coordinates for all available voxels
            M = self.mask_img.get_fdata()
            i,j,k=np.where(M>0)
            world_vox = affine_transform_mat(np.vstack((i,j,k)),self.mask_img.affine) # available voxels in world coordiantes
            linindx = np.ravel_multi_index((i,j,k),M.shape,mode='clip')
            
            # Distances between atlas coordinates and voxel coordinates 
            D = sq_eucl_distances(atlas_ind,world_vox)
            # Find voxels with substantial power under gaussian kernel 
            W = np.exp(-0.5 * D/(smooth**2))
            W[W<0.2]=0
            a,b=W.nonzero()
            # Now transfer them into a full list of voxels 
            # this is somewhat ugly and brute force 
            c = np.zeros(a.shape,dtype=int)
            c[1:]=a[0:-1]==a[1:]
            for i in range(c.shape[0]-1):
                 if c[i+1]:
                    c[i+1]=c[i+1]+c[i]
            self.vox_list=np.zeros((N,c.max()+1),dtype=np.int32)
            self.vox_weight=np.zeros((N,c.max()+1))
            self.vox_list[a,c]=linindx[b]
            self.vox_weight[a,c]=W[a,b]
            self.vox_weight = self.vox_weight / self.vox_weight.sum(axis=1,     keepdims=True)
        pass

class AtlasMapSurf(AtlasMap):
    def __init__(self, dataset, atlas, participant_id, 
                white_surf,pial_surf,mask_img):
        """AtlasMapSurf stores the mapping rules for a freesurfer-style surface (pial + white surface pair)
        Args:
            dataset_id (str): name of
            participant_id (str): Participant name
            white_surf (str): Name for white matter surface
            pial_surf (str): Name for pial surface
            mask_img (str): Name of masking image that defines the functional data space. 
        """
        super().__init__(dataset,atlas,participant_id)
        self.id = atlas.id
        self.white_surf = nb.load(white_surf)
        self.pial_surf = nb.load(pial_surf)
        self.mask_img = nb.load(mask_img)
    
    def build(self,smooth = None, depths=[0,0.2,0.4,0.6,0.8,1.0]):
        """
        Using the dataset, build creates a list of voxel indices of
        each of the nodes
        """
        n_points = len(depths)
        c1 = self.white_surf.darrays[0].data.T
        c2 = self.pial_surf.darrays[0].data.T
        n_vert = c1.shape[1]
        if c2.shape[1] != n_vert:
            raise(NameError('White and pial surfaces should have same number of vertices.'))

        # Get the indices for all the points being sampled
        indices = np.zeros((3,n_vert,n_points))
        for i in range(n_points):
            indices[:,:,i] = (1-depths[i])*c1+depths[i]*c2

        self.vox_list,good = coords_to_linvidxs(indices,self.mask_img,mask=True)
        self.vox_weight = good / good.sum(axis=1,keepdim=True)

def get_data(fnames,atlas_maps): 
    """Extracts the data for a list of fnames
    for a list of atlas_maps. This is usually called by DataSet.get_data() 
    to extract the required raw data before processing it further

    Args:
        fnames (list): list of file names to be sampled 
        atlas_maps (list): list of built atlas-map objects 
    """
    n_atlas = len(atlas_maps)
    n_files = len(fnames)
    data = []
    # Make the empty data structures
    for at in atlas_maps:
        data.append(np.full((n_files,at.vox_list.shape[0]),np.nan))
    for j,f in enumerate(fnames): 
        V = nb.load(f)
        X = V.get_fdata()
        if (X.ndim>3):
            raise(NameError('extraction right now only for 3d-niftis'))
        # Map this file into the data structures 
        X = X.ravel()
        for i,at in enumerate(atlas_maps): 
            d=X[at.vox_list] * at.vox_weight  # Expanded data
            d = np.nansum(d,axis=1)
            d[np.nansum(at.vox_weight,axis=1)==0]=np.nan
            data[i][j,:]=d
    return data
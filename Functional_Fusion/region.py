#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Region of Interest (ROI) classes 

Created on 3/30/2022 at 3:00 PM
Author: jdiedrichsen
"""
import numpy as np
from numpy.linalg import inv
import nibabel as nb
import warnings
import nitools as nt
import Functional_Fusion.atlas_map as am

class Region():
    def __init__(self):
        self.atlas = None
        self.locations = None
        self.P = 0 
        pass 

class RegionVolumetric(Region):
    def __init__(self,atlas):
        sle
        if isinstance(atlas,am.AtlasVolumetric):
            self.atlas = atlas
        else:
            raise ValueError('atlas must be of type AtlasVolumetric')


    def get_indiv_map(self,deform=None):
        """Returns the mapping from the region to the individual space
        Args:
            deform (Nifti1Image,list): Nii or list of nii that describe the deformation map. If list, the deformations are applied in order. 
            If empty, the identity map is assumed.
        Returns:
            atlas_map (am.AtlasMap): Individual Atlas map 
        """
        # Not implemented yet
        raise(NotImplementedError)

    def data_to_nifti(data):
        raise(NotImplementedError)
        


class RegionVolumetricSphere(RegionVolumetric):
    def __init__(self,atlas,coord,radius):
        self.atlas = atlas
        self.voxel = voxel


class RegionVolumetricImage(RegionVolumetric):
    def __init__(self,atlas,image,value=1):
        raise(NotImplementedError)

class RegionSurface(Region):
    def __init__(self,atlas,vertex):
        if isinstance(atlas,am.AtlasSurface):
            self.atlas = atlas
        else:
            raise ValueError('atlas must be of type AtlasSurface')

    def get_indiv_map(self,white,pial):
        """Returns the mapping from the region to the individual space
        Returns:
            atas_map (am.AtlasMap): Individual Atlas map 
        """
        raise(NotImplementedError)

    def region_to_cifti(): 
        raise(NotImplementedError)


class RegionSurfaceCircle(RegionSurface):
    def __init__(self,atlas,vertex,radius):
        """ Defines a speherical ROI on a group surface"""
        self.atlas = atlas
        self.vertex = vertex
        self.radius = radius
        raise(NotImplementedError)

class RegionSurfaceImage(RegionSurface):
    def __init__(self,atlas,img,hem=0,value=1):
        self.__super__(atlas)
        if isinstance(img, str):
            img = nb.load(img)
        if isinstance(img, nb.Cifti2Image):
            data = self.atlas.cifti_to_data(img)
        elif isinstance(img, nb.GiftiImage):
            d = img.agg_data()
            data = d[self.vertex[hem]]
    
        self.locations = np.where(data==value)[0]


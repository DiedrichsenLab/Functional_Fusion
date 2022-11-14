#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The functions of atlas definition and atlas mapping

Created on 3/30/2022 at 3:00 PM
Author: dzhi, jdiedrichsen
"""
from matplotlib.ticker import IndexLocator
import numpy as np
from numpy.linalg import inv
import nibabel as nb
import os
import warnings

import Functional_Fusion.matrix as matrix
import SUITPy as suit
import surfAnalysisPy as surf
import nitools as nt

def get_atlas(atlas_str,atlas_dir):
    """ returns an atlas from a code

    Args:
        atlas_str (str): Name of the atlas
        atlas_dir (str): directory name for the atlas
    """
    # Make the atlas object
    if (atlas_str=='SUIT3'):
        mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
        atlas = AtlasVolumetric('cerebellum',mask_img=mask)
    elif (atlas_str=='SUIT2'):
        mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-2_gmcmask.nii'
        atlas = AtlasVolumetric('cerebellum',mask_img=mask)
    elif (atlas_str=='MNISymC3'):
        mask = atlas_dir + '/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
        atlas = AtlasVolumetric('cerebellum',mask_img=mask)
    elif (atlas_str =='MNISymC2'):
        mask = atlas_dir + '/tpl-MNI152NLin2000cSymC/tpl-MNISymC_res-2_gmcmask.nii'
        atlas = AtlasVolumetric('cerebellum',mask_img=mask)
    elif atlas_str == 'fs32k':
        bm_name = ['cortex_left','cortex_right']
        mask = []
        for i,hem in enumerate(['L','R']):
            mask.append(atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-{hem}_mask.label.gii')
        atlas = AtlasSurface('fs32k', mask_gii=mask, structure=bm_name)
    elif atlas_str == 'fs32k_L':
        bm_name = ['cortex_left']
        mask = [atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii']
        atlas = AtlasSurface('fs32k', mask_gii=mask, structure=bm_name)
    elif atlas_str == 'fs32k_R':
        bm_name = ['cortex_right']
        mask = [atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii']
        atlas = AtlasSurface('fs32k', mask_gii=mask, structure=bm_name)
    else:
        raise(NameError(f'Unknown atlas string:{atlas_str}'))
    return atlas


class Atlas():
    """The Atlas class implements the general atlas functions
    for mapping from the P brain locations back to nii or gifti files
    Each Atlas is associated with a set of atlas maps
    """
    def __init__(self,name):
        self.name = name
        self.P = np.nan # Number of locations in this atlas

    def map_data(self,data):
        """Maps data back into some atlas form.
        Args:
            data (numpy.ndarray): P or N x P array
        """

class AtlasVolumetric(Atlas):
    """ Volumetric atlas with specific 3d-locations
    """
    def __init__(self,name,mask_img):
        """Atlas Volumetric class constructor

        Args:
            name (str): Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_img (str): file name of mask image defining atlas location
        """
        super().__init__(name)
        self.mask_img = nb.load(mask_img)
        Xmask = self.mask_img.get_data()
        Xmask = (Xmask>0)
        i,j,k = np.where(Xmask>0)
        self.vox = np.vstack((i,j,k))
        self.world = nt.affine_transform_mat(self.vox,self.mask_img.affine)
        self.P = self.world.shape[1]

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

    def get_brain_model_axis(self):
        """ Returns brain model axis

        Returns:
            bm (cifti2.BrainModelAxis)
        """
        bm = nb.cifti2.BrainModelAxis.from_mask(self.mask_img.get_data(),
                                            name=self.name,
                                            affine = self.mask_img.affine)
        return bm

    def data_to_nifti(self,data):
        """Transforms data in atlas space into
        3d or 4d nifti image

        Args:
            data (np.ndarray): Data to be mapped into nift

        Returns:
            Nifti1Image: NiftiImage object
        """
        if data.ndim==1:
            data = data.reshape(1,-1)
        N,p = data.shape
        if p != self.P:
            raise(NameError('Data needs to be a P vector or NxP matrix'))
        if N>1:
            X=np.zeros(self.mask_img.shape+(N,))
            X[self.vox[0],self.vox[1],self.vox[2]]=data.T
        else:
            X=np.zeros(self.mask_img.shape)
            X[self.vox[0],self.vox[1],self.vox[2]]=data
        img = nb.Nifti1Image(X,self.mask_img.affine)
        return img

    def sample_nifti(self,img,interpolation):
        """ Samples a img at the atlas locations
        The image needs to be in atlas space.

        Args:
            img (str or NiftiImage): Nifti to be sampled

        Returns:
            np.array: Data sample at the atlas position
        """
        if isinstance(img,str):
            img = nb.load(img)
        data = nt.sample_image(img,
                            self.world[0],
                            self.world[1],
                            self.world[2],
                            interpolation)
        return data

class AtlasVolumeSymmetric(AtlasVolumetric):
    """ Volumetric atlas with left-right symmetry
    The atlas behaves like AtlasVolumetrc, but provides
    mapping indices from a full representation to
    a reduced (symmetric) representation of size Psym.
    """
    def __init__(self,name,mask_img):
        """AtlasVolumeSymmeytric class constructor: Generates members
        indx_full, indx_reduced, indx_flip.
        Assume you have a
            Full: N x P array
            Left: N x Psym array
            Right: N x Psym array
        then:
            Left = Full[:,index_full[0]]
            Right = Full[:,index_full[1]]
            Avrg = (Left + Right)/2
            Full = Avrg[:,index_reduced]
        To Flip:
            flippedFull = Full[:,index_flip]
        Args:
            name (str): Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_img (str): file name of mask image defining atlas location
        """
        super().__init__(name,mask_img)
        # Find left and righgt hemispheric voxels
        indx_left = np.where(self.world[0,:]<=0)[0]
        indx_right = np.where(self.world[0,:]>=0)[0]
        # Make a version with absolute x-coordiate
        world_coord = self.world.copy()
        world_coord[0,:]=np.abs(world_coord[0,:])
        # Initialize indices
        self.Psym = indx_left.shape[0]
        self.indx_full = np.zeros((2,self.Psym),dtype=int)
        self.indx_full[0,:] = indx_left
        self.indx_reduced = np.zeros((self.P,),dtype=int)
        self.indx_reduced[indx_left] = np.arange((self.Psym))

        # Now find a ordering of the right hemisphere
        # that matches the left hemisphere
        for i in range(self.Psym):
            a=np.nonzero(np.all(world_coord[:,i:i+1]==self.world[:,indx_right],axis=0))[0]
            if len(a)!=1:
                raise(NameError('The voxels in mask do not seem to be fully symmetric along the x-axis'))
            self.indx_full[1,i]=indx_right[a[0]]
            self.indx_reduced[indx_right[a[0]]]=i
        # Generate flipping index
        indx_orig = np.arange(self.P,dtype=int)
        self.indx_flip = np.zeros((self.P,),dtype=int)
        self.indx_flip[self.indx_full[0]]=indx_orig[self.indx_full[1]]
        self.indx_flip[self.indx_full[1]]=indx_orig[self.indx_full[0]]

class AtlasSurface(Atlas):
    """Surface-based altas space
    """
    def __init__(self, name, mask_gii, structure):
        """Atlas Surface class constructor
        Args:
            name (str): Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_gii (list): gifti file name of mask image defining atlas locations
            structure (list): [cortex_left, cortex_right] gifti file name of mask image defining
            atlas locations
        Notes:
            Since this class is called 'AtlasSurface', I think we should
            only integrate surface datas in cifti brain structures so that
            the volume data shouldn't be in.   -- dzhi
        """
        super().__init__(name)

        assert len(mask_gii) == len(structure), \
            "The length of mask and brain structure should be matched!"

        self.structure = structure
        Xmask = [nb.load(mg).agg_data() for mg in mask_gii]
        self.vertex_mask = [(X>0) for X in Xmask]
        self.vertex = [np.nonzero(X)[0] for X in self.vertex_mask]
        self.P = sum([v.shape[0] for v in self.vertex])

    def data_to_cifti(self, data, row_axis=None):
        """Maps data back into a cifti image
        Args:
            data: the input data to be mapped
                (ndarray) - 1-d Numpy array of the size (P,)
                (list) - list of ndarray
            row_axis: label for row axis in cifti file, it can be
                (list) - a list of colum names
                (object) - a pandas framework object of the colum names
                (cifti2.Axis) - a cifti2 Axis object that can be directly
                                used to write header. (e.g. ScalarAxis,
                                SeriesAxis, ...)
                None - default to generate a list of column names that
                       matches the input data
        Returns:
            Cifti2Image: Cifti2Image object
        """
        if isinstance(data, list):
            # Check #1: the passing data is a list, then it should match
            # the number of brain structures in current atlasMap, and the
            # number of vertices in each structure should be aligned.
            assert len(data) == len(self.structure), \
                "The length of data and brain structure should be matched!"
            for i, dat in enumerate(data):
                assert dat.shape[1] == self.vertex[i].shape[0], \
                    f"The number of vertices doesn't match in {self.structure[i]}"

            # If list is passed, horizontal stack them.
            data = np.hstack(data)
        elif isinstance(data, np.ndarray):
            # Check #2: the passing data is np.ndarray, then check the
            # number of vertices whether is matched with data
            assert data.shape[1] == self.P, \
                "The length of data and brain structure should be matched!"
        else:
            raise ValueError('The input data must be a np.ndarray or a list of np.ndarray')

        if row_axis is None:
            row_axis = [f'row {r:03}' for r in range(data.shape[0])]
            row_axis = nb.cifti2.ScalarAxis(row_axis)
        elif hasattr(row_axis, '__iter__'):
            assert data.shape[0] == len(row_axis), \
                "The length of row_axis should match the data!"
            row_axis = nb.cifti2.ScalarAxis(row_axis)
        elif isinstance(row_axis, nb.cifti2.cifti2_axes.Axis):
            pass
        else:
            raise ValueError('The input row_axis instance type does not meet the requirement!')

        bm = self.get_brain_model_axis()
        header = nb.Cifti2Header.from_axes((row_axis, bm))
        cifti_img = nb.Cifti2Image(dataobj=data, header=header)

        return cifti_img

    def cifti_to_data(self, cifti):
        """Gets the data from a CIFTI file, checking the structure name
           and vertex information in the cifti file. If it doesn't match
           the vertex information in the atlas object, it gives a warning,
           but corrects for it by extracting the available data.

        Args:
            cifti (ciftiimage or filename): Cifti file to be used
        Returns: 
            np.ndarray: NxP in single np-array
        """
        # First check the input is a cifti image object
        if isinstance(cifti, str):
            cifti = nb.load(cifti)

        data = cifti.get_fdata()
        # Second check if the data has the same vertices index or
        # the number of vertices with the current atlasMap
        col_axis = cifti.get_header().get_axis(1)
        if np.array_equal(np.hstack(self.vertex), col_axis.vertex):
            pass
        else:
            warnings.warn('The input cifti image does not match the current '
                          'atlas_map. The atlas map attributes will be changed to'
                          'align the input Cifti2Image.')
            names, idx = np.unique(col_axis.name, return_index=True)
            self.structure = names
            self.vertex = np.split(col_axis.vertex, idx[1:])

            # Align self.vertex_mask to the cifti image
            self.vertex_mask = []
            for i, stru_nam in enumerate(self.structure):
                mask = np.full((col_axis.nvertices[stru_nam],), False, dtype=bool)
                mask[self.vertex[i]] = True
                self.vertex_mask.append(mask)

            self.P = col_axis.size

        return data

    def get_brain_model_axis(self):
        """ Returns brain model axis

        Returns:
            bm (cifti2.BrainModelAxis)
        """
        # Make the brain Structure models
        for i, name in enumerate(self.structure):
            if i == 0:
                bm = nb.cifti2.BrainModelAxis.from_mask(
                    self.vertex_mask[i],
                    name=self.structure[i])
            else:
                bm = bm + nb.cifti2.BrainModelAxis.from_mask(
                    self.vertex_mask[i],
                    name=self.structure[i])
        return bm

class AtlasVolumeParcel(Atlas):
    """ Volume-based atlas that is based on
    """
    def __init__(self,name,label_img,mask_img=None):
        """AtlasSurfaceParcel class constructor

        Args:
            name (str):
                Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_gii (str):
                gifti file name of mask image with the allowed locations - note that the label_vec is still indexing the the original surface space
        """
        self.name = name
        self.label_img = nb.load(label_img)

        if mask_img is not None:
            self.mask_img = nb.load(mask_img)
        else:
            self.mask_img = self.label_img

        Xmask = self.mask_img.get_data()
        i,j,k = np.where(Xmask>0)
        self.vox = np.vstack((i,j,k))
        self.world = nt.affine_transform_mat(self.vox,self.mask_img.affine)
        self.label_vec = suit.reslice.sample_image(self.label_img,
                    self.world[0],
                    self.world[1],
                    self.world[2],
                    0)
        # Find the number of parcels with label > 0
        self.P = np.unique(self.label_vec[self.label_vec>0]).shape[0]

    def agg_data(self,data,func=np.nanmean):
        """Aggregate data into the parcels -
        Note that by convention the lable=0 is reserved for vertices that are being ignored - the first column is label_vec=1....

        Args:
            data (nparray):
                NxP array of data
            func (function):
                Aggregation function - needs to take data,axis=1 as input arguments
        Returns:
            Aggregated data:
        """

        # get the aggregrated data - assumes consequtive labels [1...P]
        agg_data = np.empty((data.shape[0],self.P))
        for p in range(self.P):
            agg_data[:,p]=func(data[:,self.label_vec==p+1],axis=1)
        return agg_data


    def get_parcel_axis(self):
        """ Returns parcel axis

        Returns:
            bm (cifti2.ParcelAxis)
        """

        # loop over labels and create brain models
        bm_list = []
        for l in np.unique(self.label_vec):
            if l > 0:
                # Create BrainModelAxis for each label
                bm = nb.cifti2.BrainModelAxis.from_mask(self.label_vec == l,
                                                    name=self.name)
                # append a tuple containing the name of the parcel and the corresponding BrainAxisModel
                bm_list.append((f"{self.name}-{l:02}", bm))

        # create parcel axis from the list of brain models created for labels
        self.parcels = nb.cifti2.ParcelsAxis.from_brain_models(bm_list)

        return self.parcels

class AtlasSurfaceParcel(Atlas):
    def __init__(self,name,label_gii,mask_gii=None):
        """AtlasSurfaceParcel class constructor

        Args:
            name (str):
                Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_gii (str):
                gifti file name of mask image with the allowed locations - note that the label_vec is still indexing the the original surface space
        """
        self.name = name
        self.label_gii = nb.load(label_gii)

        # get the labels into an array
        self.label_vec = self.label_gii.agg_data() # this

        # Use of mask it not tested-treat with care
        if mask_gii is not None:
            self.mask_gii = nb.load(mask_gii)
            Xmask = self.mask_gii.agg_data()
            self.label_vec = np.delete(self.label_vec, Xmask==0)
        # Find the number of parcels with label > 0
        self.P = np.unique(self.label_vec[self.label_vec>0]).shape[0]

    def agg_data(self,data,func=np.nanmean):
        """Aggregate data into the parcels -
        Note that by convention the lable=0 is reserved for vertices that are being ignored - the first column is label_vec=1....

        Args:
            data (nparray):
                NxP array of data
            func (function):
                Aggregation function - needs to take data,axis=1 as input arguments
        Returns:
            Aggregated data:
        """

        # get the aggregrated data - assumes consequtive labels [1..P]
        agg_data = np.empty((data.shape[0],self.P))
        for p in range(self.P):
            agg_data[:,p]=func(data[:,self.label_vec==p+1],axis=1)
        return agg_data


    def get_parcel_axis(self):
        """ Returns parcel axis

        Returns:
            bm (cifti2.ParcelAxis)
        """

        # loop over labels and create brain models
        bm_list = []
        for l in np.unique(self.label_vec):
            if l > 0:
                # Create BrainModelAxis for each label
                bm = nb.cifti2.BrainModelAxis.from_mask(self.label_vec == l,
                                                    name=self.name)
                # append a tuple containing the name of the parcel and the corresponding BrainAxisModel
                bm_list.append((f"{self.name}-{l:02}", bm))

        # create parcel axis from the list of brain models created for labels
        self.parcels = nb.cifti2.ParcelsAxis.from_brain_models(bm_list)

        return self.parcels

class AtlasMap():
    def __init__(self, dataset, name, P , participant_id):
        """AtlasMap stores the mapping rules from a specific data set (and participant) to the desired atlas space in form of a voxel list
        Args:
            dataset_id (string): name of
            participant_id (string): Participant name
        """
        self.dataset = dataset # Reference to corresponding data set
        self.name = name
        self.P = P       #  Number of brain locations
        self.participant_id = participant_id

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
    def __init__(self, dataset, name, world, participant_id, deform_img,mask_img):
        """AtlasMapDeform stores the mapping rules for a non-linear deformation
        to the desired atlas space in form of a voxel list
        Args:
            dataset_id (str): name of
            name (str): Name of atlas map
            worlds (ndarray): 3xP ND array of world locations
            participant_id (str): Participant name
            deform_img (str/list): Name of deformation map image(s)
            mask_img (str): Name of masking image that defines the functional data space.
        """
        P = world.shape[1]
        super().__init__(dataset,name,P,participant_id)
        self.world = world
        if type(deform_img) is not list:
            deform_img = [deform_img]
        self.deform_img = []
        for di in deform_img:
            self.deform_img.append(nb.load(di))
        self.mask_img = nb.load(mask_img)

    def build(self, smooth = None, additional_mask=None):
        """
        Using the dataset, builds a list of voxel indices of
        For each of the locations. It creates:
        vox_list: List of voxels to sample for each atlas location
        vox_weight: Weight of each of these voxels to determine the atlas location
        Arg:
            smooth (double): SD of smoothing kernel (mm) or None for nearest neighbor
            additional_mask: Additional Mask image (not necessarily in functional space - only voxels with elements > 0 in that image
            will be used for the altas )
        """
        # Caluculate locations of atlas in individual (deformed) coordinates
        # Apply possible multiple deformation maps sequentially
        xyz = self.world.copy()
        for i,di in enumerate(self.deform_img):
            xyz = nt.sample_image(di,
                    xyz[0],
                    xyz[1],
                    xyz[2],1).squeeze().T
            pass
        atlas_ind = xyz
        N = atlas_ind.shape[1] # Number of locations in atlas

        # Determine which voxels are available in functional space
        # and apply additional mask if given
        M = self.mask_img.get_fdata()
        i,j,k=np.where(M>0)
        vox = np.vstack((i,j,k))
        world_vox = nt.affine_transform_mat(vox,self.mask_img.affine) # available voxels in world coordiantes
        if additional_mask is not None:
            # If file name, load the nifti image
            if isinstance(additional_mask,str):
                additional_mask = nb.load(additional_mask)
            add_mask = nt.sample_image(additional_mask,
                        world_vox[0],
                        world_vox[1],
                        world_vox[2],1)
            world_vox = world_vox[:,add_mask>0]
            vox = vox[:,add_mask>0]

        if smooth is None: # Use nearest neighbor interpolation
            linindx,good = nt.coords_to_linvidxs(atlas_ind,self.mask_img,mask=True)
            self.vox_list = linindx.reshape(-1,1)
            self.vox_weight = np.ones((linindx.shape[0],1))
            self.vox_weight[np.logical_not(good)]=np.nan
        else:              # Use smoothing kernel of specific size
            linindx = np.ravel_multi_index((vox[0,:],vox[1,:],vox[2,:]),
                                            M.shape,mode='clip')
            # Distances between atlas coordinates and voxel coordinates
            D = nt.euclidean_dist_sq(atlas_ind,world_vox)
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
            # Avoid divide by zero error:
            mw = self.vox_weight.sum(axis=1, keepdims=True)
            mw[mw==0]=np.nan
            self.vox_weight = self.vox_weight / mw
        pass

class AtlasMapSurf(AtlasMap):
    def __init__(self, dataset, name, vertex, participant_id,
                white_surf,pial_surf,mask_img):
        """AtlasMapSurf stores the mapping rules for a freesurfer-style surface (pial + white surface pair)
        Args:
            dataset_id (str): name of
            participant_id (str): Participant name
            white_surf (str): Name for white matter surface
            pial_surf (str): Name for pial surface
            mask_img (str): Name of masking image that defines the functional data space.
        """
        P = len(vertex)
        super().__init__(dataset,name,P,participant_id)
        self.vertex = vertex
        self.white_surf = nb.load(white_surf)
        self.pial_surf = nb.load(pial_surf)
        self.mask_img = nb.load(mask_img)

    def build(self,smooth = None, depths=[0,0.2,0.4,0.6,0.8,1.0]):
        """
        Using the dataset, builds a list of voxel indices of
        each of the nodes
        vox_list: List of voxels to sample for each atlas location
        vox_weight: Weight of each of these voxels to determine the atlas location
        Arg:
            depths (iterable): List of depth between pial (1) and white (0) surface that
            will be sampled
        """
        n_points = len(depths)
        c1 = self.white_surf.darrays[0].data[self.vertex,:].T
        c2 = self.pial_surf.darrays[0].data[self.vertex,:].T
        n_vert = c1.shape[1]
        if c2.shape[1] != n_vert:
            raise(NameError('White and pial surfaces should have same number of vertices.'))

        # Get the indices for all the points being sampled
        indices = np.zeros((n_points,3,n_vert))
        for i in range(n_points):
            indices[i,:,:] = (1-depths[i])*c1+depths[i]*c2

        self.vox_list,good = nt.coords_to_linvidxs(indices,self.mask_img,mask=True)
        all = good.sum(axis=0)
        print(f'{self.name} has {np.sum(all==0)} vertices without data')
        all[all==0]=1
        self.vox_weight = good / all
        self.vox_list = self.vox_list.T
        self.vox_weight = self.vox_weight.T

def get_data3D(fnames,atlas_maps):
    """Extracts the data for a list of fnames
    for a list of atlas_maps. This is usually called by DataSet.get_data()
    to extract the required raw data before processing it further

    Args:
        fnames (list): list of N file names to be sampled
        atlas_maps (list): list of K built atlas-map objects
    returns:
        data (list): List of NxP_k 2-d array data matrices (np)
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
            raise(NameError('use get_data4D for 4D data'))
        # Map this file into the data structures
        X = X.ravel()
        for i,at in enumerate(atlas_maps):
            d=X[at.vox_list] * at.vox_weight  # Expanded data
            d = np.nansum(d,axis=1)
            d[np.nansum(at.vox_weight,axis=1)==0]=np.nan
            data[i][j,:]=d
    return data

def get_data4D(vol_4D,atlas_maps):
    """Extracts the data for a list of fnames
    for a list of atlas_maps. This is usually called by DataSet.get_data()
    to extract the required raw data before processing it further

    Args:
        fnames (list): list of file names to be sampled
        atlas_maps (list): list of built atlas-map objects
    """

    # temp code:
    # convert 4D to 3D
    vol_3D_list = nb.funcs.four_to_three(vol_4D)

    n_atlas = len(atlas_maps)
    n_files = len(vol_3D_list)
    data = []
    # Make the empty data structures
    for at in atlas_maps:
        data.append(np.full((n_files,at.vox_list.shape[0]),np.nan))
    for j in range(len(vol_3D_list)):
        V = vol_3D_list[j]
        X = V.get_fdata()
        # Map this file into the data structures
        X = X.ravel()
        for i,at in enumerate(atlas_maps):
            d=X[at.vox_list] * at.vox_weight  # Expanded data
            d = np.nansum(d,axis=1)
            d[np.nansum(at.vox_weight,axis=1)==0]=np.nan
            data[i][j,:]=d
    return data


def data_to_cifti(data,atlas_maps,names=None):
    """Transforms a list of data sets and list of atlas maps
    into a cifti2image

    Args:
        data (list):
            List / array of data arrays - need to have all same shape[0]
            and a shape[1] that matches the corresponding atlas map
        atlas_maps (list):
            List / array of atlas maps
        names (list of str):
            Names for the scalar axis
    Returns:
        img: nibabel.cifti2image
            Can be saved as (*.dscalar.nii) file
    """
    # Check is a single is given
    if type(data) is not list:
        data = [data]
    if type(atlas_maps) is not list:
        atlas_maps = [atlas_maps]

    # Make the brain Structure models
    for i,atm in enumerate(atlas_maps):
        if i == 0:
            bm = atm.atlas.get_brain_model_axis()
            D = data[i]
        else:
            bm = bm+atm.atlas.get_brain_model_axis()
            D = np.c_[D,data[i]]

    # row_axis = nb.cifti2.SeriesAxis(start=0,step=1,size=D.shape[0])
    if names is None:
        names = [f'row {r:02}' for r in range(D.shape[0])]
    row_axis = nb.cifti2.ScalarAxis(names)
    header = nb.Cifti2Header.from_axes((row_axis,bm))
    cifti_img = nb.Cifti2Image(dataobj=D,header=header)
    return cifti_img


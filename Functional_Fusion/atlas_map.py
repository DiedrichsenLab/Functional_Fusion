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
import warnings
import nitools as nt
import json
import re
import os

# Need to do self import here to get Atlas directory
import Functional_Fusion.atlas_map as am
default_atlas_dir = os.path.dirname(am.__file__) + '/Atlases'

def get_atlas(atlas_str, atlas_dir=default_atlas_dir):
    """returns an atlas from a code

    Args:
        atlas_str (str): Name of the atlas
        atlas_dir (str): Atlas Directory (defaults to Functional_Fusion/Atlases)

    Returns:
        atlas (AtlasVolumetric, AtlasSurface, (or symmetric variant) object): Atlas object
        ainf (dict): Atlas info
    """
    with open(atlas_dir + "/atlas_description.json") as file:
        atlases = json.load(file)
    if atlas_str not in atlases:
        raise (NameError(f"Unknown Atlas: {atlas_str}"))
    ainf = atlases[atlas_str]

    # Make the atlas object
    At = getattr(am, ainf["type"])
    if isinstance(ainf["mask"], list):
        mask = []
        for i, m in enumerate(ainf["mask"]):
            mask.append(atlas_dir + f"/{ainf['dir']}/{m}")
    else:
        mask = f"{atlas_dir}/{ainf['dir']}/{ainf['mask']}"
    atlas = At(atlas_str, mask)
    atlas.structure = ainf["structure"]
    atlas.space = ainf["space"]
    return atlas, ainf

def get_deform(target, source="MNIAsym2",atlas_dir = default_atlas_dir):
    """Get name of group deformation map between two volumetric atlas spaces

    Args:
        target (str/atlas): Target space
        source (str): Source space
        atlas_dir (str): Atlas Directory (defaults to Functional_Fusion/Atlases)
    Returns:
        deform (str): Name of deformation map
        mask (str): Name of mask for the source space
    """
    if isinstance(target, str):
        target, _ = get_atlas(target, atlas_dir)
    with open(atlas_dir + "/atlas_description.json") as file:
        atlases = json.load(file)
    if target.name not in atlases:
        raise (NameError(f"Unknown Atlas: {target.name}"))
    if source not in atlases:
        raise (NameError(f"Unknown Atlas: {source}"))
    tar = atlases[target.name]
    src = atlases[source]
    deform = f"{atlas_dir}/{src['dir']}/tpl-{src['space']}_space-{tar['space']}_xfm.nii"
    mask = f"{atlas_dir}/{src['dir']}/{src['mask']}"
    return deform, mask

def parcel_recombine(label_vector,parcels_selected,label_id=None,label_name=None):
    """ Recombines and selects different parcels into a new label vector for ROI analysis.
    IMPORTANT: Note that each old parcel can only be mapped to one new parcel. That is the elements
    of parcels_selected must be non-overlapping.

    Args:
        label_vector (ndarrray): voxel array of labels
        parcels_selected (list):
            list of the parcels being selected. Each element should be
            a) Scalar value of lable_ids,
            b) List/array of scalar values of label_ids
            c) A regexp - for example 'D.L' would select all labels starting with D.L
        label_id (ndarrray): ndarray of parcel ids
        label_name (list): List of parcel names

    Returns:
        label_vector: new coded vector, with the first ROI label with 1. all others 0
        label_id: new label id
        label_name: new label name
    """

    # If all parcels are selected, just return the original
    if (parcels_selected is None) or (parcels_selected == "all"):
        return label_vector, label_id, label_name
    elif isinstance(parcels_selected, (list,np.ndarray)):
        label_vector_new = np.zeros(label_vector.shape,dtype=int)
        label_id_new = np.arange(len(parcels_selected)+1)
        label_name_new = ['0']
        for i,p in enumerate(parcels_selected):
            if isinstance(p,int):
                indx=[p]
                label_name_new.append(label_name[np.nonzero(label_id==p)[0][0]])
            elif isinstance(p,str):
                indx = [label_id[i] for i, x in enumerate(label_name) if re.search(p, x) is not None]
                label_name_new.append(p)
            elif isinstance(p,list):
                indx = p
                label_name_new.append(f'reg_{i}')
            label_vector_new[np.isin(label_vector,indx)]=i+1
    else:
        raise ValueError('parcels_selected must be a list')
    return label_vector_new, label_id_new, label_name_new

class Atlas:
    def __init__(self, name):
        """ The Atlas class implements the mapping from the P brain locations back to the defining
        Volumetric or surface space. It also allows for parcellation of these brain locations.
        Each Atlas is associated with a set of atlas maps

        Args:
            name (str): Unique code for the atlas (e.g. 'SUIT3')
        """
        self.name = name
        self.P = np.nan  # Number of locations in this atlas
        self.structure = 'unknown'
        self.space = 'unknown'

    def get_parcel(self, label_img):
        """adds a label_vec to the atlas that assigns each voxel / node of the atlas
        label_vec[i] = 0 mean that the ith voxel / node is unassigned

        Args:
            label_img (str/list) - filename(s) of label image(s)
        Returns:
            label_vec (list) - list of numpy array containing label values corresponding to label images (0 means no ROI)
            labels (ndarray) - List of unique labels > 0
        """

        # get label vectors for each label image
        self.label_vector = self.read_data(label_img, 0).astype(int)
        self.labels = np.unique(self.label_vector[self.label_vector > 0])
        self.n_labels = self.labels.shape[0]
        return self.label_vector, self.labels


class AtlasVolumetric(Atlas):
    def __init__(self, name, mask_img):
        """Atlas Volumetric is an atlas for a Volumetric
        space / structure.

        Args:
            name (str): Name of atlas (atlas string)
            mask_img (str): file name of mask image defining atlas location
            structure (str): the brain structure name for Cifti (thalamus, cerebellum)
        """
        super().__init__(name)
        self.mask_img = nb.load(mask_img)
        Xmask = self.mask_img.get_fdata()
        Xmask = Xmask > 0
        i, j, k = np.where(Xmask > 0)
        self.vox = np.vstack((i, j, k))
        self.mask = Xmask
        self.world = nt.affine_transform_mat(self.vox, self.mask_img.affine)
        self.P = self.world.shape[1]

    def get_brain_model_axis(self):
        """Returns brain model axis

        Returns:
            bm (cifti2.BrainModelAxis)
        """
        bm = nb.cifti2.BrainModelAxis.from_mask(
            self.mask_img.get_fdata(), name=self.structure, affine=self.mask_img.affine
        )
        return bm

    def get_parcel_axis(self):
        """Returns parcel axis based on the current setting of labels

        Returns:
            bm (cifti2.ParcelAxis)
        """

        bm_list = []
        if not hasattr(self, "labels"):
            raise (NameError("Atlas has no parcels defined yet - call atlas.getparcels(image) first")
            )
        bm_list = []
        for h, label in enumerate(self.labels):
            pname = self.structure + f"_{label:02d}"
            indx = self.label_vector == label
            bm = nb.cifti2.BrainModelAxis(
                self.structure,
                voxel=self.vox[:, indx].T,
                affine=self.mask_img.affine,
                volume_shape=self.mask_img.shape,
            )
            bm_list.append((pname, bm))
        parcel_axis = nb.cifti2.ParcelsAxis.from_brain_models(bm_list)
        return parcel_axis

    def data_to_cifti(self, data, row_axis=None):
        """Transforms data into a cifti image

        Args:
            data (ndarray/list): the input data to be mapped 1-d Numpy array of the size (P,)
            row_axis: label for row axis in cifti file, it can be:

                | (list) - a list of row names
                | (object) - a pandas framework object of the row names
                | (cifti2.Axis) - a cifti2 Axis object that can be directly used to write header. (e.g. ScalarAxis, SeriesAxis, ...)
                | (None) - default to generate a list of row names that matches the input data

        Returns:
            Cifti2Image: Cifti2Image object
        """
        if isinstance(data, list):
            data = data[0]
        assert (
            data.shape[1] == self.P
        ), "The length of data and brain structure should be matched!"
        if row_axis is None:
            row_axis = [f"row {r:03}" for r in range(data.shape[0])]
            row_axis = nb.cifti2.ScalarAxis(row_axis)
        elif hasattr(row_axis, "__iter__"):
            assert data.shape[0] == len(
                row_axis
            ), "The length of row_axis should match the data!"
            row_axis = nb.cifti2.ScalarAxis(row_axis)
        elif isinstance(row_axis, nb.cifti2.cifti2_axes.Axis):
            pass
        else:
            raise ValueError(
                "The input row_axis instance type does not meet the requirement!"
            )

        bm = self.get_brain_model_axis()
        header = nb.Cifti2Header.from_axes((row_axis, bm))
        cifti_img = nb.Cifti2Image(dataobj=data, header=header)
        return cifti_img

    def data_to_nifti(self, data):
        """Transforms data in atlas space into 3d or 4d nifti image depending on data type, the empty parts of the image will be NaNs or zeros.

        Args:
            data (np.ndarray): Data to be mapped into nifti
        Returns:
            Nifti1Image(nb.Nifti1Image): NiftiImage object
        """
        if data.ndim == 1:
            data = data.reshape(1, -1)
        N, p = data.shape
        if p != self.P:
            raise (NameError("Data needs to be a P vector or NxP matrix"))
        if N > 1:
            if np.issubdtype(data.dtype, np.floating):
                X = np.empty(self.mask_img.shape + (N,), dtype=float)
                X.fill(np.nan)
            elif np.issubdtype(data.dtype, np.integer):
                X = np.zeros(self.mask_img.shape + (N,), dtype=int)
            else:
                raise ValueError("Data type not supported")
            X[self.vox[0], self.vox[1], self.vox[2]] = data.T
        else:
            if np.issubdtype(data.dtype, np.floating):
                X = np.empty(self.mask_img.shape, dtype=float)
                X.fill(np.nan)
            elif np.issubdtype(data.dtype, np.integer):
                X = np.zeros(self.mask_img.shape, dtype=int).astype(np.int8)
            else:
                raise ValueError("Data type not supported")
        
            X[self.vox[0], self.vox[1], self.vox[2]] = data
        img = nb.Nifti1Image(X, self.mask_img.affine)
        return img

    def read_data(self, img, interpolation=0):
        """Read data from a NIFTI or CIFTI file into the volumetric atlas space

        Args:
            img (nibabel.image): Nifti or Cifti image
            interpolation (int)): nearest neighbour (0), trilinear (1)
        Returns:
            data (ndarray): (N,P) ndarray
        """
        if isinstance(img, str):
            img = nb.load(img)
        if isinstance(img, nb.Cifti2Image):
            img = nt.volume_from_cifti(img, [self.structure])
        if isinstance(img, (nb.Nifti1Image, nb.Nifti2Image)):
            data = nt.sample_image(
                img, self.world[0, :], self.world[1, :], self.world[2, :], interpolation
            )
        else:
            raise(NameError("Unknown image type"))
        return data

    def sample_nifti(self, img, interpolation):
        """Samples a img at the atlas locations
        The image needs to be in atlas space.

        Args:
            img (str or NiftiImage): Nifti to be sampled

        Returns:
            np.array: Data sample at the atlas position
        """
        warnings.DeprecationWarning(
            "sample_nifti is depreciated. Use self.read_data instead"
        )
        if isinstance(img, str):
            img = nb.load(img)
        data = nt.sample_image(
            img, self.world[0], self.world[1], self.world[2], interpolation
        )
        return data


class AtlasVolumeSymmetric(AtlasVolumetric):
    def __init__(self, name, mask_img, structure="cerebellum"):
        """Volumetric atlas with left-right symmetry. The atlas behaves like AtlasVolumetric,
        but provides mapping indices from a full representation to a
        reduced (symmetric) representation of size Psym.
        The constructor generates members indx_full, indx_reduced, indx_flip.

        | Assume you have a
        | Full: N x P array
        | Left: N x Psym array
        | Right: N x Psym array
        | then:
        | Left = Full[:,index_full[0]]
        | Right = Full[:,index_full[1]]
        | Avrg = (Left + Right)/2
        | Full = Avrg[:,index_reduced]
        | To Flip:
        | flippedFull = Full[:,index_flip]

        Args:
            name (str): Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_img (str): file name of mask image defining atlas location
        """

        super().__init__(name, mask_img, structure)
        # Find left and right hemispheric voxels
        indx_left = np.where(self.world[0, :] <= 0)[0]
        indx_right = np.where(self.world[0, :] >= 0)[0]
        # Make a version with absolute x-coordiate
        world_coord = self.world.copy()
        world_coord[0, :] = np.abs(world_coord[0, :])
        # Initialize indices
        self.Psym = indx_left.shape[0]
        self.indx_full = np.zeros((2, self.Psym), dtype=int)
        self.indx_full[0, :] = indx_left
        self.indx_reduced = np.zeros((self.P,), dtype=int)
        self.indx_reduced[indx_left] = np.arange((self.Psym))

        # Now find a ordering of the right hemisphere
        # that matches the left hemisphere
        for i in range(self.Psym):
            a = np.nonzero(
                np.all(world_coord[:, i : i + 1] == self.world[:, indx_right], axis=0)
            )[0]
            if len(a) != 1:
                raise (
                    NameError(
                        "The voxels in mask do not seem to be fully symmetric along the x-axis"
                    )
                )
            self.indx_full[1, i] = indx_right[a[0]]
            self.indx_reduced[indx_right[a[0]]] = i
        # Generate flipping index
        indx_orig = np.arange(self.P, dtype=int)
        self.indx_flip = np.zeros((self.P,), dtype=int)
        self.indx_flip[self.indx_full[0]] = indx_orig[self.indx_full[1]]
        self.indx_flip[self.indx_full[1]] = indx_orig[self.indx_full[0]]


class AtlasSurface(Atlas):
    """Surface-based altas space"""

    def __init__(self, name, mask_img, structure):
        """Atlas Surface class constructor
        Args:
            name (str): Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_gii (list): gifti file name of mask image defining atlas locations
            structure (list): [cortex_left, cortex_right] - CIFTI brain structure names
        Notes:
            Since this class is called 'AtlasSurface', I think we should
            only integrate surface datas in cifti brain structures so that
            the volume data shouldn't be in.   -- dzhi
        """
        super().__init__(name)

        assert len(mask_img) == len(
            structure
        ), "The length of mask and brain structure should be matched!"

        self.structure = structure
        Xmask = [nb.load(mg).agg_data() for mg in mask_img]
        self.mask = [(X > 0) for X in Xmask]
        self.vertex_mask = [(X > 0) for X in Xmask]
        self.vertex = [np.nonzero(X)[0] for X in self.vertex_mask]
        # Correction: needs to be number of vertices in mask (JD)
        self.P = sum([v.shape[0] for v in self.vertex])

    def data_to_cifti(self, data, row_axis=None):
        """Maps data back into a cifti image

        Args:
            data: the input data to be mapped

                | (ndarray) - 1-d Numpy array of the size (P,)
                | (list) - list of ndarray

            row_axis: label for row axis in cifti file, it can be

                | (list) - a list of colum names
                | (object) - a pandas framework object of the colum names
                | (cifti2.Axis) - a cifti2 Axis object that can be directly used to write header. (e.g. ScalarAxis, SeriesAxis, ...)
                | None - default to generate a list of column names that matches the input data

        Returns:
            Cifti2Image: Cifti2Image object
        """
        if isinstance(data, list):
            # Check #1: the passing data is a list, then it should match
            # the number of brain structures in current atlasMap, and the
            # number of vertices in each structure should be aligned.
            assert len(data) == len(
                self.structure
            ), "The length of data and brain structure should be matched!"
            for i, dat in enumerate(data):
                assert (
                    dat.shape[1] == self.vertex[i].shape[0]
                ), f"The number of vertices doesn't match in {self.structure[i]}"

            # If list is passed, horizontal stack them.
            data = np.hstack(data)
        elif isinstance(data, np.ndarray):
            # Check #2: the passing data is np.ndarray, then check the
            # number of vertices whether is matched with data
            assert (
                data.shape[1] == self.P
            ), "The length of data and brain structure should be matched!"
        else:
            raise ValueError(
                "The input data must be a np.ndarray or a list of np.ndarray"
            )

        if row_axis is None:
            row_axis = [f"row {r:03}" for r in range(data.shape[0])]
            row_axis = nb.cifti2.ScalarAxis(row_axis)
        elif hasattr(row_axis, "__iter__"):
            assert data.shape[0] == len(
                row_axis
            ), "The length of row_axis should match the data!"
            row_axis = nb.cifti2.ScalarAxis(row_axis)
        elif isinstance(row_axis, nb.cifti2.cifti2_axes.Axis):
            pass
        else:
            raise ValueError(
                "The input row_axis instance type does not meet the requirement!"
            )

        bm = self.get_brain_model_axis()
        header = nb.Cifti2Header.from_axes((row_axis, bm))
        cifti_img = nb.Cifti2Image(dataobj=data, header=header)

        return cifti_img

    def read_data(self, img, interpolation=0):
        """reads data for surface-based atlas from list of gifti [left,right]or single cifti file. Adjusts automatically for node masks.

        Args:
            img (nibabel.image): Cifti or (list of) gifti images
            interpolation (int): nearest neighbour (0), trilinear (1)
        Returns:
            data (ndarray): (N,P) ndarray
        """
        if isinstance(img, nb.Cifti2Image):
            data = self.cifti_to_data(img)
        elif isinstance(img, list):
            if len(img) != len(self.structure):
                raise (NameError("Number of images needs to match len(self.structure)"))
            data = []
            for i, im in enumerate(img):
                if isinstance(im, str):
                    im = nb.load(im)
                d = im.agg_data()
                data.append(d[self.vertex[i]])
            data = np.hstack(data)
        return data

    def cifti_to_data(self, cifti):
        """Gets the data from a CIFTI file, checking the structure name
           and vertex information in the cifti file. If it doesn't match
           the vertex information in the atlas object, it gives a warning,
           but corrects for it by extracting the available data.

        Args:
            cifti (ciftiimage or filename): Cifti file to be used
        Returns:
            data (np.ndarray): NxP np-array
        """
        # First check the input is a cifti image object
        if isinstance(cifti, str):
            cifti = nb.load(cifti)

        data = cifti.get_fdata()
        # Second check if the data has the same vertices index or
        # the number of vertices with the current atlasMap
        col_axis = cifti.header.get_axis(1)

        return_data = []
        img_stru = list(col_axis.iter_structures())
        img_names = [n[0] for n in img_stru]
        for i, stru in enumerate(self.structure):
            # Align the structure name to cifti file
            stru = nb.cifti2.BrainModelAxis.to_cifti_brain_structure_name(stru)

            if stru not in img_names:
                print(f"The input image does not contain {stru}! (Fill with NaN)")
                # if the input image doesn't have current brain structure
                # we then fill these vertices with NaN value
                this_data = np.full((data.shape[0], self.vertex[i].shape[0]), np.nan)
            else:
                idx = img_names.index(stru)
                data_idx = img_stru[idx][1]
                bm = img_stru[idx][2]

                # Restore the full data for this structure
                this_full_data = np.full((data.shape[0], bm.nvertices[stru]), np.nan)
                this_full_data[:, bm.vertex] = data[:, data_idx]
                this_data = this_full_data[:, self.vertex[i]]

            return_data.append(this_data)

        return np.hstack(return_data)

    def get_brain_model_axis(self):
        """Returns brain model axis

        Returns:
            bm (cifti2.BrainModelAxis)
        """
        # Make the brain Structure models
        for i, name in enumerate(self.structure):
            if i == 0:
                bm = nb.cifti2.BrainModelAxis.from_mask(
                    self.vertex_mask[i], name=self.structure[i]
                )
            else:
                bm = bm + nb.cifti2.BrainModelAxis.from_mask(
                    self.vertex_mask[i], name=self.structure[i]
                )
        return bm

    def get_parcel_axis(self):
        """Returns parcel axis based on the current setting of labels

        Returns:
            bm (cifti2.ParcelAxis)
        """
        bm_list = []
        if self.unite_struct:
            raise (
                NameError(
                    "Cannot create parcel axis with ROIs spanning both hemispheres. Set unite_struct to FALSE."
                )
            )
        if not hasattr(self, "labels"):
            raise (
                NameError(
                    "Atlas has no parcels defined yet - call atlas.getparcels(image) first"
                )
            )
        indx = 0
        for i, struct_name in enumerate(self.structure):
            n_vert = self.vertex[i].shape[0]
            label_vec = self.label_vector[indx : indx + n_vert]
            for h, label in enumerate(self.labels):
                pindx = label_vec == label
                if pindx.sum() > 0:
                    pname = struct_name + f"_{label:02d}"
                    bm = nb.cifti2.BrainModelAxis.from_surface(
                        self.vertex[i][pindx],
                        self.vertex_mask[i].shape[0],
                        name=struct_name,
                    )
                    bm_list.append((pname, bm))
        parcel_axis = nb.cifti2.ParcelsAxis.from_brain_models(bm_list)
        return parcel_axis

    def get_parcel(self, label_img, unite_struct=False):
        """adds a label_vec to the atlas that assigns each voxel / node of the atlas
        label_vec[i] = 0 mean that the ith voxel / node is unassigned

        Args:
            label_img (str/list) - filename(s) of label image(s)
            unite_struct (bool) - join labels across hemispheres? Default: False
        Returns:
            label_vec (list) - list of numpy array containing label values corresponding to label images
            labels (ndarray) - List of unique labels > 0
        """

        # get label vectors for each label image
        self.label_vector = self.read_data(label_img, 0).astype(int)
        self.labels = np.unique(self.label_vector[self.label_vector > 0])
        self.n_labels = self.labels.shape[0]

        # change the values of non zero labels if necessary (only for the second hemi)
        self.unite_struct = unite_struct
        self.label_list = []
        if not self.unite_struct:
            indx = 0
            for i, s in enumerate(self.structure):
                n_vert = self.vertex[i].shape[0]
                label_vec = self.label_vector[indx : indx + n_vert]
                label_vec[label_vec > 0] += self.n_labels * i
                self.label_list.append(label_vec)
                indx = indx + n_vert
        self.labels = np.unique(self.label_vector[self.label_vector > 0])
        self.n_labels = self.labels.shape[0]
        return self.label_vector, self.labels



class AtlasSurfaceSymmetric(AtlasSurface):
    """
    """

    def __init__(self, name, mask_gii, structure = ["cortex_left","cortex_right"]):
        """ Surface atlas with left-right symmetry
        The atlas behaves like AtlasSurface, but provides
        a reduced (symmetric) representation of size Psym.AtlasSurfaceSymmeytric
        Constructor: Generates members indx_full, indx_reduced, indx_flip.
        mapping indices from a full representation to

        | Assume you have a
        | Full: N x P array
        | Left: N x Psym array
        | Right: N x Psym array
        | then:
        | Left = Full[:,index_full[0]]
        | Right = Full[:,index_full[1]]
        | Avrg = (Left + Right)/2
        | Full = Avrg[:,index_reduced]
        | To Flip:
        | flippedFull = Full[:,index_flip]

        Args:
            name (str): Name of the brain structure (cortex_left, cortex_right, cerebellum)
            mask_gii (list): gifti file name of mask image defining atlas locations
            structure (list): [cortex_left, cortex_right] - CIFTI brain structure names
        """
        super().__init__(name, mask_gii, structure)
        assert np.array_equal(
            self.vertex[0], self.vertex[1]
        ), "The left and right hemisphere must be symmetric!"

        # Initialize indices
        # This is the number of vertices before masking
        n_vertex = self.vertex_mask[0].shape[0]
        self.Psym = int(self.P / 2)
        self.indx_full = np.zeros((2, self.Psym), dtype=int)
        # n_vertex = self.vertex[0].shape[0]

        # Generate full/reduce index
        self.indx_full[0, :] = np.arange(self.Psym)
        self.indx_full[1, :] = np.arange(self.Psym) + self.Psym
        self.indx_reduced = np.tile(np.arange(self.Psym), 2)

        # Generate flipping index
        indx_orig = np.arange(self.P, dtype=int)
        self.indx_flip = np.zeros((self.P,), dtype=int)
        self.indx_flip[self.indx_full[0]] = indx_orig[self.indx_full[1]]
        self.indx_flip[self.indx_full[1]] = indx_orig[self.indx_full[0]]


class AtlasMapDeform:
    def __init__(self, world, deform_img, mask_img):
        """AtlasMapDeform stores the mapping rules for a non-linear deformation
        to the desired atlas space in form of a voxel list

        Args:
            worlds (ndarray): 3xP ND array of world locations
            deform_img (str/list): Name of deformation map image(s)
            mask_img (str): Name of masking image that defines the functional data space.
        """
        self.P = world.shape[1]
        self.world = world
        if type(deform_img) is not list:
            deform_img = [deform_img]
        self.deform_img = []
        for di in deform_img:
            self.deform_img.append(nb.load(di))
        self.mask_img = nb.load(mask_img)

    def build(self, smooth=None, additional_mask=None):
        """ Using the dataset, builds a list of voxel indices of
        For each of the locations. It creates:
        vox_list: List of voxels to sample for each atlas location
        vox_weight: Weight of each of these voxels to determine the atlas location

        Args:
            smooth (double): SD of smoothing kernel (mm) or None for nearest neighbor
            additional_mask: Additional Mask image (not necessarily in functional space - only voxels with elements > 0 in that image
            will be used for the altas )
        """
        # Caluculate locations of atlas in individual (deformed) coordinates
        # Apply possible multiple deformation maps sequentially
        xyz = self.world.copy()
        for i, di in enumerate(self.deform_img):
            xyz = nt.sample_image(di, xyz[0], xyz[1], xyz[2], 1).squeeze().T
            pass
        atlas_ind = xyz
        N = atlas_ind.shape[1]  # Number of locations in atlas

        # Determine which voxels are available in functional space
        # and apply additional mask if given
        M = self.mask_img.get_fdata()
        i, j, k = np.where(M > 0)
        vox = np.vstack((i, j, k))
        # available voxels in world coordiantes
        world_vox = nt.affine_transform_mat(vox, self.mask_img.affine)
        if additional_mask is not None:
            # If file name, load the nifti image
            if isinstance(additional_mask, str):
                additional_mask = nb.load(additional_mask)
            add_mask = nt.sample_image(
                additional_mask, world_vox[0], world_vox[1], world_vox[2], 1
            )
            world_vox = world_vox[:, add_mask > 0]
            vox = vox[:, add_mask > 0]

        if smooth is None:  # Use nearest neighbor interpolation
            linindx, good = nt.coords_to_linvidxs(atlas_ind, self.mask_img, mask=True)
            self.vox_list = linindx.reshape(-1, 1)
            self.vox_weight = np.ones((linindx.shape[0], 1))
            self.vox_weight[np.logical_not(good)] = np.nan
        else:  # Use smoothing kernel of specific size
            linindx = np.ravel_multi_index(
                (vox[0, :], vox[1, :], vox[2, :]), M.shape, mode="clip"
            )
            # Distances between atlas coordinates and voxel coordinates
            # TODO: For a lot of voxels, calculating the Euclidian distance is very memory hungry. Build the atlas iteratively to avoid running into memory issues.
            D = nt.euclidean_dist_sq(atlas_ind, world_vox)
            # Find voxels with substantial power under gaussian kernel
            W = np.exp(-0.5 * D / (smooth**2))
            W[W < 0.2] = 0
            a, b = W.nonzero()
            # Now transfer them into a full list of voxels
            # this is somewhat ugly and brute force
            c = np.zeros(a.shape, dtype=int)
            c[1:] = a[0:-1] == a[1:]
            for i in range(c.shape[0] - 1):
                if c[i + 1]:
                    c[i + 1] = c[i + 1] + c[i]
            self.vox_list = np.zeros((N, c.max() + 1), dtype=np.int32)
            self.vox_weight = np.zeros((N, c.max() + 1))
            self.vox_list[a, c] = linindx[b]
            self.vox_weight[a, c] = W[a, b]
            # Avoid divide by zero error:
            mw = self.vox_weight.sum(axis=1, keepdims=True)
            mw[mw == 0] = np.nan
            self.vox_weight = self.vox_weight / mw
        pass


class AtlasMapSurf:
    def __init__(self, vertex, white_surf, pial_surf, mask_img):
        """AtlasMapSurf stores the mapping rules for a freesurfer-style surface (pial + white surface pair)

        Args:
            vertex (ndarray): Array of vertices included in the atlas map
            white_surf (str): Name for white matter surface
            pial_surf (str): Name for pial surface
            mask_img (str): Name of masking image that defines the functional data space.
        """
        self.P = len(vertex)
        self.vertex = vertex
        self.white_surf = nb.load(white_surf)
        self.pial_surf = nb.load(pial_surf)
        self.mask_img = nb.load(mask_img)

    def build(self, smooth=None, depths=[0, 0.2, 0.4, 0.6, 0.8, 1.0]):
        """ Using the dataset, builds a list of voxel indices of
        each of the nodes

            | vox_list: List of voxels to sample for each atlas location
            | vox_weight: Weight of each of these voxels to determine the atlas location

        Args:
            smooth (double): SD of smoothing kernel (mm) or None for nearest neighbor
            depths (iterable): List of depth between pial (1) and white (0) surface that
            will be sampled
        """
        n_points = len(depths)
        c1 = self.white_surf.darrays[0].data[self.vertex, :].T
        c2 = self.pial_surf.darrays[0].data[self.vertex, :].T
        n_vert = c1.shape[1]
        if c2.shape[1] != n_vert:
            raise (
                NameError(
                    "White and pial surfaces should have same number of vertices."
                )
            )

        # Get the indices for all the points being sampled
        indices = np.zeros((n_points, 3, n_vert))
        for i in range(n_points):
            indices[i, :, :] = (1 - depths[i]) * c1 + depths[i] * c2

        self.vox_list, good = nt.coords_to_linvidxs(indices, self.mask_img, mask=True)
        all = good.sum(axis=0)
        # print(f'{self.name} has {np.sum(all==0)} vertices without data')
        all[all == 0] = 1
        self.vox_weight = good / all
        self.vox_list = self.vox_list.T
        self.vox_weight = self.vox_weight.T


def get_data_nifti(fnames, atlas_maps):
    """Extracts the data for a list of fnames
    for a list of atlas_maps. This is usually called by DataSet.extract_data()
    to extract the required raw data before processing it further

    Args:
        fnames (list): List of file names to be sampled
        atlas_maps (list): List of K atlas-map or atlas objects
    returns:
        data (list): List of NxP_k 2-d array data matrices (np)
    """
    n_atlas = len(atlas_maps)

    # Deal with 4-d vs. 3d-nifti
    vols = []
    for j, f in enumerate(fnames):
        if isinstance(f, str):
            V = nb.load(f)
        else:
            V = f
        if V.ndim > 3:
            vols = vols + nb.funcs.four_to_three(V)
        else:
            vols.append(V)

    n_vols = len(vols)
    data = []
    # Make the empty data structures
    for at in atlas_maps:
        data.append(np.full((n_vols, at.P), np.nan))
    for j, V in enumerate(vols):
        X = V.get_fdata()
        X = X.ravel()
        for i, at in enumerate(atlas_maps):
            d = X[at.vox_list] * at.vox_weight  # Expanded data
            d = np.nansum(d, axis=1)
            d[np.nansum(at.vox_weight, axis=1) == 0] = np.nan
            data[i][j, :] = d
    return data


def get_data_cifti(fnames, atlases):
    """Extracts the data for a list of fnames
    for a list of atlas_maps. This is usually called by DataSet.get_data()
    to extract the required raw data before processing it further

    Args:
        fnames (list): list of file names to be sampled
        atlas_maps (list): list of K built atlas-map objects
    returns:
        data (list): List of NxP_k 2-d array data matrices (np)
    """
    n_atlas = len(atlases)
    data = [[]] * n_atlas
    # Make the empty data structures
    # Loop over files
    for j, f in enumerate(fnames):
        cifti = nb.load(f)
        for i, at in enumerate(atlases):
            if isinstance(at, AtlasMapDeform):
                V = nt.volume_from_cifti(cifti, ["cerebellum"])
                data[i].append(get_data_nifti([V], [at])[0])
            elif isinstance(at, AtlasVolumetric):
                V = nt.volume_from_cifti(cifti, [at.structure])
                data[i].append(
                    nt.sample_image(
                        V, at.world[0, :], at.world[1, :], at.world[2, :], 0
                    )
                )
            elif isinstance(at, AtlasSurface):
                data[i].append(at.cifti_to_data(cifti))
    for i in range(n_atlas):
        data[i] = np.vstack(data[i])
    return data
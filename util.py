import numpy as np 
from numpy.linalg import inv
import nibabel as nb

def affine_transform(x1, x2, x3, M):
    """
    Returns affine transform of x

    Args:
        x1 (np-array):
            X-coordinate of original
        x2 (np-array):
            Y-coordinate of original
        x3 (np-array):
            Z-coordinate of original
        M (2d-array):
            4x4 transformation matrix

    Returns:
        x1 (np-array):
            X-coordinate of transform
        x2 (np-array):
            Y-coordinate of transform
        x3 (np-array):
            Z-coordinate of transform

    """
    y1 = M[0,0]*x1 + M[0,1]*x2 + M[0,2]*x3 + M[0,3]
    y2 = M[1,0]*x1 + M[1,1]*x2 + M[1,2]*x3 + M[1,3]
    y3 = M[2,0]*x1 + M[2,1]*x2 + M[2,2]*x3 + M[2,3]
    return (y1,y2,y3)

def affine_transform_mat(x, M):
    """
    Returns affine transform of x

    Args:
        x (np-array):
            3xN array for original coordinates 
        M (2d-array):
            4x4 transformation matrix
    Returns:
        y (np-array):
            3xN array pof X-coordinate of transformed coordinaters 
    """
    y = M[0:3,0:3] @ x + M[0:3,3:]
    return (y)


def coords_to_linvidxs(coords,vol_def,mask=False):
    """
    Maps coordinates to linear voxel indices

    INPUT:
        coords (3xN matrix or Qx3xN array):
            (x,y,z) coordinates
        vol_def (nibabel object):
            Nibabel object with attributes .affine (4x4 voxel to coordinate transformation matrix from the images to be sampled (1-based)) and shape (1x3 volume dimension in voxels)
        mask (bool):
            If true, uses the mask image to restrict voxels (all outside = -1)
    OUTPUT:
        linvidxs (N-array or QxN matrix):
            Linear voxel indices
        good (bool) boolean array that tells you whether the index was in the mask
    """
    mat = inv(vol_def.affine)

    # Check that coordinate transformation matrix is 4x4
    if (mat.shape != (4,4)):
        raise(NameError('Error: Matrix should be 4x4'))

    rs = coords.shape
    if (rs[-2] != 3):
        raise(NameError('Coords need to be a (Kx) 3xN matrix'))

    # map to 3xP matrix (P coordinates)
    ijk = mat[0:3,0:3] @ coords + mat[0:3,3:]
    ijk = np.rint(ijk).astype(int)
    
    if ijk.ndim<=2:
        i = ijk[0]
        j = ijk[1]
        k = ijk[2]
    elif ijk.ndim==3:
        i = ijk[:,0]
        j = ijk[:,1]
        k = ijk[:,2]

    # Now set the indices out of range to 
    good = (i>=0) & (i<vol_def.shape[0]) & (j>=0) & (j<vol_def.shape[1]) &  (k[2]>=0) & (k[2]<vol_def.shape[2])
    
    linindx = np.ravel_multi_index((i,j,k),vol_def.shape,mode='clip')

    if mask:
        M=vol_def.get_fdata().ravel()
        good = good & (M[linindx]>0)
    
    return linindx, good

def sq_eucl_distances(coordA,coordB): 
    D = coordA.reshape(3,-1,1)-coordB.reshape(3,1,-1)
    D = np.sum(D**2,axis=0)
    return D


def volume_from_cifti(ts_cifti, struct_names=None):
        """
        Gets the 4D nifti object containing the time series
        for all the subcortical structures
        Args:
            ts_cifti (cifti obj ) - cifti object of the time series
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

            # if (struct_names is None) | (nam in struct_names): 

            # get the voxels/vertices corresponding to the current brain model
            ijk = bm.voxel
            bm_data = ts_array[:, slc]
            i  = (ijk[:,0] > -1)

            # fill in data
            subcorticals_vol[ijk[i, 0], ijk[i, 1], ijk[i, 2], :]=bm_data[:,i].T

        # save as nii
        nii_vol_4d = nb.Nifti1Image(subcorticals_vol,bmf.affine)
        # if save:
        #     ts_nifti = dir+'/sub-100307_ses-01_task-rest_space-subcortex_run-01_bold.nii'
        #     nb.save(nii_vol,ts_nifti)
        return nii_vol_4d

def surf_from_cifti(ts_cifti,
                    struct_names=['CIFTI_STRUCTURE_CORTEX_LEFT',
                                    'CIFTI_STRUCTURE_CORTEX_RIGHT']):
        """
        Gets the time series of cortical surface vertices (Left and Right)
        Args:
            ts_cifti (cifti obj) - cifti object of time series
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
            if nam in struct_names: 
                values = np.full((ts_array.shape[0],bmf.nvertices[nam]),np.nan)
                # get the values corresponding to the brain model
                values[:,bm.vertex] = ts_array[:, slc]
                ts_list.append(values)
            else:
                break
        return ts_list

def zstandarize_ts(X):
    X = X - X.mean(axis = 0, keepdims = True)
    X = X / np.sqrt(np.nansum(X**2, axis=0)/X.shape[0])
    return X 

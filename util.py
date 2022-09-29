import numpy as np
from numpy.linalg import inv
import nibabel as nb


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

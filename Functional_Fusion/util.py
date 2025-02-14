import numpy as np
from numpy.linalg import inv, pinv
import nibabel as nb
import h5py, os, subprocess
import Functional_Fusion.atlas_map as am
from pathlib import Path

default_atlas_dir = os.path.dirname(am.__file__) + '/Atlases'

def get_base_dir():
    possible_dirs = ['/Volumes/diedrichsen_data$/data/FunctionalFusion',
                     '/srv/diedrichsen/data/FunctionalFusion',
                     'Y:/data/FunctionalFusion']
    for directory in possible_dirs:
        if Path(directory).exists():
            return directory
    raise FileNotFoundError('Could not find base_dir')
    
def sq_eucl_distances(coordA,coordB):
    D = coordA.reshape(3,-1,1)-coordB.reshape(3,1,-1)
    D = np.sum(D**2,axis=0)
    return D

def nan_linear_model(X,Y,unknowns_to_nan = True):
    """ Calculates the Nansafe estimates of the regression model Y on design matrix X
    Y = X @ B
    If an entire row in Y is missing, the B-estimate will skip that observation
    If this there is no information in the data about that B anymore, the B's are set to NaN

    Args:
        X (ndarray): N x Q Design matrix
        Y (ndarray): N x P data matrix
        unknowns_to_nan (bool): Set regression coefficients without information to Nan (instead of to zero)
    Returns:
        B (ndarray): Q x P regression weights - non-estimable weights are set to NaN
    """
    indx = np.logical_not(np.isnan(Y).all(axis=1))
    B = pinv(X[indx,:]) @ Y[indx,:]
    if unknowns_to_nan:
        unknown = np.abs(X[indx,:]).sum(axis=0)==0
        B[unknown,:]=np.nan
    return B


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
                                  'CIFTI_STRUCTURE_CORTEX_RIGHT'],
                    mask_gii=None):
        """Gets the time series of cortical surface vertices (Left and Right)

        Args:
            ts_cifti (cifti obj) - cifti object of time series
            struct_names (list): the struct name of left and right cortex
            mask_gii: (list of Obj.) the mask gifti object for each hemisphere
                                     if given. Default is None,
                                     indicating no mask for return.
                      (list of str.) the file path of mask gifti for each
                                     hemisphere if given. Default is None,
                                     indicating no mask for return.
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
                if mask_gii is not None:
                    Xmask = mask_gii[struct_names.index(nam)]
                    if isinstance(Xmask, str):
                        Xmask = nb.load(Xmask).agg_data()
                    elif isinstance(Xmask, object):
                        Xmask = Xmask.agg_data()
                    else:
                        raise ValueError("The input mask_gii must be either a string of path"
                                         "or a nibabel gii image object!")
                    values = values[:, np.where(Xmask>0)[0]]
                ts_list.append(values)
            else:
                break
        return ts_list

def zstandarize_ts(X):
    X = X - X.mean(axis = 0, keepdims = True)
    X = X / np.sqrt(np.nansum(X**2, axis=0)/X.shape[0])
    return X

def correlate(X, Y):
    """ Correlate X and Y numpy arrays after standardizing them"""
    X = zstandarize_ts(X)
    Y = zstandarize_ts(Y)
    return Y.T @ X / X.shape[0]

def templateflow_xfm_h5_to_nii(filename):
    """This is a first attempt to understand the _xfm.h5 files from templateflow
    and convert them to nifti files. Not easy without documentation - so this is work
    in progress for now

    Args:
        filename (str): h5-filename
    """
    with h5py.File(filename, 'r') as f:
        P2 = np.array(f['TransformGroup']['2']['TranformParameters'])
        P2f = np.array(f['TransformGroup']['2']['TranformFixedParameters'])
        P1 = np.array(f['TransformGroup']['1']['TranformParameters'])
        P1f = np.array(f['TransformGroup']['1']['TranformFixedParameters'])

        # Wild guess on the paramaters for the second transform
        dim = P2f[:3].astype(int)
        trans = P2f[3:6]
        voxsize = P2f[6:9]
        rot = P2f[9:].reshape(3,3)
        affine = np.zeros((4,4))
        affine[:3,:3] = rot
        affine[:3,3] = trans
        affine[3,3] = 1
        P2 = P2.reshape(np.r_[dim,(3)])
        deform_img = nb.Nifti1Image(P2,affine)
        nb.save(deform_img, filename.replace('.h5','.nii'))
        pass


def get_volumes(data, atlas_name='MNISymC2'):
    """
    Projects CIFTI data from any space (e.g., SUIT3/MNISymC2) back to volume space.

    Args:
        data (list, np.ndarray): Input data. Can be a list, NumPy array.
        atlas (str): Atlas code to specify the projection space (e.g., 'SUIT3', 'MNISymC3'). Default is 'MNISymC2'.

    Returns:
        list: A list of NIfTI volume objects projected from the input data.
    """

    atlas,_ = am.get_atlas(atlas_name)

    # Determine number of volumes based on data type
    if isinstance(data, np.ndarray):
        n_vols = data.shape[0]
    elif isinstance(data, list):
        n_vols = len(data)
    else:
        raise TypeError("data must be a list, np.ndarray, or torch.Tensor")

    nii_vols = []
    
    # Project each data volume back to NIfTI volume space
    for i in range(n_vols):
        cifti_data = data[i]
        vol = atlas.data_to_nifti(cifti_data)
        nii_vols.append(vol)

    return nii_vols


def split_string(s):
    """ Split string into integer and string

    Args:
        s: input string to be splitted

    Returns:
        integer part, and string part
    """
    # Split string into integer and string
    for i, char in enumerate(s):
        if not char.isdigit():
            return int(s[:i]), s[i:]


def smooth_fs32k_data(input_file, smooth=1, kernel='gaussian',
                      return_data_only=False):
    """ Smooth cortical data in fs32k space from cifti file, then save as cifti
        format using workbench command

    Args:
        input_file: the input cifti file to be smoothed
        smooth: width of smoothing kernel in mm
        kernel: the type of smoothing kernel, either "gaussian" or "fwhm"
        return_data_only: only return the smoothed data but not store any

    Returns:
        Write in the cifti file of smoothed data
    """
    # get the surfaces for smoothing
    surf_L = default_atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_sphere.surf.gii'
    surf_R = default_atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_sphere.surf.gii'

    # Making in / out file names
    dir_path, file_name = os.path.split(input_file)
    base_name = file_name.split('.')[0]
    ext = '.' + '.'.join(file_name.split('.')[1:])

    if kernel == 'gaussian':
        smooth_suffix = f'_desc-sm{smooth}'
    elif kernel == 'fwhm':
        smooth_suffix = f'_desc-sm{smooth}fwhm'
    else:
        raise ValueError('Only gaussian and fwhm kernels are supported!')

    cifti_out = os.path.join(dir_path, base_name + smooth_suffix + ext)

    # Load data and fill nan with zeros if unsmoothed data contains any
    contain_nan = False  ## a flag
    C = nb.load(input_file)
    if np.isnan(C.get_fdata()).any():
        contain_nan = True
        mask = np.isnan(C.get_fdata())
        C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
        input_file = dir_path + f'/tmp' + ext
        nb.save(C, input_file)

    # Write in smoothed surface data (filled with 0)
    smooth_cmd = f"wb_command -cifti-smoothing {input_file} " \
                 f"{smooth} {smooth} COLUMN {cifti_out} " \
                 f"{f'-{kernel} ' if kernel == 'fwhm' else ''}" \
                 f"-left-surface {surf_L} -right-surface {surf_R} " \
                 f"-fix-zeros-surface"
    subprocess.run(smooth_cmd, shell=True)

    # Double-check if the original data contain NaN values
    if contain_nan:
        os.remove(dir_path + f'/tmp' + ext)
        # Replace 0s back to NaN (we don't want the 0s impact model learning)
        data = nb.load(cifti_out).get_fdata()
        data[mask] = np.nan
        C = nb.Cifti2Image(dataobj=data, header=C.header)
        nb.save(C, cifti_out)

    if return_data_only:
        os.remove(cifti_out)
        return data
    

def align_conditions(Ya, Yb, info_a, info_b):
    """
    Align two datasets based on shared conditions, align all conditions to the mean of shared conditions,
    then average shared conditions and append unique conditions.

    Args:
    Ya (numpy array): Dataset A (subjects x conditions x voxels) or (conditions x voxels)
    Yb (numpy array): Dataset B (subjects x conditions x voxels) or (conditions x voxels)
    info_a (pandas.DataFrame): Info file for Dataset A
    info_b (pandas.DataFrame): Info file for Dataset B

    Returns:
    combined_data (numpy array): Combined dataset
    combined_info (pandas.DataFrame): Combined info
    """

    shared_conditions = np.intersect1d(info_a['cond_code'], info_b['cond_code'])
    if len(shared_conditions) == 0:
        raise ValueError("No shared conditions between datasets.")

    # Standardize input dimensions to 3D if needed
    if len(Ya.shape) == 2:
        Ya = Ya[None, :, :]
    if len(Yb.shape) == 2:
        Yb = Yb[None, :, :]

    # Sort shared conditions and get indices
    shared_sorted = np.sort(shared_conditions)
    order_a = []  # Indices for shared conditions in Dataset A
    order_b = []  # Indices for shared conditions in Dataset B

    # Loop through each condition in the sorted shared conditions array
    for cond in shared_sorted:
        # Find the index of the condition in info_a that matches the current condition code
        idx_a = info_a[info_a['cond_code'] == cond].index[0]
        # Find the index of the condition in info_b that matches the current condition code
        idx_b = info_b[info_b['cond_code'] == cond].index[0]

        # Append the indices to the respective lists
        order_a.append(idx_a)
        order_b.append(idx_b)

    # Align shared conditions
    Ya_shared = Ya[:, order_a, :]
    Yb_shared = Yb[:, order_b, :]
    Ya_mean, Yb_mean = Ya_shared.mean(1, keepdims=True), Yb_shared.mean(1, keepdims=True)
    Ya_aligned, Yb_aligned = Ya - Ya_mean, Yb - Yb_mean

    # Average aligned shared conditions
    shared_avg = (Ya_aligned[:, order_a, :] + Yb_aligned[:, order_b, :]) / 2.0

    # Combine shared and unique conditions
    unique_a = np.setdiff1d(info_a['cond_code'], shared_sorted)
    unique_a_indices = info_a['cond_code'].isin(unique_a)
    unique_b = np.setdiff1d(info_b['cond_code'], shared_sorted)
    unique_b_indices = info_b['cond_code'].isin(unique_b)

    Ya_aligned_unique = Ya_aligned[:, unique_a_indices, :]
    Yb_aligned_unique = Yb_aligned[:, unique_b_indices, :]
    combined_data = np.concatenate([shared_avg, Ya_aligned_unique,
                                    Yb_aligned_unique], axis=1)

    # Create combined info file
    shared_info = info_a.loc[order_a, ['cond_name', 'cond_code']].copy()
    shared_info['source'] = 'averaged'
    unique_info = pd.concat([info_a[unique_a_indices], info_b[unique_b_indices]])
    unique_info['source'] = 'Novel'
    combined_info = pd.concat([shared_info, unique_info[['cond_name', 'cond_code', 'source']]], ignore_index=True)

    if Ya.shape[0] == 1:
        combined_data = combined_data[0]
    else:
        combined_data = combined_data


    return combined_data, combined_info
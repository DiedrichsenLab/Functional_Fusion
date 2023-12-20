import nibabel as nib
from pathlib import Path
import os.path as op
import preprocessing.paths as paths
import Functional_Fusion.util as util

# create seed-based correlation map
# extract voxel timeseries in native space
def extract_voxel_timeseries(data, coord):
    x, y, z = coord
    voxel_timeseries = data[x, y, z, :]
    return voxel_timeseries

# correlate with all other voxels in native space
def correlate_with_all_voxels(X, Y):
    X = util.zstandarize_ts(X)
    Y = util.zstandarize_ts(Y)
    return Y.T @ X / X.shape[0]

# z-transform correlation map
def z_transform_correlation_map(correlation_map):
    # TODO: Implement the z-transform of the correlation map
    pass
    return correlation_map


if __name__ == "__main__":
    # set paths
    base_dir = paths.set_base_dir()
    rest_dir = f'{base_dir}/../Cerebellum/super_cerebellum/resting_state/'
    subject = 's10'
    run = 1
    figure_dir = paths.set_figure_dir()

    

    # Load data
    data_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/filtered_func_data.nii.gz'
    data = nib.load(data_path).get_fdata()

    # Extract voxel timeseries
    seed_coord = (35, 25, 21)
    voxel_timeseries = extract_voxel_timeseries(data, seed_coord)
    # voxel_timeseries = voxel_timeseries.reshape((voxel_timeseries.shape[0], 1)).T

    # Load mask that excludes non-brain and non-skull background voxels
    mask_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/mask.nii.gz'
    mask = nib.load(mask_path).get_fdata()
    data = data[mask > 0]


    # Correlate with all other voxels
    correlation_map = correlate_with_all_voxels(data.T, voxel_timeseries)

    # Z-transform correlation map
    # z_transformed_map = z_transform_correlation_map(correlation_map)

    # Save correlation map
    save_path = f'{rest_dir}/fix_ica/correlation_map.nii.gz'
    correlation_img = mask.copy()
    correlation_img[mask > 0] = correlation_map
    correlation_image = nib.Nifti1Image(correlation_img, affine=nib.load(mask_path).affine)
    # Save nibabel image
    nib.save(correlation_image, save_path)
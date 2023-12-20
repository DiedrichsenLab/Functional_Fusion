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

coordinates = {
    's02': (36, 23, 20),
    's03': (38, 21, 17),
    's06': (37, 21, 18),
    's07': (39, 24, 19),
    's08': (38, 23, 19),
    's09': (39, 24, 22),
    's10': (38, 26, 21),
    's12': (40, 25, 20),
    # 's13': (37, 22, 19),
    # 's14': (37, 22, 19),
    # 's15': (37, 22, 19),
    # 's16': (37, 22, 19),
    # 's17': (37, 22, 19),
    # 's18': (37, 22, 19),
    # 's19': (37, 22, 19),
    # 's20': (37, 22, 19),
    # 's21': (37, 22, 19),
    # 's22': (37, 22, 19),
    # 's23': (37, 22, 19),
    # 's24': (37, 22, 19),
    # 's25': (37, 22, 19),
    # 's26': (37, 22, 19),
    # 's27': (37, 22, 19),
    # 's28': (37, 22, 19),
    # 's29': (37, 22, 19),
    # 's30': (37, 22, 19),
    # 's31': (37, 22, 19),

}

subjects_runs = {
    's06': [1],
    's07': [1],
    's08': [2],
    's10': [2]}
if __name__ == "__main__":
    # set paths
    base_dir = paths.set_base_dir()
    rest_dir = f'{base_dir}/../Cerebellum/super_cerebellum/resting_state/'
    # for s, subject in enumerate(rest_dir.glob('s*')):
    for s, subject in enumerate(subjects_runs.keys()):
        # for run in [1, 2]:
        for run in subjects_runs[subject]:
            for clean in [True, False]:
                # Load data
                data_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/filtered_func_data.nii.gz'
                if clean:
                    data_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/filtered_func_data_clean.nii.gz'
                data = nib.load(data_path).get_fdata()

                # Extract voxel timeseries
                seed_coord = coordinates[subject]
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
                save_path = f'{rest_dir}/fix_ica/correlation_map_{subject}-run{run:02d}.nii.gz'
                if clean:
                    save_path = f'{rest_dir}/fix_ica/correlation_map_{subject}-run{run:02d}_clean.nii.gz'
                correlation_img = mask.copy()
                correlation_img[mask > 0] = correlation_map
                correlation_image = nib.Nifti1Image(correlation_img, affine=nib.load(mask_path).affine)
                # Save nibabel image
                nib.save(correlation_image, save_path)

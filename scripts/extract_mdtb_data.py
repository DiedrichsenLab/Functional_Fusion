# Script for importing the MDTB data set from super_cerebellum to general format.
import os, subprocess, time
import pandas as pd
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import scripts.fusion_paths as paths
import Functional_Fusion.atlas_map as am
import Functional_Fusion.connectivity as conn
import Functional_Fusion.util as fut

from pathlib import Path
from Functional_Fusion.dataset import DataSetMDTB

base_dir = paths.set_base_dir()
atlas_dir = paths.set_atlas_dir(base_dir)
dname = 'MDTB'
data_dir = paths.set_fusion_dir(base_dir)
mdtb_dir = data_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'


def group_mtdb(ses_id='ses-s1', type='CondHalf', atlas='SUIT3'):
    mdtb_dataset = DataSetMDTB(mdtb_dir)
    mdtb_dataset.group_average_data(ses_id, type, atlas)


def parcel_mdtb_fs32k(res=162, ses_id='ses-s1', type='condHalf'):
    # Make the atlas object
    surf_parcel = []
    hem_name = ['cortex_left', 'cortex_right']
    # Get the parcelation
    for i, h in enumerate(['L', 'R']):
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i], gifti))

    # initialize the data set object
    mdtb_dataset = DataSetMDTB(mdtb_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    for s in T.participant_id[0:2]:
        print(f'Average {s}')
        s_dir = mdtb_dataset.base_dir.format(s)
        C = nb.load(s_dir + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
        bmf = C.header.get_axis(1)
        bmp = []
        R = []
        for idx, (nam, slc, bm) in enumerate(bmf.iter_structures()):
            D = np.asanyarray(C.dataobj[:, slc])
            X = np.zeros((D.shape[0], surf_parcel[0].label_vec.shape[0]))
            X[:, bm.vertex] = D
            R.append(surf_parcel[idx].agg_data(X))
            bmp.append(surf_parcel[idx].get_parcel_axis())
        header = nb.Cifti2Header.from_axes(
            (C.header.get_axis(0), bmp[0] + bmp[1]))
        cifti_img = nb.Cifti2Image(dataobj=np.c_[R[0], R[1]], header=header)
        nb.save(cifti_img, s_dir +
                f'/{s}_space-fs32k_{ses_id}_{type}_Iso-{res}.pscalar.nii')
        pass


def smooth_mdtb_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1, kernel='gaussian'):
    myatlas, _ = am.get_atlas('fs32k')
    mdtb_dataset = DataSetMDTB(mdtb_dir)
    T = mdtb_dataset.get_participants()

    # get the surfaces for smoothing
    surf_L = mdtb_dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.L.midthickness.surf.gii'
    surf_R = mdtb_dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.R.midthickness.surf.gii'

    for s in T.participant_id:
        print(f'Smoothing data for {s} fs32k {ses_id} in {smooth}mm {kernel} ...')

        start = time.perf_counter()
        file = mdtb_dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        fut.smooth_fs32k_data(file, surf_L, surf_R, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")

def reshape_data(data, info, cond_column='cond_num_uni', part_column='run', mean_centering=False):
    """Reshape data from (n_subjects, n_trials, n_voxels) to (n_subjects, n_runs, n_conditions, n_voxels) to comply with the decompose_pattern_into_group_indiv_noise function."""
    # Extract each run and concatenate in third dimension
    data_reshaped = np.zeros((data.shape[0], info[part_column].max(), info[cond_column].max(), data.shape[-1]))
    for i in range(1, info.run.max()+1):
        run_data = data[:, info.run == i, :]
        data_reshaped[:, i-1, :, :] = run_data

    # Set nans to number and print number of nans
    print(f'Setting {np.sum(np.isnan(data_reshaped))} nan values to zero.')
    data_reshaped = np.nan_to_num(data_reshaped)

    # Mean centering
    if mean_centering:
        # Subtract the mean across subjects and runs from each voxel
        mean_across_conditions = np.mean(data_reshaped, axis=(2))
        # Repeat the mean across conditions for each condition
        mean_across_conditions = np.repeat(mean_across_conditions[:, :, np.newaxis, :], data_reshaped.shape[2], axis=2)
        data_reshaped = data_reshaped - mean_across_conditions

    return data_reshaped

if __name__ == "__main__":
    T = pd.read_csv(
        mdtb_dir + '/participants.tsv', delimiter='\t')
    subject_subset = T.participant_id[T['ses-rest'] == 1].tolist()
    mdtb_dataset = DataSetMDTB(mdtb_dir)
    mdtb_dataset.extract_all(ses_id='ses-s1', type='TaskRun', atlas='MNISymC3', subj=subject_subset)
    mdtb_dataset.extract_all(ses_id='ses-s2', type='TaskRun', atlas='MNISymC3', subj=subject_subset)
    # mdtb_dataset.extract_all(ses_id='ses-s1', type='CondRun', atlas='MNISymDentate1', smooth=None)
    # mdtb_dataset.extract_all(ses_id='ses-s2', type='CondRun', atlas='MNISymDentate1', smooth=None)

    # smooth_mdtb_fs32k(ses_id='ses-s1', type='CondHalf', smooth=4, kernel='fwhm')
    # smooth_mdtb_fs32k(ses_id='ses-s2', type='CondHalf', smooth=4, kernel='fwhm')
    # -- Extract resting-state timeseries --
    # for s in [2,4,5,6]:
    #     mdtb_dataset.extract_all(ses_id='ses-s1', type='CondHalf', atlas='MNISymC3', smooth=s)
    #     mdtb_dataset.extract_all(ses_id='ses-s2', type='CondHalf', atlas='MNISymC3', smooth=s)

    # -- Extract task timeseries --
    # mdtb_dataset.extract_all(ses_id='ses-s1', type='FixTseries', atlas='MNISymC3', smooth=2.0, subj=[23])
    # mdtb_dataset.extract_all(ses_id='ses-s1', type='FixTseries', atlas='fs32k', smooth=2.0, subj=[23])
    # mdtb_dataset.extract_all(ses_id='ses-s1', type='FixTseries', atlas='fs32k', smooth=2.0, subj=[16])
    # mdtb_dataset.extract_all(ses_id='ses-s2', type='FixTseries', atlas='MNISymC3', smooth=2.0)
    # mdtb_dataset.extract_all(ses_id='ses-s2', type='FixTseries', atlas='fs32k', smooth=2.0, subj=[23])

    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06FixRun', space='MNISymC3', ses_id='ses-s1', subj=['sub-31'])
    # conn.get_connectivity_fingerprint(dname,
    #                                 type='Fus06FixRun', space='MNISymC3', ses_id='ses-s2')

    # # -- Get connectivity fingerprint for rest session --
    # T = pd.read_csv(
    #     mdtb_dir + '/participants.tsv', delimiter='\t')
    # subject_subset = T.participant_id[T['ses-rest'] == 1].tolist()
    # subject_indices = T.participant_id[T['ses-rest'] == 1].index.tolist() # get indices of subjects
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest', subj=subject_subset)
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Net67Run', space='SUIT3', ses_id='ses-rest', subj=subject_subset)
    # # conn.get_connectivity_fingerprint(dname,
    # #     type='Net69Run', space='MNISymC2', ses_id='ses-rest')
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Net300Run', space='MNISymC2', ses_id='ses-rest', subj=subject_subset)
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Net300Run', space='SUIT3', ses_id='ses-rest', subj=subject_subset)
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Ico42Run', space='MNISymC2', ses_id='ses-rest', subj=subject_subset)
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Ico42Run', space='SUIT3', ses_id='ses-rest', subj=subject_subset)

    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Ico162Run', space='MNISymC2', ses_id='ses-rest', subj=subject_subset)
    # # conn.get_connectivity_fingerprint(dname,
    # #                                   type='Ico162Run', space='SUIT3', ses_id='ses-rest', subj=subject_subset)


    # Exctract Rest timeseries & connectivity fingerprint
    # mdtb_dataset = DataSetMDTB(mdtb_dir)   
    # Need to first rerun GLM to get residuals using script mdtb_rerun_glm.py
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06ResRun', space='MNISymC3', ses_id='ses-s1', subj=subject_subset)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Ico642Run', space='fs32k', ses_id='ses-rest', subj=None)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06All', space='MNISymC3', ses_id='ses-s2', subj=subject_subset)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06FixRun', space='fs32k', ses_id='ses-rest', subj=subject_subset)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06FixRun', space='fs32k', ses_id='ses-s1', subj=subject_subset)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06FixRun', space='fs32k', ses_id='ses-s2', subj=subject_subset)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for importing the HCP data set to general format.
Created on 4/25/2022 at 12:18 PM
Author: dzhi
"""
import pandas as pd
import shutil, os, sys, time, subprocess
from pathlib import Path
from copy import deepcopy
import Functional_Fusion.dataset as ds
import numpy as np
from nibabel import cifti2
import numpy as np
from neuromaps import transforms
import nibabel as nb
import Functional_Fusion.atlas_map as am
import Functional_Fusion.util as ut
import mat73
import scipy.io as spio
import nitools as nt
import FusionModel.util as fut

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/data/tge/Tian/UKBB_full/imaging'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

hcp_dir = '/data/tge/Tian/HCP_img'
atlas_dir = base_dir + '/Atlases'
hem_name = ['cortex_left', 'cortex_right']


def import_FIX_extended(source_dir, dest_dir, participant_id):
    """Imports the HCP preprocessed resting state files and
       supporting files for ICA-FIX denoise (FIX_extended)
       into HCP unrelated 100 subject data folder
    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        participant_id (str): ID of participant
    """
    run_name = ['REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL']
    # bar = progressbar.ProgressBar(maxval=20, widgets=[progressbar.Bar('=', '[', ']'), ' ',
    #                                        progressbar.Percentage()])
    for run, run_n in enumerate(run_name):
        # move data into the corresponding session folder

        try:
            print(f"   copying folder /rfMRI_{run_n}...")
            start = time.perf_counter()
            shutil.copytree(source_dir + '/MNINonLinear/Results' + f'/rfMRI_{run_n}',
                            dest_dir + f'/rfMRI_{run_n}', dirs_exist_ok=True)
            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - duration {elapse}.")
        except:
            print('skipping ' + f'/rfMRI_{run_n}')


def check_timepoints(networks):
    dest_dir = networks.split('signal')[0]
    dest_dir + f'/Net69_space-fs32k.dscalar.nii'
    for p, participant_id in enumerate(T.participant_id):
        print(p, participant_id)
        for run in range(4):
            fname = Path(dataset.func_dir.format(
                participant_id)) / f'sub-{participant_id}_run-{run}_space-MSMSulc.dtseries.nii'
            img = cifti2.load(str(fname))
            data = img.get_fdata()
            if data.shape[0] != 1200:
                print(f'Wrong timepoints for {fname}')
                print(data.shape)


def ica_networks_vol2surf(networks):
    fslr = transforms.mni152_to_fslr(networks, '32k')
    lh, rh = fslr
    # Save object as cifti
    structure = ['CORTEX_LEFT', 'CORTEX_RIGHT']
    seed_names = ['Network_{}'.format(i)
                  for i in range(1, len(lh.agg_data()) + 1)]
    bpa = nb.cifti2.ScalarAxis(seed_names)
    # lh = cifti2.Cifti2Image(lh, transforms.get_cifti2_axes('32k'))
    print(f'Writing {networks} ...')

    # Remove medial wall
    lh_masked = [data[atlas.mask[0]] for data in lh.agg_data()]
    rh_masked = [data[atlas.mask[1]] for data in rh.agg_data()]

    # --- Build a connectivity CIFTI-file and save ---
    # Make the atlas object
    atlas, atlas_info = am.get_atlas('fs32k', ut.atlas_dir)
    bmc = atlas.get_brain_model_axis()

    header = nb.Cifti2Header.from_axes((bpa, bmc))
    cifti_img = nb.Cifti2Image(
        dataobj=np.c_[lh_masked, rh_masked], header=header)
    dest_dir = networks.split('signal')[0]
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    nb.save(cifti_img, dest_dir + f'/Net69_space-fs32k.dscalar.nii')


def check_vertices(networks):
    net = nb.load(networks)
    # Check number of vertices
    for n in net.header.get_axis(1).iter_structures():
        print(f'{n[0]}: {n[1]}')


# net.header.get_axis(1)
# ['CIFTI_STRUCTURE_CORTEX_LEFT']

def scan_folders_and_save_subj_id(dir):
    # Get a list of all items in the current directory
    items = os.listdir(dir)
    # Create a pandas DataFrame with the folder names
    df = pd.DataFrame(items, columns=['participant_id'])
    df.to_csv(dir + '/participants.tsv', index=False, sep='\t')

def random_sampling_HCP_unrelated_family_subjects(original_df, num_subj=40):
    df = original_df.loc[original_df['full_4_runs']==True]

    pool = np.arange(1, 458)
    np.random.shuffle(pool)
    subj_list = pd.DataFrame()
    i = 0
    while(len(subj_list) < num_subj):
        this_group = pool[i]
        this_subs = df.loc[df['group_id'] == this_group]
        current_num = len(subj_list)
        
        if current_num + len(this_subs) <= num_subj:
            subj_list = pd.concat([subj_list, this_subs])

        i += 1

    return subj_list.sort_index()
        

def import_hcp_timeseries(data_dir, derivative_dir, ses_id='ses-rest1'):
    myatlas, _ = am.get_atlas('fs32k')
    if ses_id == 'ses-rest1':
        run_name = ['REST1_LR', 'REST1_RL']
    elif ses_id == 'ses-rest2':
        run_name = ['REST2_LR', 'REST2_RL']
    else:
        raise NameError('Unknown session id!')

    T = pd.read_csv(data_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        dest_dir = f'{derivative_dir}/{s}/data'
        # if os.path.exists(dest_dir + f'/{s}_space-fs32k_{ses_id}_Tseries.dscalar.nii'):
        #     print(f"Already imported subject {s} {ses_id} Tseries, skipping...")
        # else:
        print(f"-- Start importing subject {s} {ses_id} Tseries --")
        try:
            data = []
            info = pd.DataFrame()
            start = time.perf_counter()
            for run, run_n in enumerate(run_name):
                left_data = nb.load(f'{data_dir}/{s}' + 
                                    f'/rfMRI_{run_n}_Atlas_MSMAll_hp2000_clean_lh.func.gii')
                right_data = nb.load(f'{data_dir}/{s}' + 
                                    f'/rfMRI_{run_n}_Atlas_MSMAll_hp2000_clean_rh.func.gii')

                left_ts = [x.data for x in left_data.darrays]
                left_ts = np.reshape(left_ts, (len(left_ts), len(left_ts[0])))
                left_ts = left_ts[:, myatlas.mask[0]]

                right_ts = [x.data for x in right_data.darrays]
                right_ts = np.reshape(right_ts, (len(right_ts), len(right_ts[0])))
                right_ts = right_ts[:, myatlas.mask[1]]

                this_run = np.hstack([left_ts, right_ts])
                data.append(this_run)

                this_info = pd.DataFrame({'sn': [s] * this_run.shape[0],
                    'run': [run+1] * this_run.shape[0]})
                info = pd.concat([info, this_info], ignore_index=True)

            data = np.vstack(data)
            num_tpoint = [f'T{i+1:04}' for i in range(data.shape[0])]
            info['timepoint'] = num_tpoint
            info['task'] = 'rest'
            info['time_id'] = [i+1 for i in range(data.shape[0])]
            info['names'] = num_tpoint

            C = myatlas.data_to_cifti(data, num_tpoint)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            
            nb.save(C, dest_dir +
                        f'/{s}_space-fs32k_{ses_id}_Tseries.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_{ses_id}_info-Tseries.tsv', 
                        sep='\t', index=False)

            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - time {elapse}")
        except:
            print(f"   Incomplete - {s} {ses_id} Tseries!")


def step1_make_yeo_hcp_sess_tseries(data_dir, res_dir):
    run_name = ['REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL']

    T = pd.read_csv(data_dir + '/HCP80_training+validation_set.tsv', delimiter='\t')
    for s in T.participant_id:
        dest_dir = f'{data_dir}/{s}'

        for i, run_n in enumerate(run_name):
            fmri_list_dir = res_dir + '/profile_list/training_set' + \
                            f'/sub{s}_sess{i+1}.txt'
            # if (os.path.exists(dest_dir + f'/{s}_run{i}.dtseries.nii')) and \
            #     (not os.path.exists(fmri_list_dir)):
            #     # Writing the fMRI list file
            #     with open(fmri_list_dir, 'w') as file:
            #         file.write(f'{dest_dir}/{s}_run{i}.dtseries.nii\n')
            #     print(f'Done writing fMRI list file sub {s} run {i}')
            if os.path.exists(dest_dir + f'/{s}_run{i}.dtseries.nii'):
                print(f"Already exist subject {s} run {i} Tseries, write in fMRI_list .txt")
                # Writing the fMRI list file
                with open(fmri_list_dir, 'w') as file:
                    file.write(f'{dest_dir}/{s}_run{i}.dtseries.nii\n')
                print(f'Done writing fMRI list file sub {s} run {i}')
            else:
                try:
                    start = time.perf_counter()

                    left_data = f'{dest_dir}/rfMRI_{run_n}_Atlas_MSMAll_hp2000_clean_lh.func.gii'
                    right_data = f'{dest_dir}/rfMRI_{run_n}_Atlas_MSMAll_hp2000_clean_rh.func.gii'

                    cmd = (f'wb_command -cifti-create-dense-timeseries '
                           f'{dest_dir}/{s}_run{i}.dtseries.nii '
                           f'-left-metric {left_data} -right-metric {right_data} '
                           f'-timestep 0.72 -timestart 0')
                    subprocess.run(cmd, shell=True, check=True)

                    finish = time.perf_counter()
                    elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
                    print(f"Subject {s} run {i} Tseries. Done - time {elapse}")
                except:
                    print(f"Incomplete - {s} run {i} Tseries!")


def step1_smooth_yeo_hcp_sess_tseries(bulk, data_dir, res_dir, smooth=4, kernel=None):
    run_name = ['REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL']
    surf_L = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_midthickness.surf.gii'
    surf_R = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_midthickness.surf.gii'

    T = pd.read_csv(bulk, delimiter='\t')
    for s in T.participant_id:
        dest_dir = f'{data_dir}/{s}'

        for i, run_n in enumerate(run_name):
            fmri_list_dir = res_dir + '/profile_list/training_set' + \
                            f'/sub{s}_sess{i+1}.txt'

            unsmoothed_ts = dest_dir + f'/{s}_run{i}.dtseries.nii'
            cifti_out = dest_dir + f'/{s}_run{i}_desc-sm{smooth}' \
                    f'{kernel if kernel is not None else ""}.dtseries.nii'
            
            if os.path.exists(cifti_out):
                print(f"Already exist subject {s} run {i} smoothed Tseries")
                # # Writing the fMRI list file
                # with open(fmri_list_dir, 'w') as file:
                #     file.write(f'{dest_dir}/{s}_run{i}.dtseries.nii\n')
                # print(f'Done writing fMRI list file sub {s} run {i}')
            else:
                print(f'- Smoothing data for {s} fs32k in {smooth}mm ...')
                
                start = time.perf_counter()
                contain_nan = False
                # Load the unsmoothed data
                C = nb.load(unsmoothed_ts)

                # fill nan with zeros if unsmoothed data contains any
                if np.isnan(C.get_fdata()).any():
                    contain_nan = True
                    mask = np.isnan(C.get_fdata())
                    C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
                    nb.save(C, f'{s}_tmp.dscalar.nii')
                    unsmoothed_ts = f'{s}_tmp.dscalar.nii'

                # Write in smoothed surface data (filled with 0)
                smooth_cmd = f"wb_command -cifti-smoothing {unsmoothed_ts} " \
                            f"{smooth} {smooth} COLUMN {cifti_out} " \
                            f"{f'-{kernel} ' if kernel is not None else ''}" \
                            f"-left-surface {surf_L} -right-surface {surf_R} " \
                            f"-fix-zeros-surface"
                subprocess.run(smooth_cmd, shell=True)
                
                if contain_nan:
                    os.remove(f'{s}_tmp.dscalar.nii')
                    # Replace 0s back to NaN (we don't want the 0s impact model learning)
                    C = nb.load(cifti_out)
                    data = C.get_fdata()
                    data[mask] = np.nan
                    C = nb.Cifti2Image(dataobj=data, header=C.header)
                    nb.save(C, cifti_out)

                finish = time.perf_counter()
                elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
                print(f"Done subject {s} run {i}- time {elapse}.")
            

def make_yeo_hcp_sess_rsfc(subj_list, data_dir, ses_id='ses-rest1'):
    myatlas, _ = am.get_atlas('fs32k')
    if ses_id == 'ses-rest1':
        run_name = ['sess1', 'sess2']
    elif ses_id == 'ses-rest2':
        run_name = ['sess3', 'sess4']
    
    T = pd.read_csv(subj_list, delimiter='\t')
    for s in T.participant_id:
        start = time.perf_counter()
        dest_dir = f'{data_dir}/sub{s}'
        data, col_name = [], []
        info = pd.DataFrame()
        for run, run_n in enumerate(run_name):
            d_file = f'{dest_dir}/{run_n}' + \
            f'/sub{s}_{run_n}_fs_LR_32k_roifs_LR_900.surf2surf_profile.mat'    
            this_run = mat73.loadmat(d_file)['profile_mat'].T.astype(np.int8)
            data.append(this_run)
            this_names = [f'ROI_{i+1:04}' for i in range(this_run.shape[0])]
            col_name.append(this_names)
            this_info = pd.DataFrame({'sn': [s] * this_run.shape[0],
                                      'sess': [ses_id] * this_run.shape[0],
                                      'run': [run+1] * this_run.shape[0],
                                      'half': [run+1] * this_run.shape[0],
                                      'net_id': np.arange(this_run.shape[0])+1,
                                      'names': this_names})
            info = pd.concat([info, this_info], ignore_index=True)

        data = np.vstack(data)[:, np.concatenate(myatlas.mask)]
        C = myatlas.data_to_cifti(data, sum(col_name, []))

        out_dir = hcp_dir + f'/derivatives/{s}/data'
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        
        nb.save(C, out_dir +
                    f'/{s}_space-fs32k_{ses_id}_ROI1483Run_desc-sm4fwhm_binarized.dscalar.nii')
        info.to_csv(out_dir + f'/{s}_{ses_id}_info-ROI1483Run.tsv', 
                    sep='\t', index=False)
        
        print(f'Write HCP rsfc (ROI1483) for subject {s} {ses_id}!')
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"Done - time {elapse}")

def make_yeo_hcp_custom_rsfc(subj_list, out_dir, type='ICA15Run',
                             ses_id='ses-rest1', smooth=None):
    """ Replication of MS-HBM algorithm 
        Step 1: prepare the input rsfc (training set + validation)
        This function is to generate the resting-state functional 
        connectivity profile for Kong2019 MS-HBM algorithm input only. 
    """
    myatlas, _ = am.get_atlas('fs32k')
    mask = np.concatenate(myatlas.mask)
    hcp_dataset = ds.DataSetHcpResting(hcp_dir)
    if ses_id == 'ses-rest1':
        run_name = ['sess1', 'sess2']
    elif ses_id == 'ses-rest2':
        run_name = ['sess3', 'sess4']
    else:
        raise ValueError('Unknown session ID.')
    
    T = pd.read_csv(subj_list, delimiter='\t')
    for i, s in enumerate(T.participant_id):
        start = time.perf_counter()
        data_dir = hcp_dataset.data_dir.format(s)
        data, col_name = [], []
        
        # basic info / data file names
        info_raw = pd.read_csv(data_dir
                            + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t')   
        if smooth is None or (smooth == 0):
            file_name = f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        else:
            file_name = f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        
        try:
            dat = nb.load(data_dir + file_name)
            # this_data.append(atlas.read_data(data_dir.format(s) + file_name).T)
            dat = dat.get_fdata()
            data = np.zeros((dat.shape[0], mask.size))
            data[:, mask] = dat
        except FileNotFoundError as e:
            print(f'Missing data {s} fs32k {ses_id} {type}')
            continue

        for run, run_n in enumerate(run_name):
            d_folder = f'{out_dir}/sub{i+1}/{run_n}'
            d_file = f'/sub{i+1}_{run_n}_fs_LR_32k_roifs_LR_900.surf2surf_profile.mat'
            os.makedirs(os.path.dirname(d_folder+d_file), exist_ok=True)

            ind = info_raw["run"] == run+1
            out_data = {'profile_mat': data[ind,:].T}
            spio.savemat(d_folder + d_file, out_data)

        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"Done - fs32k subject {s} all runs in {ses_id} {type}, time {elapse}")


def step1_generate_fmri_list():
    fmri_dir = '/data/tge/Tian/HCP_img/rfMRI/fix_32k'
    res_dir = '/data/tge/dzhi/workspace/res/data_list/fMRI_list'
    T = pd.read_csv(fmri_dir + '/HCP923_test_set.tsv', delimiter='\t')

    # Generate 40 lines
    for sub_id in T.participant_id:
        for ses in [1,2,3,4]:
            # Open the file in write mode
            with open(res_dir + f'/sub{sub_id}_sess{ses}.txt', 'w') as file:
                # Write each line with the format "subject<index>"
                # file.write(f'{res_dir}/sub{i}/sess4/sub{i}_sess4_fs_LR_32k_roifs_LR_900.surf2surf_profile.mat\n')
                # file.write(f'{fmri_dir}/{sub_id}/{sub_id}_run2.dtseries.nii\n'
                #         f'{fmri_dir}/{sub_id}/{sub_id}_run3.dtseries.nii\n')
                
                file.write(f'{fmri_dir}/{sub_id}/{sub_id}_run{ses-1}_desc-sm4fwhm.dtseries.nii\n')
    
    print(f"File '{res_dir}' generated with 40 lines.")

def step1_make_avrg_file():
    fmri_dir = '/data/tge/Tian/HCP_img/rfMRI/fix_32k'
    res_dir = '/data/tge/dzhi/workspace/res/profiles'
    T = pd.read_csv(fmri_dir + '/HCP80_training+validation_set.tsv', delimiter='\t')

    for i, s in enumerate(T.participant_id):
        for ses in [1,2,3,4]:
            shutil.copy(res_dir + f'/sub{s}/sess{ses}/sub{s}_sess{ses}_fs_LR_32k_roifs_LR_900.surf2surf_profile.mat', 
                        res_dir + f'/sub{i+1}/sess{ses}/sub{i+1}_sess{ses}_fs_LR_32k_roifs_LR_900.surf2surf_profile.mat')
    
    print('Done.')

def step3_generate_profile_list():
    profiles_dir = '/data/tge/dzhi/workspace/res/profiles'
    res_dir = '/data/tge/dzhi/workspace/res/profile_list/test_set'
    T = pd.read_csv('/data/tge/Tian/HCP_img/rfMRI/fix_32k/HCP923_test_set.tsv', delimiter='\t')

    # Generate 40 lines
    for ses in [1,2,3,4]:
        for sub_id in T.participant_id:
            # Open the file in write mode
            with open(res_dir + f'/sess{ses}.txt', 'a') as file:
                # Write each line with the format "subject<index>"
                file.write(f'{profiles_dir}/sub{sub_id}/sess{ses}/sub{sub_id}_sess{ses}_fs_LR_32k_roifs_LR_900.surf2surf_profile.mat\n')
    
    print(f"Done.")

def convert_yeo_prior_to_cifti(file, out_dir, col_names=None):
    atlas, _ = am.get_atlas('fs32k')
    left = spio.loadmat(file)['lh_labels'].reshape(-1)
    right = spio.loadmat(file)['rh_labels'].reshape(-1)
    left_labels = left[atlas.mask[0]]
    right_labels = right[atlas.mask[1]]
    parcel = np.concatenate([left_labels, right_labels])

    img = nt.make_label_cifti(parcel.reshape(-1,1), atlas.get_brain_model_axis(),
                              label_RGBA=col_names)
    
    nb.save(img, out_dir + f'/group.dlabel.nii')

def convert_MSHBM_prior_to_cifti(mat_file, color, col_names=None):
    import torch as pt
    import HierarchBayesParcel.evaluation as hev
    atlas, _ = am.get_atlas('fs32k_Asym')
    align = atlas.cifti_to_data('/data/tge/dzhi/Indiv_par/Kong_2019/group_prior' \
                       '/HCP_40/Kong-2019_MSHBM_HCP40_prob_prior.dscalar.nii')
    align = pt.tensor(align, dtype=pt.get_default_dtype())

    stem = os.path.dirname(mat_file)
    Prob = spio.loadmat(mat_file)['Params']['theta'][0][0]
    Prob = Prob[np.concatenate(atlas.mask),:].T
    Prob = pt.tensor(Prob, dtype=pt.get_default_dtype())

    indx = hev.matching_greedy(align, pt.softmax(Prob, dim=0))
    # Prob = Prob[indx,:]

    # C = atlas.data_to_cifti(Prob.cpu().numpy(), row_axis=col_names)
    # nb.save(C, f'{stem}/MSHBM_group_prior_HCP40training_k-17.dscalar.nii')

    parcel = Prob.argmax(axis=0) + 1
    img = nt.make_label_cifti(parcel.cpu().numpy(), atlas.get_brain_model_axis(),
                              column_names=['Kong_2019'],
                              label_names=['???'] + col_names,
                              label_RGBA=color)
    nb.save(img, f'{stem}/MSHBM_group_prior_HCP40training_k-17_hard_new_withoutcorrection.dlabel.nii')

if __name__ == "__main__":
    # if len(sys.argv) != 2:
    #     print("Usage: python import_hcp.py <i>")
    #     sys.exit(1)        

    # step1_generate_fmri_list()
    # step1_make_avrg_file()
    # step3_generate_profile_list()

    #### Making HCP training/validation/test set subject list
    # df = pd.read_csv('/data/tge/Tian/HCP_img/family_group_id.tsv', delimiter='\t')
    # training_set = random_sampling_HCP_unrelated_family_subjects(df, num_subj=40)
    # validation_set = random_sampling_HCP_unrelated_family_subjects(df.drop(index=training_set.index),
    #                                                                num_subj=40)
    # df = df.loc[df['full_4_runs']==True]
    # test_set = df.drop(index=np.concatenate([training_set.index, validation_set.index]))

    #### Convert and make Yeo prior to cifti file
    # convert_MSHBM_prior_to_cifti('/data/tge/dzhi/workspace/res/priors/Params_Final.mat', 
    #                              colors, col_names=network_names)
    # img = nt.make_label_cifti(parcel, atlas.get_brain_model_axis(),
    #                           column_names=['Kong_2019'])
    # nb.save(img, f'/test.dlabel.nii')

    # convert_yeo_prior_to_cifti('/data/tge/dzhi/workspace/res/group/group.mat',
    #                            '/data/tge/dzhi/workspace/res/group',
    #                            col_names='network_names')
    
    #### import HCP data
    # scan_folders_and_save_subj_id(hcp_dir + '/rfMRI/fix_32k')
    # import_hcp_timeseries(hcp_dir + '/rfMRI/fix_32k', hcp_dir + '/derivatives',
    #                       ses_id='ses-rest1')
    # import_hcp_timeseries(hcp_dir + '/rfMRI/fix_32k', hcp_dir + '/derivatives',
    #                       ses_id='ses-rest2')

    #### Relication of Kong 2019 MS-HBM algorithm
    # step1_make_yeo_hcp_sess_tseries(hcp_dir + '/rfMRI/fix_32k', 
    #                           '/data/tge/dzhi/workspace/res')
    step1_smooth_yeo_hcp_sess_tseries(hcp_dir + '/subj_list/test_split/' + 'HCP923_test_set_split_6.tsv',
                                       hcp_dir + '/rfMRI/fix_32k', 
                                      '/data/tge/dzhi/workspace/res', 
                                      smooth=4, kernel='fwhm')
    # make_yeo_hcp_sess_rsfc(hcp_dir + '/subj_list/HCP80_training+validation_set.tsv',
    #                        '/data/tge/dzhi/workspace/res/profiles',
    #                        ses_id='ses-rest1')
    # make_yeo_hcp_sess_rsfc(hcp_dir + '/subj_list/HCP80_training+validation_set.tsv',
    #                        '/data/tge/dzhi/workspace/res/profiles',
    #                        ses_id='ses-rest2')
    
    # make_yeo_hcp_custom_rsfc(hcp_dir + '/subj_list/HCP80_training+validation_set.tsv',
    #                        '/data/tge/dzhi/workspace/res/profiles',
    #                        type='Ico42Run', ses_id='ses-rest2', smooth='4fwhm')

    # T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    # for s in T.participant_id:
    #     print(f"-Start importing subject {s}")
    #     # old_id = s.replace('sub-','s',1)
    #     dir1 = os.path.join(orig_dir, str(s))
    #     dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
    #     import_func_resting(dir1, dir2, str(s))
    #     print(f"-Done subject {s}")

    # S3server_dir = os.path.join('X:/', 'HCP_1200')
    # to_dir = os.path.join(base_dir, 'HCP_UR100_rfMRI')
    # T = pd.read_csv(to_dir + '/participants.tsv', delimiter='\t')
    # for s in T.participant_id[1:20]:
    #     print(f"-Importing FIX_extended subject {s}")
    #     # old_id = s.replace('sub-','s',1)
    #     dir1 = os.path.join(S3server_dir, str(s))
    #     dir2 = os.path.join(to_dir, '%s/MNINonLinear/Results' % str(s))
    #     import_FIX_extended(dir1, dir2, str(s))
    #     print(f"-Done subject {s}")

    # Test dimensions of func data
    # check_timepoints()

    # Create reginfo file for all data
    # create_reginfo(log_message=True, ses_id='ses-rest1')
    # create_reginfo(log_message=True, ses_id='ses-rest2')

    # Get ICA Networks into surface space
    # ica_networks_vol2surf(networks=target_dir +
    #                       '/group_ica/dim_auto/signal/signal_components.nii.gz')

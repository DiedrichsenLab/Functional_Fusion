#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for extracting the UK Biobank data into general format.
Created on 3/8/2024 at 1:23 PM
Author: dzhi
"""
import shutil, os, sys, logging, time, subprocess
import pandas as pd
import numpy as np
import nibabel as nb
from pathlib import Path

import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as ut
# import IndividualParcellation.scripts.group_parcellation as gp

BASE_DIR = '/data/tge/Tian/UKBB_full/imaging'

atlas_dir = BASE_DIR + '/Atlases'
WORK_DIR = os.path.join(BASE_DIR, 'scripts')
IMG_DIR = os.path.join(BASE_DIR, 'rfMRI')
DERIVATIVES_DIR = os.path.join(BASE_DIR, 'derivatives')
hem_name = ['cortex_left', 'cortex_right']
SESS_DIC = {'ses-rest1': '20227_2_0', 'ses-rest2': '20227_3_0'}
ukb_dir = BASE_DIR

# Configure the logging settingsd
logging.basicConfig(filename=WORK_DIR + '/log/missing_subj_list.log',
                    level=logging.ERROR, format='%(message)s')

def is_folder_exist_and_empty(folder_path):
    ''' Check if a folder exists and empty
    Args:
        folder_path: give the folder directory
    Returns:
        whether a given folder is exist and empty
    '''
    if not os.path.exists(folder_path):
        return False

    return not any(os.scandir(folder_path))

def is_folder_exist_and_not_empty(folder_path):
    ''' Check if a folder exists and empty
    Args:
        folder_path: give the folder directory
    Returns:
        whether a given folder is exist and empty
    '''
    return os.path.exists(folder_path) and \
        any(os.scandir(folder_path))

def check_within_subj_reliability(subj_list="test.tsv", sess=['ses-rest1'], type=['ICA25All'],
                                  smooth=None):
    print(f'Start loading data: UKB - {sess} - {type} ...')
    tic = time.perf_counter()
    data, cond_vec, part_vec, subj_ind = gp.build_ukb_datasets(BASE_DIR, subj_list,
                                                        space='MNIAsymC2', ses_list=sess,
                                                        type=type, smooth=smooth)
    toc = time.perf_counter()
    print(f'Done loading. Used {toc - tic:0.4f} seconds!')


def extract_ukb_timeseries(subj, ses_id='ses-rest1', type='Tseries', 
                           atlas='MNISymC3', smooth=None):
    df = pd.DataFrame({'participant_id': [subj]})
    sub_file = WORK_DIR + f'/{subj}.tsv'
    df.to_csv(sub_file, index=False, sep='\t')

    print(f'Step 3 -- Extracting {ses_id} {type} MNI space to cifti file')
    ukb_dataset = ds.DataSetUkbResting(ukb_dir, subj_id_file=f'/scripts/{subj}.tsv')
    ukb_dataset.extract_all(ses_id, type, atlas, smooth=smooth)

    os.remove(sub_file)
    print('Done!')


def make_subj_reginfo(s, type='Tseries', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file
    Args:
        type: file type
        ses_id: session id
    Returns:
        Write in the modified tsv file
    """
    ukb_dataset = ds.DataSetUkbResting(ukb_dir)
    ses_index = list(SESS_DIC.keys())

    # Make info
    print(f'Step 2 -- Making info file for subj {s} {ses_id} {type}')
    dest_dir = ukb_dataset.func_dir.format(s)
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    img = nb.load(dest_dir + 
                f'/{s}_run-{ses_index.index(ses_id)}_space-MNIAsym2.nii.gz')

    # num_tpoint = img.get_fdata().shape[3]
    num_tpoint = 490
    info = pd.DataFrame({'sn': [s] * num_tpoint,
                        'run': [ses_index.index(ses_id)+1] * num_tpoint,
                        'timepoint': [f'T{i+1:04}' for i in range(num_tpoint)],
                        'task': ['rest'] * num_tpoint,
                        'time_id': [i+1 for i in range(num_tpoint)]})

    info.to_csv(
        dest_dir + f'/{s}_{ses_id}_reginfo.tsv', sep='\t', index=False)
    print("Done!")

def make_subj_rsfcinfo(s, type='ICA25All', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file
    Args:
        type: file type
        ses_id: session id
    Returns:
        Write in the modified tsv file
    """
    ukb_dataset = ds.DataSetUkbResting(ukb_dir)
    ses_index = list(SESS_DIC.keys())

    # Make info
    print(f'Step 4 -- Making info file for subj {s} {ses_id} {type}')
    dest_dir = ukb_dataset.data_dir.format(s)
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    img = nb.load(dest_dir + 
                f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
    
    num_network = img.get_fdata().shape[0]
    info = pd.DataFrame({'sn': [s] * num_network,
                        'sess': [ses_id] * num_network,
                        'run': [ses_index.index(ses_id)+1] * num_network,
                        'half': [ses_index.index(ses_id)+1] * num_network,
                        'task': ['rest'] * num_network,
                        'net_id': np.arange(num_network)+1,
                        'names': [f'Network_{i+1}' for i in range(num_network)]})

    info.to_csv(
        dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)
    print("Done!")

def download_rfMRI_timeseries(subj_id, img_dir=IMG_DIR, sess='ses-rest1'):
    os.chdir(WORK_DIR)
    fetch_ukb = '/data/tge/dzhi/software/ukbfetch'
    ukb_key = 'k32568r674571.ukbkey'
    instance_id = SESS_DIC[sess]

    # Strip newline characters and other whitespaces
    dest_dir = os.path.join(img_dir, str(subj_id), instance_id)
    # Creat folder for current subject / instance
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    # If the derivatives already exist, skip downloading
    if check_file_completeness(subj_id, ses_id=sess):
        print(f'{subj_id} {sess} derivative files have completed!')
        return

    if not os.path.exists(dest_dir + '/rfMRI.ica/filtered_func_data_clean.nii.gz'):
        if is_folder_exist_and_empty(dest_dir):
            pass
            # tic = time.perf_counter()
            # # Download this subject's data if not exist
            # cmd = f'{fetch_ukb} -e{subj_id} -d{instance_id} -a{ukb_key}'
            # print(f"Step 1 -- Downloading {subj_id}, {instance_id} raw rfMRI...")
    
            # try:
            #     res = subprocess.run(cmd, check=True, shell=True)
            #     subprocess.run(f'mv {subj_id}_{instance_id}.zip {dest_dir}',
            #                     shell=True, check=True)
            # except subprocess.CalledProcessError as e:
            #     print("Error:", e)
            #     print("Command output:", e.stdout)
            #     print("Command error:", e.stderr)

            # toc = time.perf_counter()
            # min, sec = divmod(toc - tic, 60)
            # print(f'Step 1 -- Done - time used {int(min):02d}:{int(sec):02d}.')

            # try:
            #     print(f"Step 1 -- Unzip {subj_id}, {instance_id} raw rfMRI...")
            #     # Unzip and move to image folder
            #     cmd = (f'unzip {dest_dir}/{subj_id}_{instance_id}.zip -d {dest_dir}'
            #             f'&& rm {dest_dir}/{subj_id}_{instance_id}.zip')
            #     subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)

            #     # Extract only necessary files we will use
            #     subprocess.run(f'mv {dest_dir}/fMRI/rfMRI.ica '
            #                     f'{dest_dir}/fMRI/rfMRI_25.dr '
            #                     f'{dest_dir}/fMRI/rfMRI_100.dr {dest_dir} '
            #                     f'&& rm -r {dest_dir}/fMRI',
            #                     shell=True, check=True)
            # except:
            #     print(f'Step 1 -- Skipping extract {subj_id}_{instance_id}, the file incomplete!')
        
        elif os.path.exists(f'{dest_dir}/{subj_id}_{instance_id}.zip'):
            try:
                print(f"Step 1 -- Unzip {subj_id}, {instance_id} raw rfMRI...")
                # Unzip and move to image folder
                cmd = (f'unzip {dest_dir}/{subj_id}_{instance_id}.zip -d {dest_dir}'
                        f'&& rm {dest_dir}/{subj_id}_{instance_id}.zip')
                subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)

                # Extract only necessary files we will use
                subprocess.run(f'mv {dest_dir}/fMRI/rfMRI.ica '
                                f'{dest_dir}/fMRI/rfMRI_25.dr '
                                f'{dest_dir}/fMRI/rfMRI_100.dr {dest_dir} '
                                f'&& rm -r {dest_dir}/fMRI',
                                shell=True, check=True)
            except:
                print(f'Step 1 -- Skipping extract {subj_id}_{instance_id}, the file incomplete!')

    else:
        print(f'Step 1 -- {subj_id}_{instance_id} raw Tseries already exist!')


def extract_ukb_MNI_Tseries_nii(s, ses_id='ses-rest1', type='Tseries'):
    ses_index = list(SESS_DIC.keys())
    data_dir = BASE_DIR + f'/rfMRI/{s}/{SESS_DIC[ses_id]}/rfMRI.ica'
    out_dir = BASE_DIR + f'/derivatives/{s}/func'
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    if not os.path.exists(f'{out_dir}/{s}_run-{ses_index.index(ses_id)}_space-MNIAsym2.nii.gz'):
        try:
            print(f'Step 2 -- Wrapping {s} {ses_id} {type} to MNI space')
            
            # step 1: register fmri data from native space to MNI space
            cmd = (f'applywarp -i {data_dir}/filtered_func_data_clean.nii.gz '
                f'-r $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz '
                f'-o {out_dir}/{s}_run-{ses_index.index(ses_id)}_space-MNIAsym2.nii.gz '
                f'-w {data_dir}/reg/example_func2standard_warp.nii.gz')
            subprocess.run(cmd, shell=True, check=True, text=True)

            print(f"Step 2 -- Done warpping fMRI data to MNI space {s} {SESS_DIC[ses_id]}!")
            return True
        except:
            print(f'Step 2 -- Skipping warp raw data {s} {SESS_DIC[ses_id]} ...')
            return False
    else:
        print(f'Already extracted {type} {s} {s} {SESS_DIC[ses_id]}')
        return True


def get_connectivity(s, ses_id='ses-rest1', type='ICA25All', atlas='MNIAsymC2', smooth=None):
    my_atlas, _ = am.get_atlas(atlas)
    groupICA_dir = BASE_DIR + f'/derivatives/group'

    data_dir = BASE_DIR + f'/derivatives/{s}/data/'
    if smooth is None:
        data_file = f'{data_dir}/{s}_space-{atlas}_{ses_id}_Tseries.dscalar.nii'
        out_file = f'{data_dir}/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii'
    else:
        data_file = f'{data_dir}/{s}_space-{atlas}_{ses_id}_Tseries_desc-sm{smooth}.dscalar.nii'
        out_file = f'{data_dir}/{s}_space-{atlas}_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'

    if os.path.exists(data_file) and not os.path.exists(out_file):
        print(f'Step 4 -- Calculating RSFC for {s} {ses_id} {type}')

        # Load network time series
        if type == 'ICA25All':
            ica_dir = BASE_DIR + f'/rfMRI/{s}/{SESS_DIC[ses_id]}/rfMRI_25.dr'
            network_ts = np.loadtxt(ica_dir + '/dr_stage1.txt')
            good_comp = np.loadtxt(groupICA_dir + '/rfMRI_GoodComponents_d25_v1.txt',
                                    dtype=int) - 1
            network_ts = network_ts[:, good_comp]
        elif type == 'ICA100All':
            ica_dir = BASE_DIR + f'/rfMRI/{s}/{SESS_DIC[ses_id]}/rfMRI_100.dr'
            network_ts = np.loadtxt(ica_dir + '/dr_stage1.txt')
            good_comp = np.loadtxt(groupICA_dir + '/rfMRI_GoodComponents_d100_v1.txt',
                                    dtype=int) - 1
            network_ts = network_ts[:, good_comp]
        else:
            raise ValueError('Currently only support ICA components.')

        # Load raw cerebellum time series
        cereb_ts = nb.load(data_file).get_fdata()
        # calculate rsfc: cereb_ts (T, P), network_ts (T, D)
        rsfc = ut.correlate(cereb_ts, network_ts)
        # Save the data
        C = my_atlas.data_to_cifti(rsfc, [f"network_{r:02}" for r in range(rsfc.shape[0])])
        nb.save(C, out_file)
    else:
        if os.path.exists(out_file):
            print(f'Step 4 -- RSFC has already been extracted {s} {ses_id} {type}')
        else:
            print(f'Step 4 -- Missing cifti file {s} {atlas} {ses_id} {type}, skipping...')


def smooth_ukb_fs32k(bulk, ses_id='ses-s1', type='Tseries', smooth=1):
    ukb_dataset = ds.DataSetUkbResting(ukb_dir)
    T = pd.read_csv(bulk, sep=' ', header=None, names=['participant_id', 'instance'])

    # get the surfaces for smoothing
    surf_L = ukb_dataset.atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_midthickness.surf.gii'
    surf_R = ukb_dataset.atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_midthickness.surf.gii'

    for s in T.participant_id.drop_duplicates():
        print(f'- Smoothing data for {s} fs32k {ses_id} in {smooth}mm ...')
        # Make smoothed file name
        dest_dir = ukb_dataset.data_dir.format(s)
        cifti_out = dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        contain_nan = False

        # Load the unsmoothed data
        input_file = ukb_dataset.data_dir.format(s) \
                     + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        C = nb.load(input_file)

        # fill nan with zeros if unsmoothed data contains any
        if np.isnan(C.get_fdata()).any():
            contain_nan = True
            mask = np.isnan(C.get_fdata())
            C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
            nb.save(C, f'{s}_tmp.dscalar.nii')
            input_file = f'{s}_tmp.dscalar.nii'

        # Write in smoothed surface data (filled with 0)
        smooth_cmd = f"wb_command -cifti-smoothing {input_file} " \
                     f"{smooth} {smooth} COLUMN {cifti_out} " \
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


def smooth_ukb_cerebellum(bulk, ses_id='ses-s1', type='Tseries', smooth=1):
    ukb_dataset = ds.DataSetUkbResting(ukb_dir)
    T = pd.read_csv(bulk, sep=' ', header=None, names=['participant_id', 'instance'])

    for s in T.participant_id.drop_duplicates():
        print(f'- Smoothing data for {s} MNIAsymC2 {ses_id} in {smooth}mm ...')
        # Make smoothed file name
        dest_dir = ukb_dataset.data_dir.format(s)
        cifti_out = dest_dir + f'/{s}_space-MNIAsymC2_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        contain_nan = False

        # Load the unsmoothed data
        input_file = ukb_dataset.data_dir.format(s) \
                     + f'/{s}_space-MNIAsymC2_{ses_id}_{type}.dscalar.nii'
        C = nb.load(input_file)
        
        # fill nan with zeros if unsmoothed data contains any
        if np.isnan(C.get_fdata()).any():
            contain_nan = True
            mask = np.isnan(C.get_fdata())
            C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
            nb.save(C, f'{s}_tmp.dscalar.nii')
            input_file = f'{s}_tmp.dscalar.nii'

        # Write in smoothed surface data (filled with 0)
        smooth_cmd = f"wb_command -cifti-smoothing {input_file} " \
                     f"{smooth} {smooth} COLUMN {cifti_out} " \
                     f"-fix-zeros-volume"
        subprocess.run(smooth_cmd, shell=True)
        
        if contain_nan:
            os.remove(f"{s}_tmp.dscalar.nii")
            # Replace 0s back to NaN (we don't want the 0s impact model learning)
            C = nb.load(cifti_out)
            data = C.get_fdata()
            data[mask] = np.nan
            C = nb.Cifti2Image(dataobj=data, header=C.header)
            nb.save(C, cifti_out)


def check_file_completeness(s, ses_id='ses-rest1'):
    data_dir = BASE_DIR + f'/derivatives/{s}/data/'
    file_paths = [data_dir + f'/{s}_space-MNIAsymC2_{ses_id}_Tseries.dscalar.nii',
                    data_dir + f'/{s}_space-MNIAsymC2_{ses_id}_ICA25All.dscalar.nii',
                    data_dir + f'/{s}_space-MNIAsymC2_{ses_id}_ICA100All.dscalar.nii',
                    data_dir + f'/{s}_space-fs32k_{ses_id}_Tseries.dscalar.nii',
                    data_dir + f'/{s}_space-fs32k_{ses_id}_ICA25All.dscalar.nii',
                    data_dir + f'/{s}_space-fs32k_{ses_id}_ICA100All.dscalar.nii']
                      
    missing_files = [file for file in file_paths if not os.path.exists(file)]
    return True if len(missing_files) == 0 else False


def check_completeness_move_ica(s, ses_id='ses-rest1'):
    instance_id = SESS_DIC[ses_id]

    if not check_file_completeness(s, ses_id=ses_id):
        print(f'Step 5 -- {s} {ses_id} file incomplete, skipping...')
        logging.error(f'{s} {instance_id}')
    else:
        raw_dir = os.path.join(IMG_DIR, str(s), instance_id)
        dest_dir = BASE_DIR + f'/derivatives/{s}/ica/'
        
        # Delete func folder
        func_dir = BASE_DIR + f'/derivatives/{s}/func/'
        if is_folder_exist_and_not_empty(func_dir):
            subprocess.run(f'rm -r {func_dir}', shell=True, check=True)
        
        # Save ICA files into derivatives folder
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        if not os.path.exists(dest_dir+'/rfMRI_25.dr') or \
            (not os.path.exists(dest_dir+'/rfMRI_100.dr')):

            try:
                subprocess.run(f'mv {raw_dir}/rfMRI_25.dr '
                            f'{raw_dir}/rfMRI_100.dr {dest_dir} '
                            f'&& rm -r {raw_dir}',
                            shell=True, check=True)
                print(f'Step 5 -- {s} {ses_id} ICA files copied to derivatives!')
            except:
                print(f'Step 5 -- No ICA files for {s} {ses_id}.')
        
        else:
            print(f'Step 5 -- {s} {ses_id} completed!')


def all_in_one(bulk, ses='ses-rest1', d_type='Tseries'):
    """This is the wrapper function of prepare the UKB resting-state
       fMRI data, including downloading, wrapping to MNI group space,
       transferring to cifti format, extracting rsFC, and check file 
       completeness
    """
    T = pd.read_csv(bulk, sep=' ', header=None, 
                    names=['participant_id', 'instance_id'])
    for sub in T.participant_id.drop_duplicates():
        print(f'------------------------ {sub} {ses} ------------------------')
        # If the derivatives already exist, skip downloading
        if check_file_completeness(sub, ses_id=ses):
            print(f'{sub} {ses} data already prepared! Skipping Step 1-4...')
        
        else:
            # Step 1: download UKB data 
            download_rfMRI_timeseries(sub, img_dir=IMG_DIR, sess=ses)

            # Step 2: Warp the rs-fMRI from native space to MNI space
            has_mni_space = extract_ukb_MNI_Tseries_nii(sub, ses_id=ses, type=d_type)
            
            if has_mni_space:
                make_subj_reginfo(sub, type=d_type, ses_id=ses)
                
                # Step 3: extract raw time series (unsmoothed) in cifti format
                extract_ukb_timeseries(sub, ses_id=ses, type=d_type, atlas='MNIAsymC2')
                extract_ukb_timeseries(sub, ses_id=ses, type=d_type, atlas='fs32k')
                
                # Step 4: Calculate the rsfc for cortex and cerebellum
                get_connectivity(sub, ses_id=ses, type='ICA25All', atlas='MNIAsymC2')
                get_connectivity(sub, ses_id=ses, type='ICA100All', atlas='MNIAsymC2')
                get_connectivity(sub, ses_id=ses, type='ICA25All', atlas='fs32k',)
                get_connectivity(sub, ses_id=ses, type='ICA100All', atlas='fs32k')
                make_subj_rsfcinfo(sub, type='ICA25All', ses_id=ses)
                make_subj_rsfcinfo(sub, type='ICA100All', ses_id=ses)
            else:
                print(f'{sub} {ses} no MNI space nifit! Skipping Step 3-4...')

        # Step 5: Check file completeness and delete un-necessary files from server
        check_completeness_move_ica(sub, ses_id=ses)


if __name__ == "__main__":
    # Run script from command line
    # if len(sys.argv) != 4:
    #     print("Usage: python extract_ukb_rsfc.py <bulk_files> <ses-id> <type>")
    #     sys.exit(1)

    # all_in_one(BASE_DIR + '/scripts/subj_list/' + sys.argv[1],
    #             ses=sys.argv[2], d_type=sys.argv[3])

    ses = 'ses-rest1'
    d_type = 'Tseries'
    bulk_file = 'discovery/participants_discovery_split_196.bulk'    
    all_in_one(BASE_DIR + '/scripts/subj_list/' + bulk_file, ses=ses, d_type=d_type)

    #  -- fs32k smoothing (cortex)
    # smooth_ukb_fs32k(BASE_DIR + bulk_file, ses_id=ses, type='ICA25All', smooth=3)
    #  -- MNIAsymC2 smoothing (cerebellum)
    # smooth_ukb_cerebellum(BASE_DIR + bulk_file, ses_id=ses, type='ICA25All', smooth=1)

    # -- Get connectivity fingerprint --
    # get_connectivity(WORK_DIR + '/subj_list/' + sys.argv[1],
    #                  ses_id=sys.argv[2], type='ICA25All', atlas=sys.argv[3])
    # get_connectivity(WORK_DIR + '/subj_list/' + sys.argv[1],
    #                  ses_id=sys.argv[2], type='ICA100All', atlas=sys.argv[3])
    # for s in [None]:
    #     get_connectivity(BASE_DIR + '/test.tsv', ses_id='ses-rest2', 
    #                      type='ICA25All', atlas='MNIAsymC2', smooth=s)
    #     get_connectivity(BASE_DIR + '/test.tsv', ses_id='ses-rest2',
    #                      type='ICA100All', atlas='MNIAsymC2', smooth=s)

    #     get_connectivity(BASE_DIR + '/test.tsv', ses_id='ses-rest1', 
    #                      type='ICA25All', atlas='fs32k', smooth=s)



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script for extracting UKB resting state data

Created on 4/1/2024 at 2:07 PM
Author: dzhi
"""
import os, sys, re, logging, subprocess
import numpy as np
import nibabel as nb
import pandas as pd
import SUITPy as suit
import matplotlib.pyplot as plt

import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
from pathlib import Path
from nibabel.processing import resample_from_to, resample_to_output

BASE_DIR = '/data/tge/Tian/UKBB_full/imaging'
if not Path(BASE_DIR).exists():
    BASE_DIR = '/Volumes/diedrichsen_data$/data'
if not Path(BASE_DIR).exists():
    BASE_DIR = '/srv/diedrichsen/data'
if not Path(BASE_DIR).exists():
    BASE_DIR = 'Y:/data'
if not Path(BASE_DIR).exists():
    print('Cannot find data folder')

ukb_dir = BASE_DIR
atlas_dir = BASE_DIR + '/Atlases'
WORK_DIR = os.path.join(BASE_DIR, 'scripts')
IMG_DIR = os.path.join(BASE_DIR, 'rfMRI')
DERIVATIVES_DIR = os.path.join(BASE_DIR, 'derivatives')
hem_name = ['cortex_left', 'cortex_right']

# Configure the logging settingsd
logging.basicConfig(filename=WORK_DIR + '/log/extract_timeseries.log',
                    format='%(asctime)s - %(levelname)s - %(message)s')


def extract_ukb_timeseries(subj_id='/participants_filtered_final.tsv',
                           ses_id='ses-rest1', type='Tseries', 
                           atlas='MNISymC3', smooth=None):
    ukb_dataset = ds.DataSetUkbResting(ukb_dir, subj_id_file=subj_id)
    ukb_dataset.extract_all(ses_id, type, atlas, smooth=smooth)


def downsample_MNI_mask(res=3):
    """Downsamples a graymatter mask to a lower functional resolution
    Args:
        res: the target resolution
    Returns:
        None. Write out the resultant nifti file
    """
    adir = atlas_dir +'/tpl-MNI152NLin6AsymC'
    img_name = adir + '/tpl-MNI152NLin6AsymC_res-1_gmcmask.nii'
    out_name = adir + f'/tpl-MNI152NLin6AsymC_res-{res:d}_gmcmask.nii'
    in_img = nb.load(img_name)
    dimension = np.ceil(np.array(in_img.shape)/res).astype(int)
    if (dimension[0] % 2) == 0:
        dimension[0] = dimension[0]+1

    mag = np.eye(4)*res
    mag[3,3]=1
    affineM = in_img.affine @ mag
    temp_img = resample_from_to(in_img,(dimension,affineM))
    X=temp_img.get_fdata()
    X = (X+np.flip(X,axis=0))/2
    X=np.array(X>0.1,dtype=np.uint8)
    out_img = nb.Nifti1Image(X,temp_img.affine)
    nb.save(out_img,out_name);


def make_subj_reginfo(type='Tseries', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file
    Args:
        type: file type
        ses_id: session id
    Returns:
        Write in the modified tsv file
    """
    ukb_dataset = ds.DataSetUKBResting(ukb_dir)
    sess_dic = {'ses-rest1': '20227_2_0', 'ses-rest2': '20227_3_0'}
    ses_index = list(sess_dic.keys())

    T = pd.read_csv(ukb_dataset.base_dir + '/participants_filtered_final.tsv', sep='\t')
    for s in T.participant_id:
        # Make info
        print(f'Making info file for subj {s} {ses_id} {type}')
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


def make_subj_rsfcinfo(type='ICA25All', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file
    Args:
        type: file type
        ses_id: session id
    Returns:
        Write in the modified tsv file
    """
    ukb_dataset = ds.DataSetUkbResting(ukb_dir)
    sess_dic = {'ses-rest1': '20227_2_0', 'ses-rest2': '20227_3_0'}
    ses_index = list(sess_dic.keys())

    T = pd.read_csv(ukb_dataset.base_dir + '/test.tsv', sep='\t')
    for s in T.participant_id:
        # Make info
        print(f'Making info file for subj {s} {ses_id} {type}')
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



def make_participant_tsv(bulk, out_dir):
    # Convert UKB default .bulk file to standard tsv file
    df = pd.read_csv(bulk, sep=' ', header=None,
                     names=['participant_id', 'instance'])
    df_unique = df['participant_id'].drop_duplicates().to_frame()
    df_unique.to_csv(out_dir + '/participants.tsv',
                     sep='\t', index=False)


def extract_ukb_ts(bulk, ses_id='ses-rest1', type='Tseries', atlas='MNIAsymC2'):
    my_atlas, _ = am.get_atlas(atlas)
    raw_ts_name = 'filtered_func_data_clean.nii.gz'
    sess_dic = {'ses-rest1': '20227_2_0', 'ses-rest2': '20227_3_0'}

    # T = pd.read_csv(bulk, sep=' ', header=None, names=['participant_id', 'instance'])
    T = pd.read_csv(bulk, sep='\t')

    for s in T.participant_id.drop_duplicates():
        data_dir = BASE_DIR + f'/rfMRI/{s}/{sess_dic[ses_id]}/rfMRI.ica'
        out_dir = BASE_DIR + f'/derivatives/{s}/data'
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        if (os.path.exists(f'{data_dir}/{raw_ts_name}') and
                not os.path.exists(f'{out_dir}/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')):
            print(f'Extracting {s} {ses_id} {type}')

            # step 1: register fmri data from native space to MNI space
            cmd = (f'applywarp -i {data_dir}/filtered_func_data_clean.nii.gz '
                   f'-r $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz '
                   f'-o {out_dir}/tmp.nii.gz '
                   f'-w {data_dir}/reg/example_func2standard_warp.nii.gz')
            result = subprocess.run(cmd, shell=True, check=True, text=True)

            if result.returncode == 0:
                # step 2: extract cerebellum data from nifti to cifti
                data = my_atlas.read_data(f'{out_dir}/tmp.nii.gz', interpolation=0)
                C = my_atlas.data_to_cifti(data.T, [f"T{r:03}" for r in range(data.T.shape[0])])
                nb.save(C, out_dir +
                        f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
                result = subprocess.run(f'rm {out_dir}/tmp.nii.gz', shell=True,
                                        check=True, text=True)
            else:
                print(f"Command failed with return code {result.returncode}")
                logging.error(f"Failed to warp fMRI data from native to MNI space "
                              f"{s} {sess_dic[ses_id]} \n {result.returncode}")

        else:
            if os.path.exists(f'{out_dir}/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii'):
                print(f'Already extracted {type} {s} {s} {sess_dic[ses_id]}')
            else:
                print(f'Missing raw data {s} {sess_dic[ses_id]}, skipping...')
                logging.error(f"Missing fMRI data {s} {sess_dic[ses_id]}")


def smooth_ukb_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1):
    ukb_dataset = ds.DataSetUkbResting(ukb_dir)
    T = ukb_dataset.get_participants()

    # get the surfaces for smoothing
    surf_L = ukb_dataset.atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_midthickness.surf.gii'
    surf_R = ukb_dataset.atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_midthickness.surf.gii'

    for s in T.participant_id:
        print(f'- Smoothing data for {s} fs32k {ses_id} in {smooth}mm ...')
        # Load the unsmoothed data and fill nan with zeros
        C = nb.load(ukb_dataset.data_dir.format(s)
                    + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
        mask = np.isnan(C.get_fdata())
        C = nb.Cifti2Image(dataobj=np.nan_to_num(
            C.get_fdata()), header=C.header)
        nb.save(C, 'tmp.dscalar.nii')

        dest_dir = ukb_dataset.data_dir.format(s)
        cifti_out = dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        # Write in smoothed surface data (filled with 0)
        smooth_cmd = f"wb_command -cifti-smoothing tmp.dscalar.nii " \
                     f"{smooth} {smooth} COLUMN {cifti_out} " \
                     f"-left-surface {surf_L} -right-surface {surf_R} " \
                     f"-fix-zeros-surface"
        subprocess.run(smooth_cmd, shell=True)
        os.remove("tmp.dscalar.nii")

        # Replace 0s back to NaN (we don't want the 0s impact model learning)
        C = nb.load(cifti_out)
        data = C.get_fdata()
        data[mask] = np.nan
        C = nb.Cifti2Image(dataobj=data, header=C.header)
        nb.save(C, cifti_out)

        # # Copy info to the corresponding /smoothed folder
        # if not Path(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv').exists():
        #     info = pd.read_csv(ukb_dataset.data_dir.format(s)
        #                        + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t')
        #     info.to_csv(
        #         dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)
            

def extract_ukb_MNI_Tseries_nii(bulk, ses_id='ses-rest1', type='Tseries'):
    raw_ts_name = 'filtered_func_data_clean.nii.gz'
    sess_dic = {'ses-rest1': '20227_2_0', 'ses-rest2': '20227_3_0'}
    ses_index = list(sess_dic.keys())

    # T = pd.read_csv(bulk, sep=' ', header=None, names=['participant_id', 'instance'])
    T = pd.read_csv(bulk, sep='\t')

    for s in T.participant_id:
        data_dir = BASE_DIR + f'/rfMRI/{s}/{sess_dic[ses_id]}/rfMRI.ica'
        out_dir = BASE_DIR + f'/derivatives/{s}/func'
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        if not os.path.exists(f'{out_dir}/{s}_run-{ses_index.index(ses_id)}_space-MNIAsym2.nii.gz'):
            
            try:
                print(f'Wrapping {s} {ses_id} {type} to MNI space')
                
                # step 1: register fmri data from native space to MNI space
                cmd = (f'applywarp -i {data_dir}/filtered_func_data_clean.nii.gz '
                    f'-r $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz '
                    f'-o {out_dir}/{s}_run-{ses_index.index(ses_id)}_space-MNIAsym2.nii.gz '
                    f'-w {data_dir}/reg/example_func2standard_warp.nii.gz')
                subprocess.run(cmd, shell=True, check=True, text=True)

                print(f"Done warp fMRI data to MNI space {s} {sess_dic[ses_id]}!")
            except:
                print(f'Missing raw data {s} {sess_dic[ses_id]}, skipping...')
                logging.error(f"Missing fMRI data {s} {sess_dic[ses_id]}")
        else:
            print(f'Already extracted {type} {s} {s} {sess_dic[ses_id]}')


def extract_ukb_ts_cifti(bulk, ses_id='ses-rest1', type='Tseries', atlas='MNIAsymC2'):
    my_atlas, _ = am.get_atlas(atlas)
    sess_dic = {'ses-rest1': '20227_2_0', 'ses-rest2': '20227_3_0'}
    ses_index = list(sess_dic.keys())

    T = pd.read_csv(bulk, sep='\t')
    for s in T.participant_id:
        out_dir = BASE_DIR + f'/derivatives/{s}/data'
        func_dir = BASE_DIR + f'/derivatives/{s}/func'
        Path(out_dir).mkdir(parents=True, exist_ok=True)

        if not os.path.exists(f'{out_dir}/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii'):
            try:
                print(f'Extracting {s} {ses_id} {type}')
                ts_img = f'{func_dir}/{s}_run-{ses_index.index(ses_id)}_space-MNIAsym2.nii.gz'
                # step 2: extract cerebellum data from nifti to cifti
                data = my_atlas.read_data(ts_img, interpolation=0)
                C = my_atlas.data_to_cifti(data.T, [f"T{r:03}" for r in range(data.T.shape[0])])
                nb.save(C, out_dir +
                        f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            except:
                print(f'Missing raw time series .nii file {s} {sess_dic[ses_id]}.')
        else:
            print(f"Already extracted {s}_space-{atlas}_{ses_id}_{type}.dscalar.nii")


if __name__ == "__main__":
    # make_participant_tsv(BASE_DIR + '/ukb674571_rsfmri-return.bulk', BASE_DIR)
    # make_subj_reginfo(type='Tseries', ses_id='ses-rest1')
    # make_subj_reginfo(type='Tseries', ses_id='ses-rest2')

    #  -- fs32k smoothing
    # make_subj_rsfcinfo(type='ICA25All', ses_id='ses-rest1')
    # smooth_ukb_fs32k(ses_id='ses-rest1', type='ICA25All', smooth=3)

    # Run script from command line
    # if len(sys.argv) != 5:
    #     print("Usage: python extract_ukb_data.py <bulk files> <ses-id> <type> <space>")
    #     sys.exit(1)

    #  -- Extract timeseries --
    # extract_ukb_ts(WORK_DIR + '/subj_list/' + sys.argv[1],
    #                        ses_id=sys.argv[2], type=sys.argv[3], atlas=sys.argv[4])
    # if len(sys.argv) != 4:
    #     print("Usage: python extract_ukb_data.py <bulk files> <ses-id> <type>")
    #     sys.exit(1)

    # extract_ukb_MNI_Tseries_nii(WORK_DIR + '/subj_list/' + sys.argv[1],
    #                             ses_id=sys.argv[2], type=sys.argv[3])
    # extract_ukb_MNI_Tseries_nii(BASE_DIR + '/participants_filtered_final.tsv',
    #                              ses_id='ses-rest2', type='Tseries')
    


    # extract_ukb_ts(BASE_DIR + '/participants_filtered.tsv',
    #                        ses_id='ses-rest1', type='Tseries', atlas='MNISymC3')
    # extract_ukb_ts_cifti(BASE_DIR + '/test.tsv', ses_id='ses-rest2',
    #                      type='Tseries', atlas='MNISymC3')


    # downsample_MNI_mask(3)

    #  -- Extract smoothed timeseries in multiple resolutions --
    for s in [None]:
        extract_ukb_timeseries(subj_id='/test.tsv', ses_id='ses-rest1',
                               type='Tseries', atlas='fs32k', smooth=s)
        
    if len(sys.argv) != 4:
        print("Usage: python extract_ukb_data.py <bulk files> <ses-id> <atlas>")
        sys.exit(1)

    for s in [2,3,4,5,6]:
        extract_ukb_timeseries(subj_id='/scripts/subj_list/' + sys.argv[1],
                               ses_id=sys.argv[2], type='Tseries',
                               atlas=sys.argv[3], smooth=s)
    # extract_ukb_timeseries(ses_id='ses-rest2', type='Tseries', atlas='MNISymC2')
    # extract_ukb_timeseries(ses_id='ses-rest2', type='Tseries', atlas='MNISymC3')


    # -- Get connectivity fingerprint --

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script for extracting UKB resting state data

Created on 4/1/2024 at 2:07 PM
Author: dzhi
"""
import os, sys, re
import numpy as np
import nibabel as nb
import pandas as pd
import SUITPy as suit
import matplotlib.pyplot as plt

import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
from pathlib import Path
from nibabel.processing import resample_from_to, resample_to_output


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

hcp_dir = base_dir + '/UKB_rfMRI'
atlas_dir = base_dir + '/Atlases'

def extract_ukb_timeseries(ses_id='ses-rest1', type='Tseries', atlas='MNISymC3'):
    hcp_dataset = ds.DataSetHcpResting(hcp_dir)
    hcp_dataset.extract_all(ses_id, type, atlas)

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

def make_info(type='Tseries', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file

    Args:
        type: file type
        ses_id: session id

    Returns:
        Write in the modified tsv file
    """
    hcp_dataset = DataSetHcpResting(hcp_dir)

    T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    for p, participant_id in enumerate(T.participant_id):
        # Make info
        dest_dir = hcp_dataset.base_dir + \
            f'/derivatives/{participant_id}/data/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        info = pd.read_csv(dest_dir + f'{participant_id}_{ses_id}_info-{type}.tsv',
                           sep='\t')

        info['run_id'] = info['run'].copy()
        if ses_id == 'ses-rest2':
            info['run_id'] = info['run_id'] + 2

        info.to_csv(
            dest_dir + f'{participant_id}_{ses_id}_info-{type}.tsv', sep='\t', index=False)




if __name__ == "__main__":
    # downsample_MNI_mask(3)

    atlas, ainf = am.get_atlas('MNIAsymC3')
    # Read data from Nifti file using linear interpolation (1)
    X = atlas.read_data('/srv/diedrichsen/data/UKB_rfMRI/fMRI/rfMRI.ica/registered_ts.nii.gz', interpolation=0)

    # make_info(type='Tseries', ses_id='ses-rest1')
    # make_info(type='Tseries', ses_id='ses-rest2')
    #  -- Extract timeseries --
    # extract_ukb_timeseries(ses_id='ses-rest1', type='Tseries', atlas='MNISymC2')
    extract_ukb_timeseries(ses_id='ses-rest1', type='Tseries', atlas='MNISymC3')
    # extract_ukb_timeseries(ses_id='ses-rest2', type='Tseries', atlas='MNISymC2')
    extract_ukb_timeseries(ses_id='ses-rest2', type='Tseries', atlas='MNISymC3')


    # -- Get connectivity fingerprint --


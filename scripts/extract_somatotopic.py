# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73, subprocess
import numpy as np
import sys, os, time
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetSomatotopic
import Functional_Fusion.util as ut
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Somatotopic'
atlas_dir = base_dir + '/Atlases'

def extract_somatotopic(ses_id='ses-01',type='CondHalf',atlas='MNISymC3'):
    dataset = DataSetSomatotopic(data_dir)
    dataset.extract_all(ses_id,type,atlas, smooth=0)


def smooth_somatotopic_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1, kernel='gaussian'):
    dataset = DataSetSomatotopic(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Smoothing data for {s} fs32k {ses_id} in {smooth}mm {kernel} ...')

        start = time.perf_counter()
        file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        ut.smooth_fs32k_data(file, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


def mask_somatotopic_fs32k(ses_id='ses-s1', type='CondHalf', high_percent=0.1, low_percent=0.1,
                           smooth=None, z_transfer=False, binarized=False):
    myatlas, _ = am.get_atlas('fs32k')
    dataset = DataSetSomatotopic(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Mask data for {s} fs32k {ses_id} in high {high_percent} low {low_percent} ...')

        start = time.perf_counter()
        if smooth is not None:
            file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        else:
            file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'

        ut.mask_fs32k_data(file, high_percent=high_percent, low_percent=low_percent,
                           z_transfer=z_transfer, binarized=binarized)

        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


if __name__ == "__main__":
    smooth_somatotopic_fs32k(ses_id='ses-motor', type='CondHalf', smooth=4, kernel='fwhm')

    # --- Extracting Estimates ---
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')
    for s in [4,6,8,10]:
        print(f'Doing processing for {s}fwhm ...')
        mask_somatotopic_fs32k(ses_id='ses-motor', type=f'CondHalf', high_percent=0.1,
                         low_percent=0.1, smooth=f'{s}fwhm', z_transfer=True, binarized=False)

    # --- Group Average ---
    dataset = DataSetSomatotopic(data_dir)
    dataset.extract_all(type='CondAll', ses_id='ses-motor', atlas='MNISymC3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')

    
    # --- Show group average ---
    # dataset.plot_cerebellum(subject='group', savefig=True, colorbar=True)
    pass

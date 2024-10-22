# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys, os, time
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetSomatotopic
import Functional_Fusion.util as fut
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

    # get the surfaces for smoothing
    surf_L = dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.L.midthickness.surf.gii'
    surf_R = dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.R.midthickness.surf.gii'

    for s in T.participant_id:
        print(f'Smoothing data for {s} fs32k {ses_id} in {smooth}mm {kernel} ...')

        start = time.perf_counter()
        file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        fut.smooth_fs32k_data(file, surf_L, surf_R, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")

if __name__ == "__main__":
    smooth_somatotopic_fs32k(ses_id='ses-motor', type='CondHalf', smooth=4, kernel='fwhm')

    # --- Extracting Estimates ---
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')


    # --- Group Average ---
    dataset = DataSetSomatotopic(data_dir)
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')

    
    # --- Show group average ---
    dataset.plot_cerebellum(subject='group', savefig=True, colorbar=True)
    pass

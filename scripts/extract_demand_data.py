# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73, time
import numpy as np
import sys, subprocess
import Functional_Fusion.atlas_map as am
import Functional_Fusion.util as ut
from Functional_Fusion.dataset import DataSetDemand
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Demand'
atlas_dir = base_dir + '/Atlases'

def extract_demand(ses_id='ses-01', type='CondHalf', atlas='MNISymC3'):
    de_dataset = DataSetDemand(data_dir)
    de_dataset.extract_all(ses_id, type, atlas)


def smooth_demand_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1, kernel='gaussian'):
    dataset = DataSetDemand(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Smoothing data for {s} fs32k {ses_id} in {smooth}mm {kernel} ...')

        start = time.perf_counter()
        file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        ut.smooth_fs32k_data(file, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


def mask_demand_fs32k(ses_id='ses-s1', type='CondHalf', high_percent=0.1, low_percent=0.1,
                      smooth=None, z_transfer=False, binarized=False):
    myatlas, _ = am.get_atlas('fs32k')
    dataset = DataSetDemand(data_dir)
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
    # extract_demand(ses_id='ses-01', type='CondHalf', atlas='MNISymC2')
    de_dataset = DataSetDemand(data_dir)
    # de_dataset.group_average_data(ses_id='ses-01',type='CondHalf',atlas='MNISymC2')
    # de_dataset.group_average_data(ses_id='ses-01',type='CondHalf',atlas='MNISymC3')
    # de_dataset.group_average_data(ses_id='ses-01',type='CondHalf',atlas='fs32k')

    # for s in [2,3,5,7,9]:
    #     smooth_demand_fs32k(ses_id='ses-01', type='CondHalf', smooth=s, kernel='fwhm')

    for s in [4, 6, 8, 10]:
        print(f'Doing processing for {s}fwhm ...')
        mask_demand_fs32k(ses_id='ses-01', type=f'CondHalf', high_percent=0.1,
                          low_percent=0.1, smooth=f'{s}fwhm', z_transfer=True, binarized=False)

    T = de_dataset.get_participants()
    for s in T.participant_id:
        for ses in de_dataset.sessions:
            file = de_dataset.data_dir.format(s) + f'/{s}_{ses}_info-CondHalf.tsv'
            new_file = de_dataset.data_dir.format(s) + f'/{s}_{ses}_CondHalf.tsv'
            smooth_cmd = f"mv {file} {new_file}"
            subprocess.run(smooth_cmd, shell=True)

    de_dataset.plot_cerebellum(savefig=True, atlas='MNISymC3', colorbar=True)
    pass

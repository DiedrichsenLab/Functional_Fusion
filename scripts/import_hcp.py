#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for importing the HCP data set to general format.
Created on 4/25/2022 at 12:18 PM
Author: dzhi
"""
import pandas as pd
import shutil
from pathlib import Path
import os
import sys
import time
from copy import deepcopy
import dataset as ds
import numpy as np
from nibabel import cifti2
import numpy as np
from neuromaps import transforms
import nibabel as nb
import atlas_map as am
import util as ut

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'
if not Path(base_dir).exists():
    base_dir = 'Y:/data'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

orig_dir = os.path.join(base_dir, 'HCP_UR100_rfMRI')
target_dir = os.path.join(base_dir, 'FunctionalFusion/HCP')


def create_reginfo(log_message=False, ses_id='ses-rest1'):
    dataset = ds.DataSetHcpResting(str(target_dir))

    # Import general info
    info = pd.read_csv(
        target_dir + f'/{ses_id}_reginfo.tsv', sep='\t')

    T = dataset.get_participants()
    for _, id in T.iterrows():
        print(f'Creating reginfo for {id.participant_id}')

        # Ammend the reginfo.tsv file from the general file
        reginfo = deepcopy(info)
        reginfo.insert(loc=0, column='sn', value=[
            id.participant_id] * info.shape[0])

        # Make folder
        dest = target_dir + \
            f'/derivatives/{id.participant_id}/func/{id.participant_id}_{ses_id}_reginfo.tsv'
        Path(dest).parent.mkdir(parents=True, exist_ok=True)

        # Save reginfo.tsv file
        reginfo.to_csv(dest, sep='\t', index=False)


def import_func_resting(source_dir, dest_dir, participant_id):
    """Imports the HCP preprocessed resting state files
       into a BIDS/derivative structure
    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        participant_id (str): ID of participant
    """
    run_name = ['REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL']
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    for run, run_n in enumerate(run_name):

        # move data into the corresponding session folder
        src = (f'/rfMRI_{run_n}/rfMRI_{run_n}_Atlas_hp2000_clean.dtseries.nii')
        dest = (f'/sub-{participant_id}_run-{run}_space-MSMSulc.dtseries.nii')

        try:
            shutil.copyfile(source_dir + '/MNINonLinear/Results' + src,
                            dest_dir + dest)
        except:
            print('skipping ' + src)


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

if __name__ == "__main__":

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
    check_vertices(networks=target_dir +
                   '/group_ica/dim_auto/Net69_space-fs32k.dscalar.nii')

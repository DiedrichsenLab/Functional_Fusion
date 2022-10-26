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

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'
if not Path(base_dir).exists():
    base_dir = 'Y:/data'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

orig_dir = os.path.join(base_dir, 'HCP_UR100_rfMRI')
target_dir = os.path.join(base_dir, 'FunctionalFusion/HCP')

def import_func_resting(source_dir, dest_dir, participant_id):
    """Imports the HCP preprocessed resting state files
       into a BIDS/derivative structure
    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        participant_id (str): ID of participant
    """
    run_name=['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    for run,run_n in enumerate(run_name):

        # move data into the corresponding session folder
        src=(f'/rfMRI_{run_n}/rfMRI_{run_n}_Atlas_hp2000_clean.dtseries.nii')
        dest=(f'/sub-{participant_id}_run-{run}_space-MSMSulc.dtseries.nii')

        try:
            shutil.copyfile(source_dir+'/MNINonLinear/Results'+src, 
                    dest_dir+dest)
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
    run_name=['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    # bar = progressbar.ProgressBar(maxval=20, widgets=[progressbar.Bar('=', '[', ']'), ' ',
    #                                        progressbar.Percentage()])
    for run,run_n in enumerate(run_name):
        # move data into the corresponding session folder

        try:
            print(f"   copying folder /rfMRI_{run_n}...")
            start = time.perf_counter()
            shutil.copytree(source_dir+'/MNINonLinear/Results'+f'/rfMRI_{run_n}',
                    dest_dir+f'/rfMRI_{run_n}', dirs_exist_ok=True)
            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - duration {elapse}.")
        except:
            print('skipping ' + f'/rfMRI_{run_n}')


if __name__ == "__main__":
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        print(f"-Start importing subject {s}")
        # old_id = s.replace('sub-','s',1)
        dir1 = os.path.join(orig_dir, str(s))
        dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
        import_func_resting(dir1, dir2, str(s))
        print(f"-Done subject {s}")

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
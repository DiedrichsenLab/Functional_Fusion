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

base_dir = '/Volumes/diedrichsen_data$/data'
if sys.platform == "win32":  # for windows user
    base_dir = 'Y:/data'

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


if __name__ == "__main__":
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        print(f"-Start importing subject {s}")
        # old_id = s.replace('sub-','s',1)
        dir1 = os.path.join(orig_dir, str(s))
        dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
        import_func_resting(dir1, dir2, str(s))
        print(f"-Done subject {s}")

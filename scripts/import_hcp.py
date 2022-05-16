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
    for ss in [1, 2]:
        src = []
        dest = []
        # add a folder for session
        dest_sess_dir = dest_dir + f"/ses-{ss:02}"
        # Make the destination directory
        Path(dest_sess_dir).mkdir(parents=True, exist_ok=True)

        # move data into the corresponding session folder
        src.append(f'/rfMRI_REST{ss:01}_LR/rfMRI_REST{ss:01}_LR_Atlas_hp2000_clean.dtseries.nii')
        dest.append(f'/sub-{participant_id}_ses-{ss:02}_space-fsLR32k_run-01.dtseries.nii')

        src.append(f'/rfMRI_REST{ss:01}_RL/rfMRI_REST{ss:01}_RL_Atlas_hp2000_clean.dtseries.nii')
        dest.append(f'/sub-{participant_id}_ses-{ss:02}_space-fsLR32k_run-02.dtseries.nii')
        for i in range(len(src)):
            try:
                shutil.copyfile(source_dir+'/MNINonLinear/Results'+src[i], dest_sess_dir+dest[i])
            except:
                print('skipping ' + src[i])


if __name__ == "__main__":
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        print(f"-Start importing subject {s}")
        # old_id = s.replace('sub-','s',1)
        dir1 = os.path.join(orig_dir, str(s))
        dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
        import_func_resting(dir1, dir2, str(s))
        print(f"-Done subject {s}")

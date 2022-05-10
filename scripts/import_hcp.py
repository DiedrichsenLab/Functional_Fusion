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
# base_dir = 'Y:/data'
base_dir = '/Volumes/diedrichsen_data$/data'
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
   
    src=[]
    dest =[]

    for ss in [1, 2]:
         # add a folder for session
         dest_dir = dest_dir + f"ses-{ss:02}"
         # Make the destination directory
         Path(dest_dir).mkdir(parents=True, exist_ok=True)

         # move data into the corresponding session folder
         src.append(f'/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii')
         dest.append(f'/sub-{participant_id}_ses-{ss:02}_task-rest_space-fsLR32k_run-01_dtseries.nii')

         src.append(f'/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii')
         dest.append(f'/sub-{participant_id}_ses-{ss:02}_task-rest_space-fsLR32k_run-02_dtseries.nii')
         for i in range(len(src)):
             try:
                 shutil.copyfile(source_dir+'/MNINonLinear/Results'+src[i], dest_dir+dest[i])
             except:
                 print('skipping ' + src[i])


if __name__ == "__main__":
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        print(f"-Doing subject {s}")
        # old_id = s.replace('sub-','s',1)
        dir1 = os.path.join(orig_dir, str(s))
        dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
        import_func_resting(dir1, dir2, str(s))

if __name__ == "__main__":
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        print(f"-Doing subject {s}")
        # old_id = s.replace('sub-','s',1)
        dir1 = os.path.join(orig_dir, str(s))
        dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
        import_func_resting(dir1, dir2, str(s))
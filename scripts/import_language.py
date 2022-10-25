#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for importing the Language data set to general format.
Created on 4/25/2022 at 12:18 PM
Author: dzhi
"""
import pandas as pd
import shutil
from pathlib import Path
import os
import sys

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Language'
if sys.platform == "win32":  # for windows user
    base_dir = 'Y:/data'

orig_dir = os.path.join(base_dir, 'raw')
target_dir = os.path.join(base_dir, 'derivatives')

def import_func_task(source_dir, dest_dir, participant_id):
    """Imports the Language preprocessed task state files
       into a BIDS/derivative structure
    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        participant_id (str): ID of participant
    """
    task_name = ['langloc_S-N', 'MDloc_H-E',
                'ProdE1_NProd', 'ProdE1_SComp', 'ProdE1_SProd', 'ProdE1_VisEvSem', 'ProdE1_WComp', 'ProdE1_WProd',
                'ProdE2_SProd', 'ProdE2_WProd',
                'ProdE3_NProd', 'ProdE3_SComp', 'ProdE3_SProd', 'ProdE3_VisEvSem', 'ProdE3_WComp', 'ProdE3_WProd']
    sessions = ['ses-01', 'ses-01',
                'ses-02', 'ses-02', 'ses-02', 'ses-02', 'ses-02', 'ses-02',
                'ses-03', 'ses-03',
                'ses-04', 'ses-04', 'ses-04', 'ses-04', 'ses-04', 'ses-04']
    
    for task,task_n in enumerate(task_name):
        dest_dir = os.path.join(dest_dir, sessions[task])
        

        if 'loc' in task_n:
            new_task_n = task_n.split('_')[0]
        elif 'ProdE' in task_n:
            new_task_n = task_n.split('_')[1]

        task_id = sessions[task]

        # move data into the corresponding session folder
        src = (f'/{participant_id}_{task_n}_t.nii')
        dest = (
            f'/{participant_id}_run-99_con-{task_id:02}_space-MNI.nii')

        # Make the destination directory
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        try:
            shutil.copyfile(source_dir+src, 
                    dest_dir+dest)
        except:
            print('skipping ' + src)


if __name__ == "__main__":
    T = pd.read_csv(base_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        print(f"-Start importing subject {s}")
        # old_id = s.replace('sub-','s',1)
        source_dir = os.path.join(orig_dir, '%s/' % str(s))
        dest_dir = os.path.join(target_dir, '%s/estimates/' % str(s))
        import_func_task(source_dir, dest_dir, str(s))
        print(f"-Done subject {s}")

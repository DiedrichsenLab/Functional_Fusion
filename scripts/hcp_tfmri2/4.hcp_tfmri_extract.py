"""
Script for importing the Language dataset to general format.

Created Dec 2023
Author: Bassel Arafat
"""

import pandas as pd
from pathlib import Path
# import mat73
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetHcpTask


base_dir = 'Y:/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'


data_dir = base_dir + '/HCPur100'
atlas_dir = base_dir + '/Atlases'


types = ['CondHalf']
atlases  = ['fs32k']
session_list = ['ses-task2']


dataset = DataSetHcpTask(data_dir)
for ses in session_list:
    print(f'extracting session {ses}')
    participants_tsv = pd.read_csv(f'{data_dir}/participants.tsv',sep = '\t')
    subj_list = participants_tsv['participant_id'].tolist()[:14]

    # remove 'sub-105115'
    # subj_list.remove('sub-105115')
    subj_list = ['sub-101309','sub-103111','sub-103414','sub-103818','sub-105014','sub-110411','sub-117122','sub-118730','sub-123117','sub-124422','sub-127630','sub-128632','sub-130013']

    for type in types:
        print(f'extracting type {type}')
        for atlas in atlases:
            print(f'extracting atlas: {atlas}')
            dataset.extract_all(ses_id = ses,type = type, atlas = atlas, smooth=None, subj = subj_list)
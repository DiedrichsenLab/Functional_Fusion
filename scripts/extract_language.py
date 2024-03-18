"""
Script for importing the Language dataset to general format.

Created Dec 2023
Author: Bassel Arafat
"""

import pandas as pd
from pathlib import Path
# import mat73
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetLanguage


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'


data_dir = base_dir + '/Language'
atlas_dir = base_dir + '/Atlases'


types = ['CondRun']
atlases  = ['MNISymC2']
session_list = ['ses-02']


LL_dataset = DataSetLanguage(data_dir)
for ses in session_list:

    print(f'extracting session {ses}')
    participants_tsv = pd.read_csv(f'{data_dir}/participants.tsv',sep = '\t')
    if ses  == 'ses-01':
        filtered_participants = participants_tsv[participants_tsv['ses'].str.contains('loc')]
        subj_list = filtered_participants['participant_id'].tolist()

    elif ses == 'ses-02':
        filtered_participants = participants_tsv[participants_tsv['ses'].str.contains('sen')]
        subj_list = filtered_participants['participant_id'].tolist()
    else:
        raise Exception('wrong session values')


    for type in types:
        print(f'extracting type {type}')
        for atlas in atlases:
            print(f'extracting atlas: {atlas}')
            LL_dataset.extract_all(ses_id = ses,type = type, atlas = atlas, subj= subj_list)
 

"""
Script for importing the Language dataset to general format.

Created Dec 2023
Author: Bassel Arafat
"""

import pandas as pd
from pathlib import Path
# import mat73
from Functional_Fusion.dataset import DataSetLanguage


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion_new'


data_dir = base_dir + '/Language'
atlas_dir = base_dir + '/Atlases'


types = ['CondRun','CondHalf','CondAll']
atlases  = ['fs32k','MNISymC3']
session_list = ['ses-localizer','ses-localizerfm']

subj = 'all'

LL_dataset = DataSetLanguage(data_dir)
for ses in session_list:
    print(f'extracting session {ses}')
    for type in types:
        print(f'extracting type {type}')
        for atlas in atlases:
            print(f'extracting atlas: {atlas}')
            LL_dataset.extract_all(ses_id = ses,type = type, atlas = atlas, smooth=None, subj=subj)
 

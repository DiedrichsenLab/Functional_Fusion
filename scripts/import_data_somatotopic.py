"""
Script for importing the Pontine dataset to general format.

Created Sep 2022
Author: caro nettekoven
"""

import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
from import_data import *
from Functional_Fusion.dataset import DataSetSomatotopic

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'




if __name__ == '__main__':
    
    src_base_dir = base_dir + '/Cerebellum//Somatotopic'
    dest_base_dir = base_dir + '/FunctionalFusion/Somatotopic'

    dataset = DataSetSomatotopic(dest_base_dir)
    T = dataset.get_participants()


    for p, participant_id in enumerate(T.participant_id):
        print(dataset.get_data_fnames(participant_id))
        # '01', '03', '04', '07', '95', '96', '97',

        # --- Importing Estimates ---
        #source_dir = '{}/GLM_firstlevel_2/S{}/'.format(src_base_dir, #participant_id)
        # dest_dir = '{}/derivatives/sub-{}/estimates/ses-01'.format(dest_base_dir, participant_id)
        #subj_id = 'sub-{}'.format(participant_id)
        #ses_id = 'ses-01'
        #import_spm_glm(source_dir, dest_dir, subj_id, ses_id)



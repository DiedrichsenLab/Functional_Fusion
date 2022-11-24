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

if __name__ == '__main__':
    src_base_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Somatotopic/'
    dest_base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Somatotopic/'
    for participant_id in [ '02']:
        # '01', '03', '04', '07', '95', '96', '97',

        # # # --- Importing SUIT ---
        # source_dir = '{}/suit/anatomicals/S{}'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/suit'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_suit(source_dir,dest_dir,anat_name,'sub-' + participant_id)

        # # --- Importing ANAT ---
        # source_dir = '{}/anatomicals/S{}'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_anat(source_dir,dest_dir,anat_name,participant_id)

        # --- Importing Freesurfer ---
        source_dir = '{}/surfaceWB/data/S{}/'.format(src_base_dir, participant_id)
        dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        old_id = 'S{}'.format(participant_id)
        new_id = 'sub-{}'.format(participant_id)
        import_freesurfer(source_dir,dest_dir,old_id,new_id)

        # --- Importing Estimates ---
        #source_dir = '{}/GLM_firstlevel_2/S{}/'.format(src_base_dir, #participant_id)
        # dest_dir = '{}/derivatives/sub-{}/estimates/ses-01'.format(dest_base_dir, participant_id)
        #subj_id = 'sub-{}'.format(participant_id)
        #ses_id = 'ses-01'
        #import_spm_glm(source_dir, dest_dir, subj_id, ses_id)



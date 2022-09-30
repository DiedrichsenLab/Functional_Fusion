"""
Script for importing IBC dataset to general format.

Created Sep 2022
"""

import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
from import_data import *

if __name__ == '__main__':
    src_base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/ibc/raw'
    dest_base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/ibc/'
    T = pd.read_csv(dest_base_dir + 'participants.tsv',delimiter='\t')
    
    for pt in T.participant_id:

        # # --- Importing SUIT ---
        # source_dir = '{}/suit/anatomicals/S{}'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/suit'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_suit(source_dir,dest_dir,anat_name,'sub-' + participant_id)
        # pass

        # # --- Importing ANAT ---
        source_dir = '{}/{}/anat/'.format(src_base_dir, pt)
        dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, pt)
        anat_name = f'{pt}_T1w'
        import_anat(source_dir,dest_dir,anat_name,pt)
        pass


        # # --- Importing ANAT ---
        # source_dir = '{}/anatomicals/S{}'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_anat(source_dir,dest_dir,anat_name,participant_id)

        # # --- Importing Freesurfer ---
        # source_dir = '{}/surfaceFreesurfer/S{}/surf'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        # new_id = 'sub-{}'.format(participant_id)
        # import_freesurfer(source_dir,dest_dir,new_id)

        # --- Importing Estimates ---
        #source_dir = '{}/GLM_firstlevel_2/S{}/'.format(src_base_dir, #participant_id)
        # dest_dir = '{}/derivatives/sub-{}/estimates/ses-01'.format(dest_base_dir, participant_id)
        #subj_id = 'sub-{}'.format(participant_id)
        #ses_id = 'ses-01'
        #import_spm_glm(source_dir, dest_dir, subj_id, ses_id)


# /Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/GLM_firstlevel_2/S01/beta_0048.nii
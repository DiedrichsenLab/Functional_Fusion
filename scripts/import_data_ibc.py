"""
Script for importing IBC dataset to general format.

Authors: Joern Diedrichsen, Ana Luisa Pinho

Created: September 2022
Last update: October 2022

"""

import os
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
from import_data import *


if __name__ == '__main__':
    base_dir = os.path.join(os.path.expanduser('~'), 'diedrichsen_data/data')
    src_base_dir = os.path.join(base_dir, 'ibc')
    dest_base_dir = os.path.join(base_dir, 'FunctionalFusion/ibc')
    T = pd.read_csv(os.path.join(dest_base_dir, 'participants.tsv'),
                    delimiter='\t')
    for pt in T.participant_id:

        # # --- Importing SUIT ---
        source_dir = '{}/raw/{}/suit'.format(src_base_dir, pt)
        dest_dir = '{}/derivatives/{}/suit'.format(dest_base_dir, pt)
        anat_name = f'{pt}_T1w'
        import_suit(source_dir, dest_dir, anat_name, pt)
        pass

        # # --- Importing ANAT ---
        source_dir = '{}/raw/{}/anat'.format(src_base_dir, pt)
        dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, pt)
        anat_name = f'{pt}_T1w'
        import_anat(source_dir, dest_dir, anat_name, pt)
        pass

        # # --- Importing Freesurfer ---
        source_dir = '{}/surfaceWB/data/{}'.format(src_base_dir, pt)
        dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, pt)
        sub_id = '{}'.format(pt)
        import_freesurfer(source_dir, dest_dir, sub_id, sub_id)

        # --- Importing Estimates ---
        #source_dir = '{}/GLM_firstlevel_2/S{}/'.format(src_base_dir, #participant_id)
        # dest_dir = '{}/derivatives/sub-{}/estimates/ses-01'.format(dest_base_dir, participant_id)
        #subj_id = 'sub-{}'.format(participant_id)
        #ses_id = 'ses-01'
        #import_spm_glm(source_dir, dest_dir, subj_id, ses_id)

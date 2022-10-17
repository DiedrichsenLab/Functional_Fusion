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
# import mat73
import numpy as np
# import scipy.io as sio
from import_data import import_suit, import_anat, import_freesurfer


# ######################### FUNCTIONS ###################################


def import_ibc_anatderivatives(anat_type, source_basedir, destination_basedir,
                               participant):
    anat_name = f'{pt}_T1w'
    sub_id = '{}'.format(pt)
    if anat_type == 'suit':
        source_dir = '{}/raw/{}/suit'.format(source_basedir, pt)
        dest_dir = '{}/derivatives/{}/suit'.format(destination_basedir,
                                                   participant)
        import_suit(source_dir, dest_dir, anat_name, participant)
    elif anat_type == 'anatomical':
        source_dir = '{}/raw/{}/anat'.format(source_basedir, pt)
        dest_dir = '{}/derivatives/{}/anat'.format(destination_basedir,
                                                   participant)
        import_anat(source_dir, dest_dir, anat_name, participant)
    else:
        assert anat_type == 'freesurfer'
        source_dir = '{}/surfaceWB/data/{}'.format(source_basedir, pt)
        dest_dir = '{}/derivatives/{}/anat'.format(destination_basedir,
                                                   participant)
        import_freesurfer(source_dir, dest_dir, sub_id, sub_id)


def import_ibc_glm(participant_id,sess_id):
    # --- Importing Estimates ---
    source_dir = '{}//S{}/'.format(src_base_dir, participant_id)
    dest_dir = '{}/derivatives/sub-{}/estimates/ses-01'.format(dest_base_dir, participant_id)

    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    # Load reginfo file from the source dir: 
    info_name = (source_dir + f'/{participant_id}_{ses_id}_reginfo.tsv')
    info = pd.read_csv(info_name,delimiter='\t')
    N = info.shape[0]

    n_runs = np.max(info.run) 
    # Ensure that run number is an integer value
    info.to_csv(dest_dir + f'/{sub_id}_{sess_id}_reginfo.tsv', sep='\t')

    # Prepare beta files for transfer
    src = []
    dest = []
    for i,r in info.iterrows():
        src.append(f'run-{r.run:02d}/beta_{r.reg_num:04d}.nii')
        dest.append(f'/{sub_id}_{sess_id}_run-{r.run:02d}_' +
                    'reg-{r.reg_id:02d}_beta.nii')
    
    for r in np.unique(info.run):
        # Mask
        run_id = f'run-{r:02d}'
        src.append('/{run_id}/mask.nii')
        dest.append(f'/{sub_id}_{sess_id}_{run_id}_mask.nii')
        src.append('/{run_id}/resms.nii')
        dest.append(f'/{sub_id}_{sess_id}_{run_id}_resms.nii')

    # Average Mask and resms.nii across runs and write them out as 
    # f'{dest_dir}/{sub_id}_{sess_id}_resms.nii')
    # f'{dest_dir}/{sub_id}_{sess_id}_mask.nii')


    # Copy those files over
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i], dest_dir+dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])

        #import_spm_glm(source_dir, dest_dir, subj_id, ses_id)


# ######################### INPUTS ######################################

base_dir = os.path.join(os.path.expanduser('~'), 'diedrichsen_data/data')
src_base_dir = os.path.join(base_dir, 'ibc')
dest_base_dir = os.path.join(base_dir, 'FunctionalFusion/IBC')
session_names = ['archi']

# ########################## RUN ########################################

if __name__ == '__main__':
    T = pd.read_csv(os.path.join(dest_base_dir, 'participants.tsv'),
                    delimiter='\t')
    for pt in T.participant_id:

        # # --- Importing SUIT ---
        import_ibc_anatderivatives('suit', src_base_dir, dest_base_dir, pt)

        # # --- Importing ANAT ---
        import_ibc_anatderivatives('anatomical', src_base_dir, dest_base_dir,
                                   pt)

        # # --- Importing Freesurfer ---
        import_ibc_anatderivatives('freesurfer', src_base_dir, dest_base_dir,
                                   pt)

        # # --- Importing Estimates ---


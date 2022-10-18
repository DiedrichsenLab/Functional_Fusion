"""
Script for importing IBC dataset to general format.

Authors: Joern Diedrichsen, Ana Luisa Pinho

Created: September 2022
Last update: October 2022

"""

import os
import glob
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


def import_ibc_glm(source_basedir, destination_basedir, participant_id,
                   sess_id):
    # --- Importing Estimates ---
    source_dir = '{}/derivatives/{}/estimates/ses-{}'.format(
        source_basedir, participant_id, sess_id)
    dest_dir = '{}/derivatives/{}/estimates/ses-{}'.format(
        destination_basedir, participant_id, sess_id)

    # Create destination file if it does not exist
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    info_name = pt + '_ses-' + sess_id + '_reginfo.tsv'
    for info_path in glob.glob(os.path.join(source_dir, info_name)):
        # Copy reginfo files
        shutil.copyfile(info_path, os.path.join(dest_dir, info_name))
        # Load reginfo file from the source dir: 
        info = pd.read_csv(info_path, sep='\t')
        n_runs = np.max(info.run) 
        # Prepare beta files for transfer
        src = []
        dest = []
        for i, r in info.iterrows():
            src.append(f'run-{r.run:02d}/beta_{r.reg_id:04d}.nii')
            dest_fname = participant_id + '_ses-' + sess_id + \
                '_run-%02d' % r.run + '_reg-%02d' % r.reg_num + '_beta.nii'
            dest.append(dest_fname)
        # Delete any pre-existing .nii file from destination folder
        if glob.glob(dest_dir + '/*.nii'):
            for ni in glob.glob(dest_dir + '/*.nii'):
                os.remove(ni)
        # Copy those files over
        for i in range(len(src)):
            try:
                shutil.copyfile(os.path.join(source_dir, src[i]),
                                os.path.join(dest_dir, dest[i]))
            except FileNotFoundError:
                print('skipping ' + os.path.join(source_dir, src[i]))

        # Saves SPM.mat file as a npy file
        # import_spm_glm(source_dir, dest_dir, subj_id, ses_id)

    # # for r in np.unique(info.run):
    # #     # Mask
    # #     run_id = f'run-{r:02d}'
    # #     src.append('/{run_id}/mask.nii')
    # #     dest.append(f'/{sub_id}_{sess_id}_{run_id}_mask.nii')
    # #     src.append('/{run_id}/resms.nii')
    # #     dest.append(f'/{sub_id}_{sess_id}_{run_id}_resms.nii')

    # # Average Mask and resms.nii across runs and write them out as 
    # # f'{dest_dir}/{sub_id}_{sess_id}_resms.nii')
    # # f'{dest_dir}/{sub_id}_{sess_id}_mask.nii')


# ######################### INPUTS ######################################

base_dir = os.path.join(os.path.expanduser('~'), 'diedrichsen_data/data')
src_base_dir = os.path.join(base_dir, 'ibc')
dest_base_dir = os.path.join(base_dir, 'FunctionalFusion/IBC')
sessions = ['archi']

# ########################## RUN ########################################

if __name__ == '__main__':
    T = pd.read_csv(os.path.join(dest_base_dir, 'participants.tsv'),
                    delimiter='\t')
    # for pt in T.participant_id:
    for pt in ['sub-01']:

        # # # --- Importing SUIT ---
        # import_ibc_anatderivatives('suit', src_base_dir, dest_base_dir, pt)

        # # # --- Importing ANAT ---
        # import_ibc_anatderivatives('anatomical', src_base_dir, dest_base_dir,
        #        ^                    pt)

        # # # --- Importing Freesurfer ---
        # import_ibc_anatderivatives('freesurfer', src_base_dir, dest_base_dir,
        #                            pt)

        # # --- Importing Estimates ---
        for session in sessions:
            import_ibc_glm(src_base_dir, dest_base_dir, pt, session)

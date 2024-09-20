"""
Script for importing the Language localizer dataset to general format.

Created Dec 2023
Author: Bassel Arafat
"""

import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
import Functional_Fusion.import_data as id

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'


def import_spm_glm(source_dir, dest_dir, sub_id, sess_id):
    """
    Imports the output of the SPM GLM with an SPM_info.mat
    structure into BIDS deriviatie (Functional Fusion) framework.
    It assumes that a single GLM corresponds to single session.

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that
                           subject / session
        new_id (_type_): New name for the subject
        info_dict (_type_): Dictionary with the old field names and the
			                new field names for the information
    """

    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    # Generate new dictionary from SPM info
    D = pd.read_csv(source_dir + f'/{sub_id}_reginfo.tsv',sep='\t')

    N = len(D)

    if 'reg_id' not in D.columns:
        n = sum(D['run'] == 1)
        D['reg_num'] = np.arange(N)
        D['reg_id'] = D['reg_num'] % n

    D.to_csv(dest_dir + f'/{sub_id}_{sess_id}_reginfo.tsv', sep='\t')

    # Prepare beta files for transfer
    src = []
    dest = []
    for i in range(N):
        src.append(f'/beta_{i+1:04d}.nii')
        dest.append(f'/{sub_id}_{sess_id}_run-{D.run[i]:02}_' +
                    f'reg-{D.reg_id[i]:02d}_beta.nii')
    # Mask
    src.append('/mask.nii')
    dest.append(f'/{sub_id}_{sess_id}_mask.nii')

    # ResMS
    src.append('/resms.nii')
    dest.append(f'/{sub_id}_{sess_id}_resms.nii')

    # Copy those files over
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i], dest_dir+dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])

if __name__ == '__main__':
    # src_base_dir = base_dir + '/Cerebellum/Language/Language_7T'
    # dest_base_dir = base_dir + '/FunctionalFusion/Language'
    # ses = 'ses-02'
    # for participant_id in ['sub-06']:
    #     # '01', '03', '04', '07', '95', '96', '97',
    #     # # --- Importing SUIT ---
    #     source_dir = '{}/suit/anatomicals/{}'.format(src_base_dir, participant_id)
    #     dest_dir = '{}/derivatives/{}/suit'.format(dest_base_dir, participant_id)
    #     anat_name = 'anatomical'
    #     id.import_suit(source_dir,dest_dir,anat_name, participant_id)

    #     # # # # --- Importing ANAT ---
    #     source_dir = '{}/anatomicals/{}'.format(src_base_dir, participant_id)
    #     dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, participant_id)
    #     anat_name = 'anatomical'
    #     id.import_anat(source_dir,dest_dir,anat_name,participant_id)

    #     # # --- Importing Freesurfer ---
    #     source_dir = '{}/surfaceWB/data/{}/'.format(src_base_dir, participant_id)
    #     dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, participant_id)
    #     old_id = '{}'.format(participant_id)
    #     new_id = '{}'.format(participant_id)
    #     id.import_freesurfer(source_dir,dest_dir,old_id,new_id)

    #     # --- Importing Estimates ---
        
    #     subj_id = '{}'.format(participant_id)
    #     ses_id = 'ses-02'
    #     id.import_spm_glm(source_dir, dest_dir, subj_id, ses_id)

    #     # # Importing design matrix
    #     source_dir = f'{src_base_dir}/GLM_firstlevel/glm_21/{participant_id}'
    #     dest_dir = f'{dest_base_dir}/derivatives/{participant_id}/estimates/{ses}'
    #     id.import_spm_designmatrix(source_dir, dest_dir, participant_id, ses)


    # session='rest'
    session='localizer_cond'
    session_orig = '04' if session == 'rest' else '01'

    dest_dir = base_dir + '/FunctionalFusion/Language/derivatives/{sub}/estimates/' + f'ses-{session}' '/{sub}_' + f'ses-{session}'
    T = pd.read_csv(base_dir + '/FunctionalFusion/Language/participants.tsv', delimiter='\t')
    participants = T[T[f'ses-rest'] == 1].participant_id
    runs =[f'{run:02d}' for run in np.arange(1, 8)    ]

    fix=False
    if fix:
        src_stem = base_dir + '/Cerebellum/Language/Language_7T/imaging_data_fix/{sub}/' + f'ses-{session_orig}' + '/{sub}_' + f'ses-{session_orig}'
        file_ending = '_run-{run}_fix.nii'
    else:
        src_stem = base_dir + '/Cerebellum/Language/Language_7T/imaging_data/{sub}/' + f'ses-{session_orig}' + '/r{sub}_' + f'ses-{session_orig}'
        file_ending = '_run-{run}.nii'
        
    for s in participants:
        src = src_stem.format(sub=s) + file_ending
        dest = dest_dir.format(sub=s) + file_ending
        mask_file = base_dir + '/Cerebellum/Language/Language_7T/imaging_data/{sub}/ses-01/rmask_noskull.nii'.format(sub=s)
        id.import_tseries(src, dest, s, f'ses-{session_orig}', runs, mask_file=mask_file)

        pass



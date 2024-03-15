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
from Functional_Fusion.import_data import *

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data'


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
    src_base_dir = base_dir + '/Cerebellum/Language/Language_7T'
    dest_base_dir = base_dir + '/FunctionalFusion/Language'
    participant_list = ['sub-01','sub-02','sub-03','sub-04','sub-06','sub-07','sub-08']


    for participant_id in participant_list:
        print(f'importing data for {participant_id}')

        print('getting suit stuff')
        # # --- Importing SUIT ---
        source_dir = '{}/suit/anatomicals/{}'.format(src_base_dir, participant_id)
        dest_dir = '{}/derivatives/{}/suit'.format(dest_base_dir, participant_id)
        anat_name = 'anatomical'
        import_suit(source_dir,dest_dir,anat_name, participant_id)

        print('getting anat stuff')
        # # # # # --- Importing ANAT ---
        source_dir = '{}/anatomicals/{}'.format(src_base_dir, participant_id)
        dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, participant_id)
        anat_name = 'anatomical'
        import_anat(source_dir,dest_dir,anat_name,participant_id)

        print('getting surface stuff')
        # # # --- Importing Freesurfer ---
        source_dir = '{}/surfaceWB/data/{}/'.format(src_base_dir, participant_id)
        dest_dir = '{}/derivatives/{}/anat'.format(dest_base_dir, participant_id)
        old_id = '{}'.format(participant_id)
        new_id = '{}'.format(participant_id)
        import_freesurfer(source_dir,dest_dir,old_id,new_id)

        
        tsv = pd.read_csv(f'{dest_base_dir}/participants.tsv', sep='\t')
        participant_ses= tsv[tsv.participant_id == participant_id].ses.iloc[0]
        if participant_ses == 'loc':
            ses_list = ['ses-01']
            glm_list = ['glm_11']
        elif participant_ses == 'sen':
            ses_list = ['ses-02']
            glm_list = ['glm_22']
        elif participant_ses == 'loc_sen':
            ses_list = ['ses-01','ses-02']
            glm_list = ['glm_11','glm_22']


        for i,ses in enumerate(ses_list):
            glm = glm_list[i]
            print(f'getting glm for {ses}')
            # # --- Importing Estimates ---
            source_dir = '{}/GLM_firstlevel/{}/{}/'.format(src_base_dir, glm, participant_id)
            dest_dir = '{}/derivatives/{}/estimates/{}'.format(dest_base_dir, participant_id,ses)
            subj_id = '{}'.format(participant_id)
            import_spm_glm(source_dir, dest_dir, subj_id, ses)

            print(f'getting design matrix for {ses}')
            # # # Importing design matrix
            source_dir = f'{src_base_dir}/GLM_firstlevel/{glm}/{participant_id}'
            dest_dir = f'{dest_base_dir}/derivatives/{participant_id}/estimates/{ses}'
            import_spm_designmatrix(source_dir, dest_dir, participant_id, ses)


        pass



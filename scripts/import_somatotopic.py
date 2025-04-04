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
from Functional_Fusion.dataset import DataSetSomatotopic
import Functional_Fusion.atlas_map as am
import shutil
import nibabel as nb
from copy import deepcopy
import subprocess

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'


src_base_dir = Path(base_dir + '/Cerebellum//Somatotopic')
dest_base_dir = Path(base_dir + '/FunctionalFusion/Somatotopic')


def create_reginfo(log_message=False):
    dataset = DataSetSomatotopic(str(dest_base_dir))

    # Import general info
    info = pd.read_csv(dest_base_dir / 'ses-motor_reginfo.tsv', sep='\t')

    # Import scan-specific info (missing runs)
    missing_scans = pd.read_csv(src_base_dir / 'missing_scans.txt', names=['scan'], header=None)
    missing_scans[['orig_id', 'session', 'run']] = missing_scans.scan.str.split(
        "_", expand=True)
    
    T = dataset.get_participants()
    for _, id in T.iterrows():
        print(f'Creating reginfo for {id.participant_id}')
        
        # Ammend the reginfo.tsv file from the general file
        reginfo = deepcopy(info)
        reginfo.insert(loc=0, column='sn', value=[
            id.participant_id] * info.shape[0])
        
        for session in info.orig_ses.unique():
            missing = (id.orig_id == missing_scans.orig_id) & (
                f'sess{session:02d}' == missing_scans.session)
            if np.any(missing):
                for _, miss in missing_scans[missing].iterrows():
                    run = int(miss.run.split('MOTOR')[1])
                    if log_message:
                        print(f'Missing scan {reginfo[reginfo.orig_run == run]} removed from reginfo')
                    reginfo = reginfo.drop(reginfo[reginfo.orig_run == run].index)
                
        # Make folder
        dest = dest_base_dir / f'derivatives/{id.participant_id}/estimates/ses-motor/{id.participant_id}_ses-motor_reginfo.tsv'
        dest.parent.mkdir(parents=True, exist_ok=True)

        # Save reginfo.tsv file
        reginfo.to_csv(dest, sep='\t', index=False)

def import_estimates(log_message=False):
    dataset = DataSetSomatotopic(str(dest_base_dir))
    T = dataset.get_participants()

    for _, id in T.iterrows():
        print(f'Importing {id.participant_id}')

        # Load reginfo.tsv file
        reginfo = pd.read_csv(dest_base_dir /
                              f'derivatives/{id.participant_id}/estimates/ses-motor/{id.participant_id}_ses-motor_reginfo.tsv', sep='\t')

        resms = []
        for _, info in reginfo.iterrows():
            if log_message:
                print(
                    f'Run {info.run:02d} Beta {info.reg_num:02d}')
            src = src_base_dir / \
                f'raw/Functional/{id.orig_id}/{id.orig_id}_sess{info.orig_ses:02d}_MOTOR{info.orig_run}/pe{info.orig_reg}.nii.gz'
            dest = dest_base_dir / \
                f'derivatives/{id.participant_id}/estimates/ses-motor/{id.participant_id}_ses-motor_run-{info.run:02d}_reg-{info.reg_id:02d}_beta.nii.gz'
            

            # Copy func file to destination folder and rename
            if  ~dest.exists():
                try:
                    shutil.copyfile(src,
                            dest)
                except:
                    print('skipping ' + str(src))
            
            # Unzip because file name ends in .nii.gz
            subprocess.call(
                ['gunzip', '-f', dest])

            
            src = src_base_dir / \
                f'raw/Functional/{id.orig_id}/{id.orig_id}_sess{info.orig_ses:02d}_MOTOR{info.orig_run}/sigmasquareds.nii.gz'
            resms_img = nb.load(src)
            resms.append(resms_img.get_fdata())

        resms = np.mean(resms,axis=0)
        nifti_img = nb.Nifti1Image(dataobj=resms, affine=resms_img.affine)
        outname = dest_base_dir / \
            f'derivatives/{id.participant_id}/estimates/ses-motor/{id.participant_id}_ses-motor_resms.nii'
        nb.save(nifti_img, outname)
            

def import_anat(source_dir, dest_dir, anat_name, participant_id):
    """
    Imports a anatomy folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory (anatomical)
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file
        participant_id (str): ID of participant
    """

    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    src = f'/{anat_name}.nii.gz'
    dest = f'/{participant_id}_T1w.nii'

    try:
        shutil.copyfile(source_dir + src, dest_dir + dest)
        # Unzip because file name ends in .nii.gz
        subprocess.call(
            ['gunzip', '-f', source_dir + src])
    except FileNotFoundError:
        print('skipping ' + src)

    

if __name__ == '__main__':
    # --- Create reginfo ---
    # create_reginfo()
    
    # --- Importing Estimates ---
    # import_estimates()

    # --- Importing Freesurfer ---
    # for participant_id in ['01', '02', '03', '04', '05', '06', '07', '08']:

    #     source_dir = '{}/raw/surfaceWB/data/sub-{}/'.format(
    #         src_base_dir, participant_id)
    #     dest_dir = '{}/derivatives/sub-{}/anat'.format(
    #         dest_base_dir, participant_id)
    #     old_id = 'sub-{}'.format(participant_id)
    #     new_id = 'sub-{}'.format(participant_id)
    #     import_freesurfer(source_dir, dest_dir, old_id, new_id)

    # --- Importing Freesurfer & Anat ---
    
    for participant_id in ['01', '02', '03', '04', '05', '06', '07', '08']:

        source_dir = '{}/raw/sub-{}/anat/'.format(
            src_base_dir, participant_id)
        dest_dir = '{}/derivatives/sub-{}/anat'.format(
            dest_base_dir, participant_id)
        id = 'sub-{}'.format(participant_id)
        anat_name = f'{id}_T1w'
        import_anat(source_dir, dest_dir, anat_name, id)




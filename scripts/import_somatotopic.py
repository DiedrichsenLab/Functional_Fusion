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
    info = pd.read_csv(dest_base_dir / 'ses-all_reginfo.tsv', sep='\t')

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
        
        for session in np.arange(1,5):
            missing = (id.orig_id == missing_scans.orig_id) & (
                f'sess{session:02d}' == missing_scans.session)
            if np.any(missing):
                for _, miss in missing_scans[missing].iterrows():
                    run = int(miss.run.split('MOTOR')[1])
                    if log_message:
                        print(f'Missing scan {reginfo[reginfo.run == run]} removed from reginfo')
                    reginfo = reginfo.drop(reginfo[reginfo.run == run].index)
                

            # Save reginfo.tsv file
            reginfo.to_csv(dest_base_dir /
                           f'derivatives/{id.participant_id}/estimates/ses-{session:02d}/{id.participant_id}_ses-{session:02d}_reginfo.tsv', sep='\t', index=False)

def import_data(log_message=False):
    dataset = DataSetSomatotopic(str(dest_base_dir))
    T = dataset.get_participants()

    for _, id in T.iterrows():
        print(f'Importing {id.participant_id}')

        for session in np.arange(1,5):
            # Load reginfo.tsv file
            reginfo = pd.read_csv(dest_base_dir /
                        f'derivatives/{id.participant_id}/estimates/ses-{session:02d}/{id.participant_id}_ses-{session:02d}_reginfo.tsv', sep='\t')
            
            resms = []
            for _, info in reginfo.iterrows():
                if log_message:
                    print(
                    f'Session {session:02d} Run {info.run:02d} Beta {info.reg_num:02d}')
                src = src_base_dir / \
                    f'raw/Functional/{id.orig_id}/{id.orig_id}_sess{session:02d}_MOTOR{info.run}/pe{info.pe_id}.nii.gz'
                dest = dest_base_dir / \
                    f'derivatives/{id.participant_id}/estimates/ses-{session:02d}/{id.participant_id}_ses-{session:02d}_run-{info.run:02d}_reg-{info.reg_num:02d}_beta.nii.gz'
                
                # Make folder
                dest.parent.mkdir(parents=True, exist_ok=True)

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
                    f'raw/Functional/{id.orig_id}/{id.orig_id}_sess{session:02d}_MOTOR{info.run}/sigmasquareds.nii.gz'
                resms_img = nb.load(src)
                resms.append(resms_img.get_fdata())

            resms = np.mean(resms,axis=0)
            nifti_img = nb.Nifti1Image(dataobj=resms, affine=resms_img.affine)
            outname = dest_base_dir / \
                f'derivatives/{id.participant_id}/estimates/ses-{session:02d}/{id.participant_id}_ses-{session:02d}_resms.nii'
            nb.save(nifti_img, outname)
            



        



if __name__ == '__main__':
    # --- Create reginfo ---
    # create_reginfo()
    
    # --- Importing Estimates ---
    import_data()



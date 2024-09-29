# Script for getting predicted timesries for the MDTB data set from super_cerebellum 
import os
import pandas as pd
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetMDTB
import Functional_Fusion.dataset as ds
import nibabel as nb
from nitools import spm
import scripts.fusion_paths as paths
import Functional_Fusion.connectivity as conn
import matplotlib.pyplot as plt


base_dir = paths.set_base_dir()
atlas_dir = paths.set_atlas_dir(base_dir)
dname = 'MDTB'
data_dir = paths.set_fusion_dir(base_dir)
mdtb_dir = data_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

from pathlib import Path
import nibabel as nb

def save_data(dataset, subject, session, space, atlas, data, data_info, data_type, binfo=None):
    """
    Function to save the data from re-running the GLM as CIFTI files

    Parameters:
    - dataset (object): Dataset object.
    - subject (str): Subject identifier.
    - session (str): Session identifier.
    - space (str): The space identifier (e.g., MNI).
    - atlas (object): Atlas object to convert data to cifti format.
    - data (array): Data to be saved (e.g., betas, residuals, predicted timeseries).
    - data_info (object): Data information object to save as accompanying .tsv file.
    - data_type (str): Type of data being saved (e.g., FixBeta, FixResiduals, FixPredicted, FixAdjusted).
    - binfo (optional object): Optional information object for the beta data.
    """
    print(f'Saving {data_type} for {subject}')
    
    # Convert data to CIFTI format
    C = atlas.data_to_cifti(data, data_info.names)
    
    # Create the destination directory
    dest_dir = dataset.data_dir.format(subject)
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    
    # Save the .nii file
    nii_filename = f'{dest_dir}/{subject}_space-{space}_{session}_{data_type}.dscalar.nii'
    nb.save(C, nii_filename)
    
    # Save the .tsv file
    tsv_filename = f'{dest_dir}/{subject}_{session}_info-{data_type}.tsv'
    data_info.to_csv(tsv_filename, sep='\t', index=False)
    


def rerun_glms(dataset, glm_path, subjects, session, space):
    """Function to rerun GLMs for a given session and space, and save betas, residuals, predicted, and adjusted timeseries.

    Parameters:
    - dataset (object): Dataset object.
    - glm_path (str): Path to the GLM directory.
    - subjects (list): List of subject identifiers.
    - session (int): Session identifier.
    - space (str): The space identifier (e.g., MNI).
    """
    # Get the atlas object
    atlas, _ = am.get_atlas(space)
    
    for subj in range(len(subjects)):
        # Loop over subjects, load SPM.mat file, and get betas, residuals, and adj data
        subject = subjects[subj]
        subject_name_orig = f"s{subject.split('-')[1]}"
        subject_path = glm_path.format(subject_name_orig=subject_name_orig)
        
        # Load SPM.mat file
        subject_spm = spm.SpmGlm(subject_path)
        subject_spm.get_info_from_spm_mat()

        # Load functional data
        data_raw, data_info = dataset.get_data(ses_id=session, type='FixTseries', space=space, subj=[subj])

        # Rerun GLM
        beta, binfo, _, data_hat, data_adj, residuals = subject_spm.rerun_glm(data_raw[0])       

        # Save Betas
        binfo['cond_name'] = binfo['reg_name'].str.split('*').str[0]
        beta = beta[binfo['cond_name'] != 'Instruct', :] # Remove instructions from data
        binfo = binfo[binfo['cond_name'] != 'Instruct'] # Remove instructions from info
        beta_info = dataset.get_info(ses_id=session, type='CondRun', subj=[subject])
        if binfo.cond_name.tolist() == beta_info.cond_name.tolist():
            print('Condition names match')
            binfo = beta_info
        else:
            Warning('Condition names do not match. Saving beta info from SPM.mat file without any additional info')
        
        save_data(dataset, subject, session, space, atlas, beta, binfo, 'FixBeta', binfo)

        # Save Residuals
        save_data(dataset, subject, session, space, atlas, residuals, data_info, 'FixResiduals')

        # Save predicted timeseries
        save_data(dataset, subject, session, space, atlas, data_hat, data_info, 'FixPredicted')

        # Save adjusted timeseries
        save_data(dataset, subject, session, space, atlas, data_adj, data_info, 'FixAdjusted')


if __name__ == "__main__":
    mdtb_dataset = DataSetMDTB(mdtb_dir)

    # -- Get betas and residual timeseries for task session from Fix-cleaned data --
    T = pd.read_csv(
        mdtb_dir + '/participants.tsv', delimiter='\t')
    subject_subset = T.participant_id[T['ses-rest'] == 1].tolist()

    # Setttings
    space='MNISymC3'
    session_idx = 2
    session=f'ses-s{session_idx}'
    glm_path = base_dir + f'Cerebellum/super_cerebellum/sc{session_idx}/' 'GLM_firstlevel_7/{subject_name_orig}/'
    
    rerun_glms(mdtb_dataset, glm_path, subject_subset, session, space)
      
    pass
    
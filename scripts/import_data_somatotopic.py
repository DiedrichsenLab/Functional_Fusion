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
from fsl.wrappers import fnirt
from fsl.wrappers import flirt
# from fsl.wrappers.fnirt import applywarp
from fsl.wrappers.flirt import applyxfm
from fsl.utils.platform import platform
from fsl.transform.flirt import readFlirt


def norm_anat(source_dir, dest_dir, anat_name, participant_id):
    """
    Prepares the anatomicals for reconstruction by resampling T1w to MNI.
    If functional data has been analyzed in MNI space, we need to resamle
    T1w anatomicals to MNI space and then reconstruct cortical surfaces
    on the MNI-resampled T1w images.

    Args:
        source_dir (_type_): Directory of the anatomical
        dest_dir (_type_): Destination directory for the anatomical
        anat_name (_type_): Name of anatomical
        old_id (_type_): Old name of the subject
        new_id (_type_): New name for the subject

    """
    anat = Path(source_dir) / anat_name
    outlin = Path(source_dir) / f'sub-{participant_id}_T1w_mni_flirt'
    outnonlin = Path(source_dir) / f'sub-{participant_id}_T1w_mni_fnirt'
    mat = f'{outlin}.mat'
    mni = Path(platform.fsldir) / 'data' / \
        'standard' / 'MNI152_T1_1mm_brain.nii.gz'
    # get initial linear transform
    flirt(anat, mni, out=outlin, omat=mat)
    # get nonlinear transform with initial linear transform
    fnirt(src=anat, ref=mni, aff=mat, iout=outnonlin, cout=f'{outnonlin}_warp')
    pass

    
    
    
    


    pass

def resample_anat(source_dir, dest_dir, anat_name, participant_id):
    """
    Prepares the anatomicals for reconstruction by resampling T1w to MNI.
    If functional data has been analyzed in MNI space, we need to resamle
    T1w anatomicals to MNI space and then reconstruct cortical surfaces
    on the MNI-resampled T1w images.

    Args:
        source_dir (_type_): Directory of the anatomical
        dest_dir (_type_): Destination directory for the anatomical
        anat_name (_type_): Name of anatomical
        old_id (_type_): Old name of the subject
        new_id (_type_): New name for the subject

    """
    anat = Path(source_dir) / anat_name
    out = Path(source_dir) / f'sub-{participant_id}_T1w_mni'
    # warp = Path(source_dir) / \
    #     '{}to_mni_FNIRT.mat.nii.gz'.format(anat_name.split('reorient')[0])
    mat = Path(source_dir) / \
        '{}brain_to_mni.mat'.format(anat_name.split('reorient')[0])
    mni = Path(platform.fsldir) / 'data' / \
        'standard' / 'MNI152_T1_1mm_brain.nii.gz'
    # applywarp(anat, mni, out, warp=warp)
    applyxfm(anat, mni, mat, out)
    

    pass


if __name__ == '__main__':
    src_base_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Somatotopic'
    dest_base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Somatotopic'
    for participant_id in [ '02']:
        # '01', '03', '04', '07', '95', '96', '97',

        # # # --- Importing SUIT ---
        # source_dir = '{}/suit/anatomicals/S{}'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/suit'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_suit(source_dir,dest_dir,anat_name,'sub-' + participant_id)

        # --- Importing ANAT ---
        source_dir = '{}/raw/sub-{}/anat'.format(src_base_dir, participant_id)
        dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        anat_name = 'S{:01d}_mpr_reorient_defaced.nii.gz'.format(int(participant_id))
        norm_anat(source_dir, dest_dir, anat_name, participant_id)
        resample_anat(source_dir, dest_dir, anat_name, participant_id)
        import_anat(source_dir,dest_dir,anat_name,participant_id)

        # # --- Importing Freesurfer ---
        # source_dir = '{}/surfaceWB/data/S{}/'.format(src_base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        # old_id = 'S{}'.format(participant_id)
        # new_id = 'sub-{}'.format(participant_id)
        # import_freesurfer(source_dir,dest_dir,old_id,new_id)

        # --- Importing Estimates ---
        #source_dir = '{}/GLM_firstlevel_2/S{}/'.format(src_base_dir, #participant_id)
        # dest_dir = '{}/derivatives/sub-{}/estimates/ses-01'.format(dest_base_dir, participant_id)
        #subj_id = 'sub-{}'.format(participant_id)
        #ses_id = 'ses-01'
        #import_spm_glm(source_dir, dest_dir, subj_id, ses_id)



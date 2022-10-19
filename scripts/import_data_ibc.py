"""
Script for importing IBC dataset to general format.

Authors: Joern Diedrichsen, Ana Luisa Pinho

Created: September 2022
Last update: October 2022

"""

import os
import glob
import re
import shutil

from pathlib import Path

import numpy as np
import pandas as pd

from import_data import import_suit, import_anat, import_freesurfer

import mat73
import scipy.io as sio
import nibabel as nb


# ######################### FUNCTIONS ###################################


def delete_old(ddir, ext):
    if glob.glob(ddir + '/*.' + ext):
        for f in glob.glob(ddir + '/*.' + ext):
            os.remove(f)


def copy_betas(info_path, subject, sess, sdir, ddir):
    # Load reginfo file from the source dir
    reginfo = pd.read_csv(info_path, sep='\t')
    # Prepare beta files for transfer
    src = []
    dest = []
    for i, r in reginfo.iterrows():
        src.append(f'run-{r.run:02d}/beta_{r.reg_id:04d}.nii')
        dest_fname = subject + '_ses-' + sess + \
            '_run-%02d' % r.run + '_reg-%02d' % r.reg_num + '_beta.nii'
        dest.append(dest_fname)
    # Copy those files over
    for i in range(len(src)):
        try:
            shutil.copyfile(os.path.join(sdir, src[i]),
                            os.path.join(ddir, dest[i]))
        except FileNotFoundError:
            print('skipping ' + os.path.join(sdir, src[i]))


def compute_mean_nifti(vol_paths):
    vols = [nb.load(vol_path) for vol_path in vol_paths]
    X = [vol.get_fdata() for vol in vols]
    Y = np.mean(X, axis=0)
    mean_vol = nb.Nifti1Image(Y, vols[0].affine)

    return mean_vol


def import_spm_ibc_dmtx(sdir, ddir, sub_id, sess_id):
    """
    Imports the SPM design matrix for optimal contrast recombination
    at a later stage. Because python gives some errors when trying to
    read an SPM.mat structure, this requires the design matrix
    information to be extracted from the SPM.mat before, using the
    following matlab code (for every subject):
 
    load('SPM.mat');
    X = SPM.xX.xKXs.X
    save design_matrix.mat X

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that
                           subject / session
        sub_id (_type_): New name for the subject
        sess_id (_type_): ID of the session to import
    """

    # Create new files
    for dm in glob.glob(sdir + '/run-*/design_matrix_unf.mat'):
        X = sio.loadmat(dm)
        DM = X['X']
        runtag = re.match('.*/(.*)/design_matrix_unf.mat', dm).groups()[0]
        fname = sub_id + '_ses-' + sess_id + '_' + runtag + \
            '_designmatrix_unf.npy'
        fpath = os.path.join(ddir, fname)
        # Save current file in the path
        np.save(fpath, DM)


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


def import_ibc_glm(source_basedir, destination_basedir, participant,
                   session_id):

    # Define source and destination paths
    source_dir = '{}/derivatives/{}/estimates/ses-{}'.format(
        source_basedir, participant, session_id)
    destination_dir = '{}/derivatives/{}/estimates/ses-{}'.format(
        destination_basedir, participant, session_id)

    # Create destination file if it does not exist
    Path(destination_dir).mkdir(parents=True, exist_ok=True)

    # Clean destination directory
    delete_old(destination_dir, 'tsv')
    delete_old(destination_dir, 'nii')
    delete_old(destination_dir, 'npy')

    # Reginfo file of the session
    info_name = pt + '_ses-' + session_id + '_reginfo.tsv'
    info_path = os.path.join(source_dir, info_name)

    # Copy reginfo file
    shutil.copyfile(info_path, os.path.join(destination_dir, info_name))

    # Copy beta files to derivatives folder and rename them
    copy_betas(info_path, participant, session_id, source_dir, destination_dir)

    # Compute mean of masks and save in destination folder
    masks_paths = glob.glob(source_dir + '/run-*/mask.nii')
    mean_mask = compute_mean_nifti(masks_paths)
    mean_mask_path = os.path.join(
        destination_dir, participant + '_ses-' + session_id + '_mask.nii')
    nb.save(mean_mask, mean_mask_path)

    # Compute mean of Residual-Sum-of-Squares and save in destination folder
    resmss_paths = glob.glob(source_dir + '/run-*/ResMS.nii')
    mean_resms = compute_mean_nifti(resmss_paths)
    mean_resms_path = os.path.join(
        destination_dir, participant + '_ses-' + session_id + '_resms.nii')
    nb.save(mean_resms, mean_resms_path)

    # Saves SPM.mat file as a .npy file
    import_spm_ibc_dmtx(source_dir, destination_dir, participant, session_id)


# ######################### INPUTS ######################################

base_dir = os.path.join(os.path.expanduser('~'), 'diedrichsen_data/data')
src_base_dir = os.path.join(base_dir, 'ibc')
dest_base_dir = os.path.join(base_dir, 'FunctionalFusion/IBC')

session_group1 = ['archi', 'hcp1', 'hcp2', 'rsvp-language']
session_group2 = ['mtt1', 'mtt2', 'preference', 'tom', 'enumeration', 'self',
                  'clips4', 'lyon1', 'lyon2', 'mathlang',
                  'spatial-navigation']
sessions = session_group1 + session_group1
# sessions = ['archi']

# ########################## RUN ########################################

if __name__ == '__main__':
    T = pd.read_csv(os.path.join(dest_base_dir, 'participants.tsv'),
                    delimiter='\t')
    for pt in T.participant_id:
    # for pt in ['sub-01']:

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
            if pt == 'sub-02' and session in session_group2:
                continue
            else:
                import_ibc_glm(src_base_dir, dest_base_dir, pt, session)

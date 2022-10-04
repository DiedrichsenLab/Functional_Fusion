import shutil
import pandas as pd
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio


def import_suit(source_dir, dest_dir, anat_name, participant_id):
    """
    Imports a suit folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file (submitted to )
        participant_id (str): ID of participant
    """

    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    src.append(f'/c1{anat_name}.nii')
    dest.append(f'/{participant_id}_label-GMc_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/{participant_id}_label-WMc_probseg.nii')
    src.append('/maskbrainSUITGrey.nii')
    dest.append(f'/{participant_id}_desc-cereb_mask.nii')
    src.append(f'/y_{anat_name}_suitdef.nii')
    dest.append(f'/{participant_id}_space-SUIT_xfm.nii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir + src[i], dest_dir + dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])


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

    src = []
    dest = []
    src.append(f'/{anat_name}.nii')
    dest.append(f'/{participant_id}_T1w.nii')
    src.append(f'/c1{anat_name}.nii')
    dest.append(f'/{participant_id}_label-GM_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/{participant_id}_label-WM_probseg.nii')
    src.append(f'/c3{anat_name}.nii')
    dest.append(f'/{participant_id}_label-CSF_probseg.nii')
    src.append(f'/y_{anat_name}.nii')
    dest.append(f'/{participant_id}_space-MNI_xfm.nii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir + src[i], dest_dir + dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])

def import_freesurfer(source_dir,dest_dir,old_id,new_id):
    """Imports the output of a freesurfer reconstruction (and subsequent workbench import).
    Args:
        source_dir (str): Directory of the SPM GLM
        dest_dir (str): Destination directory for that subject / session
        old_id (str): Old subject name
        new_id (str): New name for the subject
    """
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    src.append(f'/{old_id}.L.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_pial.surf.gii')
    src.append(f'/{old_id}.L.white.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_white.surf.gii')
    src.append(f'/{old_id}.R.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_pial.surf.gii')
    src.append(f'/{old_id}.R.white.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_white.surf.gii')
    src.append(f'/{old_id}.L.sulc.32k.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_sulc.shape.gii')
    src.append(f'/{old_id}.R.sulc.32k.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_sulc.shape.gii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i], dest_dir+dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])


def import_spm_glm(source_dir, dest_dir, sub_id, sess_id, info_dict):
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
        info_dict (_type_): Dictonary with the old field names and the
			              new field names for the information
    """

    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    # Generate new dictionary from SPM info
    D = mat73.loadmat(source_dir + '/SPM_info.mat')
    T = {}
    for i in info_dict.items():
        series=D[i[0]]
        if isinstance(series[0],list):
            series=[series[i][0] for i in range(len(series))]
        T[i[1]]=series
    N = len(T[i[1]])
    if 'reg_id' not in T.keys():
        n = sum(T['run'] == 1)
        T['reg_num'] = np.arange(N)
        T['reg_id'] = T['reg_num'] % n

    # Ensure that run number is an integer value
    T['run'] = [int(j) for j in T['run']]
    D = pd.DataFrame(T)
    D.to_csv(dest_dir + f'/{sub_id}_{sess_id}_reginfo.tsv', sep='\t')

    # Prepare beta files for transfer
    src = []
    dest = []
    for i in range(N):
        src.append(f'/beta_{i+1:04d}.nii')
        dest.append(f'/{sub_id}_{sess_id}_run-{D.run[i]:02}_' +
                    'reg-{D.reg_id[i]:02d}_beta.nii')
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


def import_spm_designmatrix(source_dir, dest_dir, sub_id, sess_id):
    """
    Imports the SPM design matrix for optimal contrast recombination
    at a later stage. Because python gives some errors when trying to
    read an SPM.mat structure, this requires the design matrix
    information to be extracted from the SPM.mat before, using the
    following matlab code (for every subject):
	"""
    load('SPM.mat');
    X = SPM.xX.xKXs.X
    save design_matrix.mat X

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that subject / session
        sub_id (_type_): New name for the subject
        sess_id (_type_): ID of the session to import
    """

    X = sio.loadmat(source_dir + '/design_matrix_unf.mat')
    DM = X['X']
    filename = dest_dir + f'/{sub_id}_{sess_id}_designmatrix_unf.npy'
    np.save(filename, DM)

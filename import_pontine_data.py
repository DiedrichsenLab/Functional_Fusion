import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio


def import_suit(source_dir,dest_dir,anat_name,participant_id):
    """ Imports a suit folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file (submitted to )
        participant_id (str): ID of participant
    """
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src=[]
    dest =[]
    src.append(f'/c1{anat_name}.nii')
    dest.append(f'/sub-{participant_id}_label-GMc_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/sub-{participant_id}_label-WMc_probseg.nii')
    src.append(f'/maskbrainSUITGrey.nii')
    dest.append(f'/sub-{participant_id}_desc-cereb_mask.nii')
    src.append(f'/u_a_c_{anat_name}_seg1.nii')
    dest.append(f'/sub-{participant_id}_space-SUIT_xfm.nii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except:
            print('skipping ' + src[i])


def import_anat(source_dir,dest_dir,anat_name,participant_id):
    """ Imports a anatomy folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory (anatomical)
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file
        participant_id (str): ID of participant
    """
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    src=[]
    dest =[]
    src.append(f'/{anat_name}.nii')
    dest.append(f'/sub-{participant_id}_T1w.nii')
    src.append(f'/c1{anat_name}.nii')
    dest.append(f'/sub-{participant_id}_label-GM_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/sub-{participant_id}_label-WM_probseg.nii')
    src.append(f'/c3{anat_name}.nii')
    dest.append(f'/sub-{participant_id}_label-CSF_probseg.nii')
    src.append(f'/u_a_c_{anat_name}_seg1.nii')
    dest.append(f'/sub-{participant_id}_space-MNI_xfm.nii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except:
            print('skipping ' + src[i])

def import_freesurfer(source_dir,dest_dir,new_id):
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src=[]
    dest =[]
    src.append(f'/lh.pial.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_pial.surf.gii')
    src.append(f'/lh.white.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_white.surf.gii')
    src.append(f'/rh.pial.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_pial.surf.gii')
    src.append(f'/rh.pial.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_white.surf.gii')
    src.append(f'/lh.sulc.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_sulc.shape.gii')
    src.append(f'/rh.sulc.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_sulc.shape.gii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except:
            print('skipping ' + src[i])

def import_spm_glm(source_dir,dest_dir,sub_id,sess_id,info_dict):
    """Imports the output of the SPM GLM with an SPM_info.mat
    structure into BIDS deriviatie (Functional Fusion) framework.
    It assumes that a single GLM corresponds to single session.

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that subject / session
        new_id (_type_): New name for the subject
        info_dict (_type_): Dictonary with the old field names and the new field names for the information
    """
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src=[]
    dest =[]
    # Generate new dictionary from SPM info
    D = mat73.loadmat(source_dir+'/SPM_info.mat')
    T={}
    for i in info_dict.items():
        series=D[i[0]]
        if type(series[0]) is list:
            series=[series[i][0] for i in range(len(series))]
        T[i[1]]=series

    N = len(T[i[1]])
    if 'reg_id' not in T.keys():
        n = sum(T['run']==1)
        T['reg_num'] = np.arange(N)
        T['reg_id'] = T['reg_num'] % n

    # Ensure that run number is an integer value
    T['run'] = [int(j) for j in T['run']]
    D = pd.DataFrame(T)
    D.to_csv(dest_dir + f'/{sub_id}_{sess_id}_reginfo.tsv',sep='\t')

    # Prepare beta files for transfer
    src=[]
    dest =[]
    for i in range(N):
        src.append(f'/beta_{i+1:04d}.nii')
        dest.append(f'/{sub_id}_{sess_id}_run-{D.run[i]:02}_reg-{D.reg_id[i]:02d}_beta.nii')
    # Mask
    src.append(f'/mask.nii')
    dest.append(f'/{sub_id}_{sess_id}_mask.nii')

    # ResMS
    src.append(f'/resms.nii')
    dest.append(f'/{sub_id}_{sess_id}_resms.nii')

    # Copy those files over
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except:
            print('skipping ' + src[i])

def import_spm_designmatrix(source_dir,dest_dir,sub_id,sess_id):
    """Imports the output of the SPM GLM with an SPM_info.mat
    structure into BIDS deriviatie (Functional Fusion) framework.
    It assumes that a single GLM corresponds to single session.

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that subject / session
        sub_id (_type_): New name for the subject
        sess_id (_type_): ID of the session to import
    """
    X = sio.loadmat(source_dir+'/design_matrix.mat')
    DM = X['X']
    filename = dest_dir + f'/{sub_id}_{sess_id}_designmatrix.npy'
    np.save(filename,DM)


if __name__ == '__main__':
    base_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Pontine7T/'
    dest_base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Pontine7T/'
    for participant_id in ['01', '03', '04', '07', '95', '96', '97', '98']:

        # # --- Importing SUIT ---
        # source_dir = '{}/suit/anatomicals/S{}'.format(base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/suit'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_suit(source_dir,dest_dir,anat_name,participant_id)

        # # --- Importing ANAT ---
        # source_dir = '{}/anatomicals/S{}'.format(base_dir, participant_id)
        # dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        # anat_name = 'anatomical'
        # import_anat(source_dir,dest_dir,anat_name,participant_id)

        # --- Importing Freesurfer ---
        source_dir = '{}/surfaceFreesurfer/S{}/surf'.format(base_dir, participant_id)
        dest_dir = '{}/derivatives/sub-{}/anat'.format(dest_base_dir, participant_id)
        new_id = 'sub-{}'.format(participant_id)
        import_freesurfer(source_dir,dest_dir,new_id)




# Missing:
# skipping /maskbrainSUITGrey.nii

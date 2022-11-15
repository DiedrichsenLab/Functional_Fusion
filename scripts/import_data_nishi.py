# from msilib import PID_CREATE_DTM
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
import os
import import_data as im


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
    dest.append(f'/{participant_id}_label-GMc_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/{participant_id}_label-WMc_probseg.nii')
    src.append(f'/maskbrainSUITGrey.nii')
    dest.append(f'/{participant_id}_desc-cereb_mask.nii')
    src.append(f'/y_{anat_name}_suitdef.nii')
    dest.append(f'/{participant_id}_space-SUIT_xfm.nii')
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
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except:
            print('skipping ' + src[i])


def import_spm_glm(source_dir,dest_dir,sub_id,sess_id):
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
    # D = mat73.loadmat(source_dir+'/SPM_info.mat')
    T={}
    # for i in info_dict.items():
    #     series=D[i[0]]
    #     if type(series[0]) is list:
    #         series=[series[i][0] for i in range(len(series))]
    #     T[i[1]]=series

    # N = len(T[i[1]])
    # if 'reg_id' not in T.keys():
    #     n = sum(T['run']==1)
    #     T['reg_num'] = np.arange(N)
    #     T['reg_id'] = T['reg_num'] % n

    # Ensure that run number is an integer value
    # T['run'] = [int(j) for j in T['run']]
    # D = pd.DataFrame(T)
    # D.to_csv(dest_dir + f'/{sub_id}_{sess_id}_reginfo.tsv',sep='\t')
    
    D = pd.read_csv(source_dir + f'/{sub_id}_{sess_id}_reginfo.tsv',sep='\t')

    # Prepare beta files for transfer
    src=[]
    dest =[]
    for i in range(len(D.index)):
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



if __name__ == "__main__":
    base_source_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Nishimoto/raw/'
    base_dest_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Nishimoto/derivatives/'
    
    subj_list = []
    for i in range(1, 7):
        subj_list.append(f"sub-{i:02d}")


    for sub in subj_list:
        
            # # --- Import GLM data ---
            # source_dir = os.path.join(base_source_dir, sub, 'estimates', 'glm01', f'ses-{ss:02d}')
            # dest_dir = os.path.join(base_dest_dir,  sub, 'estimates', f"ses-{ss:02d}")
            # import_spm_glm(source_dir,dest_dir,sub , f"ses-{ss:02d}")
            # # import_spm_designmatrix(source_dir,dest_dir,sub,f'ses-{ss:02d}')

            # --- Import suit data ---
            source_dir = os.path.join(
                base_source_dir, sub, 'suit', 'anat')
            dest_dir = os.path.join(
                base_dest_dir, sub, 'suit')
            import_suit(source_dir, dest_dir, f'{sub}_T1w_lpi', sub)

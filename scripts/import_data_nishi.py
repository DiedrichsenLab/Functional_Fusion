# from msilib import PID_CREATE_DTM
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
import os
import import_data as im



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

    T={}
    
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
    base_source_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/Nishimoto/raw/'
    base_dest_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Nishimoto/derivatives/'
    
    subj_list = []
    for i in range(1, 7):
        subj_list.append(f"sub-{i:02d}")


    for sub in subj_list:

        # --- Import Freesurfer data ---
        source_dir = '{}/surfaceWB/data/{}/'.format(
            base_source_dir, sub)
        dest_dir = '{}/{}/anat'.format(
            base_dest_dir, sub)
        old_id = 'S{}'.format(sub)
        new_id = 'sub-{}'.format(sub)
        im.import_freesurfer(source_dir, dest_dir, old_id=sub, new_id=sub)


        # for ss in [1, 2]:
            
        #     # # --- Import GLM data ---
        #     # source_dir = os.path.join(base_source_dir, sub, 'estimates', 'glm01', f'ses-{ss:02d}')
        #     # dest_dir = os.path.join(base_dest_dir,  sub, 'estimates', f"ses-{ss:02d}")
        #     # import_spm_glm(source_dir,dest_dir,sub , f"ses-{ss:02d}")
        #     # # import_spm_designmatrix(source_dir,dest_dir,sub,f'ses-{ss:02d}')

        #     # --- Import suit data ---
        #     source_dir = os.path.join(
        #         base_source_dir, sub, 'suit', 'anat')
        #     dest_dir = os.path.join(
        #         base_dest_dir, sub, 'suit')
        #     import_suit(source_dir, dest_dir, f'{sub}_T1w_lpi', sub)


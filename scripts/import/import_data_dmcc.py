import os
import glob
import re
import shutil

from pathlib import Path
from this import d

import numpy as np
import pandas as pd

import Functional_Fusion.import_data as im

import scipy.io as sio
import nibabel as nb

# code to transfer dmcc data into functional fusion framework
from pathlib import Path

from Functional_Fusion.import_data import *

raw_dir = Path('/srv/diedrichsen/data/Cerebellum/DMCC/derivatives/fmriprep-1.3.2') # folder where raw data is stored
dest_dir = Path('/srv/diedrichsen/data/FunctionalFusion/DMCC/') # folder in Functional Fusion
subjects = ['sub-f1027ao', 'sub-f1031ax', 'sub-f1342ku', 'sub-f1550bc',
            'sub-f1552xo', 'sub-f1659oa', 'sub-f1670rz', 'sub-f1828ko',
            'sub-f2157me', 'sub-f2499cq', 'sub-f2593wi', 'sub-f2648qw',
            'sub-f2709ul', 'sub-f2968sp', 'sub-f3300jh', 'sub-f3387yq',
            'sub-f3469wa', 'sub-f3526dz', 'sub-f3680fb', 'sub-f3681wf',
            'sub-f3996sp', 'sub-f4138ge', 'sub-f4310gw', 'sub-f4354bs',
            'sub-f4467ur', 'sub-f4796rs', 'sub-f4831tn', 'sub-f5001ob',
            'sub-f5004cr', 'sub-f5094na', 'sub-f5094ya', 'sub-f5386yx',
            'sub-f5416zj', 'sub-f5445nh', 'sub-f5635rv', 'sub-f5650zm',
            'sub-f5930vp', 'sub-f6188io', 'sub-f6318if', 'sub-f6464bf',
            'sub-f6950qp', 'sub-f7227ag', 'sub-f7688lh', 'sub-f7951pz',
            'sub-f8113do', 'sub-f8194sp', 'sub-f8270up', 'sub-f8294bu',
            'sub-f8298ds', 'sub-f8570ui', 'sub-f8710qa', 'sub-f8979ai',
            'sub-f9057kp', 'sub-f9206gd', 'sub-f9271ex']        

def wrapper_suit(subj_list):

    # loop over subjects
    for subj_id in subj_list:
        # get the directory for the subject with respect to the source
        # and destination directories
        src_dir_subj = f"{raw_dir}/{subj_id}/suit"
        dest_dir_subj = f"{dest_dir}/derivatives/{subj_id}/suit"

        # get the anatomical file name
        anat_name = f'{subj_id}_space-MNI152NLin2009cAsym_desc-preproc_T1w'

        # import the suit files
        im.import_suit(source_dir = src_dir_subj, 
                    dest_dir = dest_dir_subj, 
                    anat_name = anat_name, 
                    participant_id = subj_id)
    return

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
    src.append(f'/{anat_name}_desc-preproc_T1w.nii') #sub-f4796rs_space-MNI152NLin2009cAsym_desc-preproc_T1w
    dest.append(f'/{participant_id}_T1w.nii')
    src.append(f'/{anat_name}_label-GM_probseg.nii')
    dest.append(f'/{participant_id}_label-GM_probseg.nii')
    src.append(f'/{anat_name}_label-WM_probseg.nii')
    dest.append(f'/{participant_id}_label-WM_probseg.nii')
    src.append(f'/{anat_name}_label-CSF_probseg.nii')
    dest.append(f'/{participant_id}_label-CSF_probseg.nii')
    src.append(f'/y_T1w-MNI152NLin2009cAsym_xfm.nii')
    dest.append(f'/{participant_id}_space-MNI_xfm.nii')
    for i in range(len(src)):
        try:
            print(f"** copying xfm over")
            shutil.copyfile(source_dir + src[i], dest_dir + dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])

def wrapper_anat(subj_list):

    # loop over subjects
    for subj_id in subj_list:
        print(f"- importing {subj_id}")
        # get the directory for the subject with respect to the source
        # and destination directories
        src_dir_subj = f"{raw_dir}/{subj_id}/anat"
        dest_dir_subj = f"{dest_dir}/derivatives/{subj_id}/anat"

        # get the anatomical file name
        anat_name = f'{subj_id}_space-MNI152NLin2009cAsym'

        # import the suit files
        import_anat(source_dir = src_dir_subj, 
                    dest_dir = dest_dir_subj, 
                    anat_name = anat_name, 
                    participant_id = subj_id)
    return

def wrapper_freesurfer(subj_list):
    raw_dir = Path('/srv/diedrichsen/data/Cerebellum/DMCC/derivatives/fmriprep-1.3.2') # folder where raw data is stored
    dest_dir = Path('/srv/diedrichsen/data/FunctionalFusion/DMCC/') # folder in Functional Fusion


    for subj_id in subj_list:
        src_dir_subj = f"{raw_dir}/surfaceWB/data/{subj_id}/"
        dest_dir_subj = f"{dest_dir}/derivatives/{subj_id}/anat"
        old_id_subj = subj_id
        new_id_subj = subj_id
        im.import_freesurfer(src_dir_subj, dest_dir_subj, old_id_subj, new_id_subj)
    return

def import_spm_glm(source_dir,dest_dir,sub_id,sess_id, task_name):
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

    D = pd.read_csv(source_dir + f'/{sub_id}_ses-wave1bas_{task_name}_reginfo.tsv',sep='\t')

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
    src.append(f'/ResMS.nii')
    dest.append(f'/{sub_id}_{sess_id}_resms.nii')

    # Copy those files over
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except:
            print('skipping ' + src[i])

def wrapper_spm_glm(subj_list):
    # betas are stored in separate folders for each task
    # in the functional fusion, I have two options:
    # 1. store all betas in one folder and use the spm_info file to identify the task
    # 2. store betas in separate folders for each task, in other words having separate sessions per task
    # option 2 would be easier down the line when I'm writing the dataset class for dmcc
    # For selective recruitment I can loop through tasks (sessions in option 2) and make up a big dataframe
    # and then do the final regression
    # option 2 will have a 'bas' in the name of the session as well, since I might add data
    # from other sessions (pro and rea) later on: ses-axcpt_bas
    # in option 1 all the files will be moved to a folder called ses-bas
    # I will go with option 2 for now

    # list of tasks to loop over
    task_list = ['Axcpt', 'Cuedts', 'Stern', 'Stroop']

    # first get the source dir
    # loop over subjects
    for subj_id in subj_list:

        # loop over tasks and copy them over
        for task in task_list:
            # get the source directory for the subj_id
            src_dir_subj = f"{raw_dir}/{subj_id}/ses-wave1bas/estimates/glm01/{task}"
            # get the destination directory for the subj_id
            dest_dir_subj = f"{dest_dir}/derivatives/{subj_id}/estimates/ses-{task.lower()}_bas"

            import_spm_glm(source_dir = src_dir_subj, 
                           dest_dir = dest_dir_subj, 
                           sub_id = subj_id, 
                           sess_id =  f"ses-{task.lower()}_bas",
                           task_name = task)
            
            # import the design matrix you have saved 
            im.import_spm_designmatrix(src_dir_subj, dest_dir_subj, subj_id, f"ses-{task.lower()}_bas")
    return

def wrapper_info_tsv(subj_list):
    # list of tasks to loop over
    task_list = ['Axcpt', 'Cuedts', 'Stern', 'Stroop']

    # first get the source dir
    # loop over subjects
    for subj_id in subj_list:
        for task in task_list:
            # get the source directory for the subj_id
            src_dir_subj = f"{raw_dir}/{subj_id}/ses-wave1bas/estimates/glm01/{task}"
            # get the destination directory for the subj_id
            dest_dir_subj = f"{dest_dir}/derivatives/{subj_id}/estimates/ses-{task.lower()}_bas"

            filename_src = f"{src_dir_subj}/{subj_id}_ses-wave1bas_{task}_reginfo.tsv"
            filename_dest = f"{dest_dir_subj}/{subj_id}_ses-{task.lower()}_bas_reginfo.tsv"

            # Copy those files over
            try:
                shutil.copyfile(filename_src,filename_dest)
            except:
                print('skipping ' + filename_src)

    return

if __name__ == "__main__":
    wrapper_suit(subj_list = subjects[0:20])
    wrapper_anat(subj_list = subjects[0:20])
    wrapper_freesurfer(subj_list = subjects[0:20])
    wrapper_spm_glm(subj_list = subjects[0:20])
    wrapper_info_tsv(subj_list = subjects[0:20])
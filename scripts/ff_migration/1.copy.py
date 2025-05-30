import numpy as np
import os 
import shutil
import pandas as pd
def copy_dataset(source_dir, dest_dir, dataset, subject_list=None):
    """
    Copy specific subjects or entire dataset from source to destination.

    Parameters:
    - source_dir (str): Path to FunctionalFusion
    - dest_dir (str): Path to FunctionalFusion_new
    - dataset (str): Dataset name (e.g. 'HCP_tfMRI')
    - subject_list (list or None): List of subject folder names to copy (e.g. ['sub-101309', 'group'])
                                   If None or empty, copy full dataset like it is
    """
    src_path = os.path.join(source_dir, dataset)
    dst_path = os.path.join(dest_dir, dataset)

    # Copy participant.tsv file 
    os.makedirs(dst_path, exist_ok=True)
    tsv_src = os.path.join(src_path, 'participants.tsv')
    tsv_dst = os.path.join(dst_path, 'participants.tsv')
    shutil.copy2(tsv_src, tsv_dst)
    
    if not subject_list:
        T = pd.read_csv(tsv_src, sep='\t')
        subject_list = T['participant_id'].tolist() 
    
    for subj in subject_list:
        # Copy files in estimates folder
        src = os.path.join(src_path, 'derivatives', subj,'estimates')
        dst = os.path.join(dst_path, 'derivatives', 'ffimport',subj,'func')
        save_copy(src, dst)
        # Copy files in anatomical folder
        src = os.path.join(src_path, 'derivatives', subj,'anat')
        dst = os.path.join(dst_path, 'derivatives', 'ffimport',subj,'anat')
        save_copy(src, dst)
        # Copy files in anatomical folder
        src = os.path.join(src_path, 'derivatives', subj,'suit')
        dst = os.path.join(dst_path, 'derivatives', 'ffimport',subj,'anat')
        save_copy(src, dst)
        # Copy files in data folder to ffextract
        src = os.path.join(src_path, 'derivatives', subj,'data')
        dst = os.path.join(dst_path, 'derivatives', 'ffextract',subj)
        save_copy(src, dst)
    # Copy files in group folder
    src = os.path.join(src_path, 'derivatives', 'group','data')
    dst = os.path.join(dst_path, 'derivatives', 'ffextract','group')
    save_copy(src, dst)


def save_copy(src, dst):
    """
    Save a copy of the source directory to the destination directory.
    
    Parameters:
    - src (str): Source directory path
    - dst (str): Destination directory path
    """
    if os.path.exists(src):
        shutil.copytree(src, dst, dirs_exist_ok=True)
        print(f"Copied {src}")
    else:
        print(f"{src} does not exist")


if __name__=='__main__':
    # dirs
    datashare_dir = 'Y:/data'
    if not os.path.exists(datashare_dir):
        datashare_dir = '/cifs/diedrichsen/data'
    if not os.path.exists(datashare_dir):
        datashare_dir = '/Volumes/diedrichsen_data$/data'

    source_dir = f'{datashare_dir}/FunctionalFusion/'
    dest_dir = f'{datashare_dir}/FunctionalFusion_new/'

    # what to copy
    dataset = 'Language'
    subject_list = ['sub-103111']

    # copy files
    copy_dataset(source_dir, dest_dir, dataset, subject_list=None)


    pass
import numpy as np
import os 
import shutil

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

    if not subject_list:
        shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
        print(f"Entire dataset '{dataset}' copied.")
    else:
        os.makedirs(os.path.join(dst_path, 'derivatives'), exist_ok=True)
        for subj in subject_list:
            src = os.path.join(src_path, 'derivatives', subj)
            dst = os.path.join(dst_path, 'derivatives', subj)
            if os.path.exists(src):
                shutil.copytree(src, dst, dirs_exist_ok=True)
                print(f"Copied {subj}")
            else:
                print(f"{subj} not found in derivatives source")

        # copy participant.tsv 
        tsv_src = os.path.join(src_path, 'participants.tsv')
        tsv_dst = os.path.join(dst_path, 'participants.tsv')
        shutil.copy2(tsv_src, tsv_dst)
        print(f"Copied participants.tsv for {dataset}")



if __name__=='__main__':
    # dirs
    datashare_dir = 'Y:/data/'
    if not os.path.exists(datashare_dir):
        datashare_dir = '/cifs/diedrichsen/data/'

    source_dir = f'{datashare_dir}/FunctionalFusion/'
    dest_dir = f'{datashare_dir}/FunctionalFusion_new/'

    # what to copy
    dataset = 'HCP_tfMRI'
    subject_list = ['sub-101309']

    # copy files
    copy_dataset(source_dir, dest_dir, dataset, subject_list)


    pass
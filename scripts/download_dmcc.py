import numpy as np
import pandas as pd
import openneuro as on

import os

from pathlib import Path


target_dir = Path("/srv/diedrichsen/data/Cerebellum/DMCC")
data_onid = "ds003465"

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

# subjects = ['sub-f1027ao', 'sub-f1031ax', 'sub-f1342ku', 'sub-f1550bc',
#             'sub-f1552xo', 'sub-f1659oa', 'sub-f1670rz', 'sub-f1828ko']

def get_download_list(subj_id = subjects):
    """
    returns a list of anatomical, functional, event file, 
    and freesurfer files to be downloaded
    Args:


    """
    # the folowing is a list of anatomicals (native and MNI152NLin2009cAsym) to be downloaded from the derivatives folder
    anatomicals = ["_desc-preproc_T1w",\
                   "_desc-brain_mask",\
                    "_label-CSF_probseg",\
                        "_label-GM_probseg",\
                            "_label-WM_probseg",\
                                "_space-MNI152NLin2009cAsym_desc-preproc_T1w",\
                                    "_space-MNI152NLin2009cAsym_desc-brain_mask",\
                                        "_space-MNI152NLin2009cAsym_label-CSF_probseg",\
                                            "_space-MNI152NLin2009cAsym_label-GM_probseg",\
                                                "_space-MNI152NLin2009cAsym_label-WM_probseg"]
    
    # make a list of anatomicals to be downloaded
    anat_files = []
    for subj in subj_id:
        for anat in anatomicals:
            anat_files.append(f"derivatives/fmriprep-1.3.2/{subj}/anat/{subj}{anat}.nii.gz")
    

    # get the list of freesurfer files
    surf_list = ["inflated", "midthickness", "pial", "smoothwm"]
    hemis = ["L", "R"]
    fs_files = []
    for subj in subj_id:
        print(subj)
        for hemi in hemis:
            for surf in surf_list:
                fs_files.append(f"derivatives/fmriprep-1.3.2/{subj}/anat/{subj}_hemi-{hemi}_{surf}.surf.gii")


    # the folowing is a list of functionals (native and MNI152NLin2009cAsym) to be downloaded from the derivatives folder
    # NOTE that functionals are under ses-wave1bas folder
    # example:  sub-f1027ao_ses-wave1bas_task-Axcpt_acq-mb4AP_run-1_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz
    # sub-f1027ao_ses-wave1bas_task-Axcpt_acq-mb4AP_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
    tasks = ["Axcpt", "Cuedts", "Stern", "Stroop"]
    run_enc = ["mb4AP_run-1", "mb4PA_run-2"]
    bolds = ["space-MNI152NLin2009cAsym_desc-preproc_bold", \
             "space-MNI152NLin2009cAsym_desc-brain_mask", \
             "space-MNI152NLin2009cAsym_desc-aparcaseg_dseg"]
    func_files = []
    event_files = []
    for subj in subj_id:
        for task in tasks:
            for run in run_enc:
                event_files.append(f"{subj}/ses-wave1bas/func/{subj}_ses-wave1bas_task-{task}_acq-{run}_events.tsv")
                for bold in bolds:
                    func_files.append(f"derivatives/fmriprep-1.3.2/{subj}/ses-wave1bas/func/{subj}_ses-wave1bas_task-{task}_acq-{run}_{bold}.nii.gz")
    
    return anat_files, fs_files, func_files, event_files

def download_filelist(download_dir = target_dir, download_include = subjects):
    """
    takes in a list of files to be downloaded and download them into the directory specified.
    if the list is just a list of subject names, the raw data will be downloaded for each subject in the list
    if you want to download specific files, you need to create a list contining full path to the data and pass it on
    """
    # download the subjects passed on in subj_id
    on.download(dataset=data_onid, target_dir=download_dir, include=download_include)

def move_event_tsv():
    """
    in case you have downloaded event files to the raw directory instead of derivatives
    """
    tasks = ["Axcpt", "Cuedts", "Stern", "Stroop"]
    run_enc = ["mb4AP_run-1", "mb4PA_run-2"]

    source_dir = target_dir
    dest_dir = target_dir / "derivatives" / "fmriprep-1.3.2"
    for subj in subj_id:
        for task in tasks:
            for run in run_enc:
                event_name = f"{subj}/ses-wave1bas/func/{subj}_ses-wave1bas_task-{task}_acq-{run}_events.tsv"
                os.rename(f"{source_dir}/{subj}/ses-wave1bas/func/{event_name}", f"{dest_dir}/{subj}/ses-wave1bas/func/{event_name}")                
    return


if __name__ == "__main__":
    anat_list, func_list = get_download_list(subj_id = subjects[0:2])
    download_filelist(download_dir = target_dir, download_include = anat_list)
    pass






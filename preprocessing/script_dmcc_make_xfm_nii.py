# code to generate deformation field nifti image from fmriprep outputs
# This was specifically written for the DMCC data but can be used with other data coming out
# of fmriprep with proper naming based on the atlas space
# see this link for more details on how to get xfm nifti from .h5 transformation saved in fmriprep:
# https://neurostars.org/t/h5-to-affine-warpfield-nifti/7276

import numpy as np
import pandas as pd
import os
import subprocess

from pathlib import Path

base_dir = Path('/srv/diedrichsen/data/Cerebellum/DMCC/derivatives/fmriprep-1.3.2') # folder where raw data is stored

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

def dissemble_xfm_h5(subj_list = subjects, prefix = 'y'):
    """
    this is written for when you have outputs from fmriprep
    """
    for subj_id in subj_list:
        print(f"- dissembling {subj_id}")
        # cd to the subject directory
        os.chdir(f"{base_dir}/{subj_id}/anat/")
        # make a string for the command
        # sub-f1027ao_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
        # xfm_filename = f"{subj_id}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
        xfm_filename = f"{subj_id}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"

        cmd_dissemble = f"CompositeTransformUtil --disassemble {xfm_filename} {prefix}"
        
        # run the command using subprocess
        # the command will append the prefix to the dissmebled .mat and .nii.gz files
        # and saves it in the same directory as xfm_filename
        subprocess.check_output(cmd_dissemble, shell=True)

        # some cleaning up:
        # running this command, two files are created:
        # 00_y_AffineTransform.mat
        # 01_y_DisplacementFieldTransform.nii.gz
        # rename and unzip the nii.gz file
        old_name = f"01_y_DisplacementFieldTransform.nii"
        # old_name = f"00_y_DisplacementFieldTransform.nii"
        # first gunzip it
        cmd_gzip = f"gunzip {old_name}.gz"
        # run the gunzip command
        subprocess.check_output(cmd_gzip, shell = True)

        # make the mv command to rename the file
        new_name = f"y_T1w-MNI152NLin2009cAsym_xfm.nii"
        cmd_rename = f"mv {old_name} {new_name}"
        subprocess.check_output(cmd_rename, shell=True)

    return





# Script for importing the MDTB data set from super_cerebellum to general format.
import shutil
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
import nibabel as nb
import nitools as nt
import sys

ERIS_DIR = '/home/dzhi/eris_mount'
if not Path(ERIS_DIR).exists():
    ERIS_DIR = '/data/tge'
if not Path(ERIS_DIR).exists():
    raise (NameError('Could not find ERIS_DIR'))

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = ERIS_DIR + '/Tian/UKBB_full/imaging'
if not Path(base_dir).exists():
    print('data server not mounted')

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

def resample_atlases(atlas_names):
    # create and calculate the atlas map for each participant
    for at in atlas_names:
        print(f'at')
        xfm_name = atlas_dir + '/tpl-MNI152NLin6AsymC/tpl-MNI152NLin6AsymC_from-SUIT_mode-image_xfm.nii'
        deform = nb.load(xfm_name)
        source_name = atlas_dir + f'/tpl-SUIT/atl-{at}_space-SUIT_dseg.nii'
        source = nb.load(source_name)
        nifti = nt.deform_image(source,deform,0)
        out_name = atlas_dir + f'/tpl-MNI152NLin6AsymC/atl-{at}_space-MNI152NLin6AsymC_dseg.nii'
        nb.save(nifti,out_name)

if __name__ == "__main__":
    resample_atlases(['Anatom','Buckner7','Buckner17','Ji10','MDTB10'])
    # make_mdtb_suit()
    pass


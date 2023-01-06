# Script for testing parcel axis 
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetMDTB
from dataset import DataSetHcpResting
import nibabel as nb
from matrix import indicator
import sys
import os
import sys
sys.path.append("~/Documents/Projects/Functional_Fusion")

base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if sys.platform == "win32":  # for windows user
    base_dir = 'Y:/data/FunctionalFusion'

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'


def test_mdtb_suit():
    
    # create an instance of SUIT3 volumetric atlas
    mask_image = atlas_dir + '/tpl-SUIT' + '/tpl-SUIT_res-3_gmcmask.nii'
    name = 'SUIT3'
    SUIT3_atlas = am.AtlasVolumetric(name, mask_img=mask_image, structure="cerebellum")

    # create parcel axis for mdtb10 parcellation
    label_img = atlas_dir + '/tpl-SUIT' + '/atl-MDTB10_space-SUIT_dseg.nii'
    SUIT3_atlas.get_parcel(label_img=label_img)

    return SUIT3_atlas

def test_mdtb_fs32k():
    
    # create an instance of SUIT3 volumetric atlas
    mask_image = []
    for hemi in ['L', 'R']:
        mask_image.append(atlas_dir + '/tpl-fs32k' + f'/tpl-fs32k_hemi-{hemi}_mask.label.gii')
        name = 'fs32k'
    fs_atlas = am.AtlasSurface(name, mask_img=mask_image, structure=["cortex_left", 'cortex_right'])

    # create parcel axis for mdtb10 parcellation
    label_img = []
    for hemi in ['L', 'R']:
        label_img.append(atlas_dir + '/tpl-fs32k' + f'/Icosahedron-42_Sym.32k.{hemi}.label.gii')
    
    fs_atlas.get_parcel(label_img=label_img)

    return fs_atlas


if __name__ == "__main__":
    a = test_mdtb_suit()
    a.get_parcel_axis()

    b = test_mdtb_fs32k()
    b.get_parcel_axis()

    print("hovering")
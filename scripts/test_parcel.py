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
import os
import sys
from nibabel import cifti2
sys.path.append("~/Documents/Projects/Functional_Fusion")

base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if sys.platform == "win32":  # for windows user
    base_dir = 'Y:/data/FunctionalFusion'

project_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'


def test_mdtb_suit():
    
    # create an instance of SUIT3 volumetric atlas
    mask_image = atlas_dir + '/tpl-SUIT' + '/tpl-SUIT_res-3_gmcmask.nii'
    name = 'SUIT3'
    SUIT3_atlas = am.AtlasVolumetric(name, mask_img=mask_image, structure="cerebellum")

    # create parcel axis for mdtb10 parcellation
    label_img = atlas_dir + '/tpl-SUIT' + '/atl-MDTB10_space-SUIT_dseg.nii'
    SUIT3_atlas.get_parcel(label_img=label_img)

    SUIT3_atlas.get_parcel_axis()

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
    fs_atlas.get_parcel_axis()

    return fs_atlas

# extracting average within a parcel:
# 1. load the extracted data DONE
# 2. use the label image to create a parcel axis for the space DONE
# 3. calculate mean for each parcel 
# 4. save! 

def extract_parcel_data_suit():
    """
    Example code to extract mean across parcel for mdtb suit
    """
    # perform all the steps in test_mdtb_suit() to create a parcel axis for the parcellation
    atlas_parcel = test_mdtb_suit()

    # load in the data for sub-02
    data_dir = project_dir + '/derivatives' + '/sub-02' +  '/data' + '/sub-02_space-SUIT3_ses-s1_CondAll.dscalar.nii'
    data_img = nb.load(data_dir)
    # get data for the structure
    data = data_img.get_fdata()

    # load the tsv file
    info_file = project_dir + '/derivatives' + '/sub-02' +  '/data' + '/sub-02_ses-s1_info-CondAll.tsv'
    info_tsv = pd.read_csv(info_file, sep = '\t')

    # create a matrix for aggregating data (cannot use dataset.agg_data now! Need to make changes)
    C = indicator(atlas_parcel.label_vector[0],positive=True)

    # get the mean across parcel
    data = np.nan_to_num(data)
    data_parcel = (data @ C)/np.sum(C, axis = 0)

    # create a dscale cifti with parcelAxis labels and data_parcel
    row_axis = info_tsv.cond_name
    row_axis = nb.cifti2.ScalarAxis(row_axis)
        

    bm = atlas_parcel.get_brain_model_axis()
    pa = atlas_parcel.parcel_axis
    # HEAD = cifti2.Cifti2Header.from_axes((row_axis,bm,pa))
    header = nb.Cifti2Header.from_axes((row_axis, pa))
    cifti_img = nb.Cifti2Image(dataobj=data_parcel, header=header)
    nb.save(cifti_img, 'test_mdtb_10.pscalar.nii')
    A = nb.load('test_mdtb_10.pscalar.nii')
    return data_parcel

def extract_parcel_data():
    """
    Example code to extract mean across parcel for tessellation
    """
    # perform all the steps in test_mdtb_suit() to create a parcel axis for the parcellation
    atlas_parcel = test_mdtb_fs32k()

    # load in the data for sub-02
    data_dir = project_dir + '/derivatives' + '/sub-02' +  '/data' + '/sub-02_space-fs32k_ses-s1_CondAll.dscalar.nii'
    data_img = nb.load(data_dir)
    # get data for the structure
    data = data_img.get_fdata()

    # load the tsv file
    info_file = project_dir + '/derivatives' + '/sub-02' +  '/data' + '/sub-02_ses-s1_info-CondAll.tsv'
    info_tsv = pd.read_csv(info_file, sep = '\t')

    if len(atlas_parcel.label_vector) == 1: # for cerebellum (not divided by hemi)
        vector = atlas_parcel.label_vector[0]
    else: # for cortex, divided into left and right
        vector = np.concatenate(atlas_parcel.label_vector, axis = 0)

    # create a matrix for aggregating data (cannot use dataset.agg_data now! Need to make changes)
    C = indicator(vector,positive=True)

    # get the mean across parcel
    data = np.nan_to_num(data)
    data_parcel = (data @ C)/np.sum(C, axis = 0)

    # create a dscale cifti with parcelAxis labels and data_parcel
    row_axis = info_tsv.cond_name
    row_axis = nb.cifti2.ScalarAxis(row_axis)
        

    bm = atlas_parcel.get_brain_model_axis()
    pa = atlas_parcel.parcel_axis
    # HEAD = cifti2.Cifti2Header.from_axes((row_axis,bm,pa))
    header = nb.Cifti2Header.from_axes((row_axis, pa))
    cifti_img = nb.Cifti2Image(dataobj=data_parcel, header=header)
    nb.save(cifti_img, 'test_mdtb_10.pscalar.nii')
    A = nb.load('test_mdtb_10.pscalar.nii')
    return data_parcel

if __name__ == "__main__":
    # a = test_mdtb_suit()
    # b = test_mdtb_fs32k()
    # extract_parcel_data_suit()
    extract_parcel_data()
    print("hovering")
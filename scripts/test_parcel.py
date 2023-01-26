# Script for testing parcel axis
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetMDTB
from Functional_Fusion.dataset import DataSetHcpResting
import nibabel as nb
from Functional_Fusion.matrix import indicator
import os
import sys
from nibabel import cifti2

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:\data\FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/Users/callithrix/Documents/Projects/Functional_Fusion/'
if not Path(base_dir).exists():
    base_dir = '/Users/jdiedrichsen/Data/FunctionalFusion/'
if not Path(base_dir).exists():
    raise (NameError('Could not find base_dir'))

project_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

def test_mdtb_suit():

    # create an instance of SUIT3 volumetric atlas
    SUIT3_atlas,_ = am.get_atlas('SUIT3',atlas_dir)

    # create parcel axis for mdtb10 parcellation
    label_img = atlas_dir + '/tpl-SUIT' + '/atl-MDTB10_space-SUIT_dseg.nii'
    SUIT3_atlas.get_parcel(label_img=label_img)

    pa = SUIT3_atlas.get_parcel_axis()

    return SUIT3_atlas

def test_mdtb_fs32k():

    # create an instance of SUIT3 volumetric atlas
    fs_atlas,_ = am.get_atlas('fs32k',atlas_dir)

    # create parcel axis for mdtb10 parcellation
    label_img = []
    for hemi in ['L', 'R']:
        label_img.append(atlas_dir + '/tpl-fs32k' + f'/Icosahedron-42_Sym.32k.{hemi}.label.gii')

    fs_atlas.get_parcel(label_img=label_img,unite_struct=False)
    pa = fs_atlas.get_parcel_axis()

    return fs_atlas


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
    a = test_mdtb_suit()
    b = test_mdtb_fs32k()
    # extract_parcel_data_suit()
    # extract_parcel_data()
    print("hovering")
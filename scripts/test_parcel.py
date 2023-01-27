# Script for testing parcel axis
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
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

def test_parcel_data(atlas):
    """
    Example code to extract mean across parcel for mdtb suit
    """
    # load in the data for sub-02
    mdtb = ds.get_dataset_class(base_dir,'MDTB')
    data, info = mdtb.get_data(space=atlas.name,ses_id='ses-s1',subj=[0],type="CondHalf")

    # create a matrix for aggregating data (cannot use dataset.agg_data now! Need to make changes)
    data_parcel = ds.agg_parcels(data[0],atlas.label_vector)

    # create a dscale cifti with parcelAxis labels and data_parcel
    row_axis = nb.cifti2.ScalarAxis(info.cond_name)
    pa = atlas.get_parcel_axis()

    # HEAD = cifti2.Cifti2Header.from_axes((row_axis,bm,pa))
    header = nb.Cifti2Header.from_axes((row_axis, pa))
    cifti_img = nb.Cifti2Image(dataobj=data_parcel, header=header)
    nb.save(cifti_img, f'test_{atlas.name}.pscalar.nii')


if __name__ == "__main__":
    a = test_mdtb_suit()
    test_parcel_data(a)
    b = test_mdtb_fs32k()
    test_parcel_data(b)

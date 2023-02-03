# Script for testing parcel axis
import pandas as pd
from pathlib import Path
import numpy as np
import sys
sys.path.append('..')
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import nibabel as nb
import os


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
conn_dir = '/srv/diedrichsen/data/Cerebellum/connectivity/MDTB/train'

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

def plot_weights(method = "L2Regression", 
                 cortex = "Icosahedron-1002_Sym.32k", 
                 cerebellum = "NettekovenSym34", 
                 log_alpha = 8, 
                 dataset_name = "MDTB", 
                 ses_id = "ses-s1", 
                 ):
    # get the file containing best weights
    filename = os.path.join('/srv/diedrichsen/data/Cerebellum/connectivity/MDTB/train', f'{cortex}_{ses_id}_{method}_logalpha_{log_alpha}_best_weights.npy')
    weights = np.load(filename)

    # get atlases and create parcels/parcel labels
    atlas_cereb, ainfo = am.get_atlas('SUIT3',atlas_dir)
    atlas_cortex, ainfo = am.get_atlas('fs32k', atlas_dir)

    # get label files for cerebellum and cortex
    # NOTE: To average over cerebellum or cortex, pass on masks as label files
    label_cereb = atlas_dir + '/tpl-SUIT' + f'/atl-{cerebellum}_space-SUIT_dseg.nii'
    label_cortex = []
    for hemi in ['L', 'R']:
        label_cortex.append(atlas_dir + '/tpl-fs32k' + f'/{cortex}.{hemi}.label.gii')

    # get parcel for both atlases
    atlas_cereb.get_parcel(label_cereb)
    atlas_cortex.get_parcel(label_cortex, unite_struct = False)
    pa = atlas_cereb.get_parcel_axis()

    # get the maps for left and right hemispheres
    surf_map = []
    for label in atlas_cortex.label_list:
        # loop over regions within the hemisphere
        label_arr = np.zeros([atlas_cereb.n_labels, label.shape[0]])
        for p in np.arange(1, atlas_cereb.n_labels+1):
            vox = atlas_cereb.label_vector == p
            # get average connectivity weights
            weight_region = np.nanmean(weights[vox, :], axis = 0)
            for i in np.unique(label):            
                np.put_along_axis(label_arr[p-1, :], np.where(label==i)[0], weight_region[i-1], axis=0)
        surf_map.append(label_arr)

    cifti_img = atlas_cortex.data_to_cifti(surf_map, row_axis=pa.name)

    # save weight map
    nb.save(cifti_img,f'./Test.dscalar.nii')
    return
if __name__ == "__main__":
    plot_weights()
    # a = test_mdtb_suit()
    # test_parcel_data(a)
    # b = test_mdtb_fs32k()
    # test_parcel_data(b)

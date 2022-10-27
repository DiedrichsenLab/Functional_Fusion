# Script for getting all the HCP data for cerebellar-cortical connectivity
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetHcpResting
import nibabel as nb
import SUITPy as suit
import os
import sys

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

hcp_dir = base_dir + '/HCP'
atlas_dir = base_dir + '/Atlases'
hem_name = ['cortex_left','cortex_right']

def extract_hcp_suit(ses_id='ses-s1', type='CondHalf', atlas='MNISymC3'):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    hcp_dataset.extract_all_suit(ses_id,type,atlas)

def extract_hcp_data(res=162):
    # Make the atlas object
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)

    # initialize the data set object
    hcp_dataset = DataSetHcpResting(hcp_dir)

    # Get the deformation map from MNI to SUIT
    mni_atlas = atlas_dir + '/tpl-MNI152NLin6AsymC'
    deform = mni_atlas + '/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
    mask = mni_atlas + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
    atlas_map = am.AtlasMapDeform(hcp_dataset, suit_atlas, 'group',deform,mask)
    atlas_map.build(smooth=2.0)

    # Get the parcelation
    surf_parcel =[]
    for i,h in enumerate(['L','R']):
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i],gifti))

    T = hcp_dataset.get_participants()
    for s in T.participant_id:
        print(f'Extract {s}')
        coef = hcp_dataset.get_cereb_connectivity(s,atlas_map, surf_parcel)
        # Average across runs
        coef = np.nanmean(coef,axis=0)

        # Build a connectivity CIFTI-file and save
        bmc = suit_atlas.get_brain_model_axis()
        bpa = surf_parcel[0].get_parcel_axis() + surf_parcel[1].get_parcel_axis()
        header = nb.Cifti2Header.from_axes((bpa,bmc))
        cifti_img = nb.Cifti2Image(dataobj=coef,header=header)
        dest_dir = hcp_dataset.data_dir.format(s)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        nb.save(cifti_img, dest_dir + f'/sub-{s}_tessel-{res}.dpconn.nii')


def get_hcp_fs32k(res=162, index=range(0, 100), refix=False):
    # Make the atlas object
    mask_L = atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii'
    mask_R = atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii'
    fs32k_L_atlas = am.AtlasSurface('CORTEX_LEFT', mask_gii=mask_L)
    fs32k_R_atlas = am.AtlasSurface('CORTEX_RIGHT', mask_gii=mask_R)

    # initialize the data set object
    hcp_dataset = DataSetHcpResting(hcp_dir)

    # Get the parcelation
    surf_parcel = []
    for i, h in enumerate(['L', 'R']):
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i], gifti))

    T = hcp_dataset.get_participants()
    for i,s in enumerate(T.participant_id):
        data_dir = hcp_dataset.data_dir.format(s)
        Ci = nb.load(data_dir + f'/sub-{s}_tessel-{res}.dpconn.nii')
        if i==0:
            R = np.empty((T.shape[0],Ci.shape[0],Ci.shape[1]))

        R[i,:,:]=np.asanyarray(Ci.dataobj)
    Rm = np.nanmean(R,axis=0)
    Cm = nb.Cifti2Image(Rm,Ci.header)
    if refix:
        nb.save(Cm, hcp_dir + f'/group_tessel-{res}-ReFIX.dpconn.nii')
    else:
        nb.save(Cm, hcp_dir + f'/group_tessel-{res}.dpconn.nii')

    pass

def avrg_hcp_dpconn_cortex(res=162, index=range(0,100), refix=False):
    # initialize the data set object
    hcp_dataset = DataSetHcpResting(hcp_dir)
    T = hcp_dataset.get_participants()
    subjects_id = T.participant_id[index]

    for h in range(2):
        for i,s in enumerate(subjects_id):
            data_dir = hcp_dataset.data_dir.format(s)
            if refix:
                Ci = nb.load(data_dir + f'/sub-{s}_tessel-{res}_{hem_name[h]}-ReFIX.dpconn.nii')
            else:
                Ci = nb.load(data_dir + f'/sub-{s}_tessel-{res}_{hem_name[h]}.dpconn.nii')

            if i==0:
                R = np.empty((subjects_id.shape[0],Ci.shape[0],Ci.shape[1]))
            R[i,:,:]=np.asanyarray(Ci.dataobj)
        Rm = np.nanmean(R,axis=0)
        Cm = nb.Cifti2Image(Rm,Ci.header)
        if refix:
            nb.save(Cm, hcp_dir + f'/group_tessel-{res}_{hem_name[h]}-ReFIX.dpconn.nii')
        else:
            nb.save(Cm, hcp_dir + f'/group_tessel-{res}_{hem_name[h]}.dpconn.nii')

    pass

def parcel_hcp_dpconn(dpconn_file):
    """
    Args:
        dpconn_file (_type_): _description_
    """
    C = nb.load(dpconn_file)
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    label = atlas_dir + '/tpl-SUIT/atl-MDTB10_space-SUIT_dseg.nii'
    mdtb_atlas = am.AtlasVolumeParcel('MDTB10',label_img=label,mask_img=mask)
    R = mdtb_atlas.agg_data(np.asanyarray(C.dataobj))
    bm_cortex = C.header.get_axis(0)
    names = [f'MDTB {r+1:02}' for r in range(R.shape[1])]
    row_axis = nb.cifti2.ScalarAxis(names)
    cifti_img = nb.Cifti2Image(R.T,[row_axis,bm_cortex])
    return cifti_img

def parcel_hcp_dpconn_cortex(dpconn_file):
    """
    Args:
        dpconn_file (_type_): _description_
    """
    label_file = [atlas_dir + '/tpl-fs32k/ROI.32k.L.label.gii',
                  atlas_dir + '/tpl-fs32k/ROI.32k.R.label.gii',]
    cifti_img = []
    mask = [atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii',
            atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii']

    for h in range(2):
        C = nb.load(dpconn_file[h])
        roi_atlas = am.AtlasSurfaceParcel('ROI',label_gii=label_file[1-h], mask_gii=mask[h])
        R = roi_atlas.agg_data(np.asanyarray(C.dataobj))
        bm_cortex = C.header.get_axis(0)
        names = [name.label for name in roi_atlas.label_gii.labeltable.labels[1:]]
        row_axis = nb.cifti2.ScalarAxis(names)
        this_cifti_img = nb.Cifti2Image(R.T,[row_axis,bm_cortex])
        cifti_img.append(this_cifti_img)

    return cifti_img

def indv_hcp_pscalar(res=162, index=range(0,100), refix=False):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    T = hcp_dataset.get_participants()
    subjects_id = T.participant_id[index]
    for s in T.participant_id[index]:
        data_dir = hcp_dataset.data_dir.format(s)
        if refix:
            C = parcel_hcp_dpconn(data_dir + f'/sub-{s}_tessel-{res}-ReFIX.dpconn.nii')
            nb.save(C, data_dir + f'/sub-{s}_tessel-{res}-ReFIX.pscalar.nii')
        else:
            C = parcel_hcp_dpconn(data_dir + f'/sub-{s}_tessel-{res}.dpconn.nii')
            nb.save(C, data_dir + f'/sub-{s}_tessel-{res}.pscalar.nii')

        print(f"-Saved scalar file for subject {s}, ReFIX={refix}")

if __name__ == "__main__":
    extract_hcp_suit(ses_id='ses-s1', type='CondAll', atlas='MNISymC3')
    extract_hcp_suit(ses_id='ses-s2', type='CondAll', atlas='MNISymC3')
    # extract_hcp_data()
    # avrg_hcp_dpconn()
    # C=parcel_hcp_dpconn(hcp_dir + '/group_tessel-162.dpconn.nii')
    # nb.save(C,hcp_dir + '/group_tessel-162.pscalar.nii')

# Script for importing the MDTB data set from super_cerebellum to general format.
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

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if sys.platform == "win32":  # for windows user
    base_dir = 'Y:/data/FunctionalFusion'

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

def make_mdtb_suit():
    # Make the atlas object
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-2_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('SUIT',mask_img=mask)
    # initialize the data set object
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    for s in T.participant_id:
        deform = mdtb_dataset.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
        mask = mdtb_dataset.suit_dir.format(s) + f'/{s}_desc-cereb_mask.nii'
        atlas_map = am.AtlasMapDeform(suit_atlas.world,deform, mask)
        atlas_map.build(smooth=2.0)
        data,info,str = mdtb_dataset.get_data(s,[atlas_map],'ses-s1')
        #a=mdtb_dataset.get_data_fnames(s,'ses-s1')
        pass

def make_mdtb_fs32k():
    # Make the atlas object
    atlas =[]
    bm_name = ['cortex_left','cortex_right']
    for i,hem in enumerate(['L','R']):
        mask = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-{hem}_mask.label.gii'
        atlas.append(am.AtlasSurface(bm_name[i],mask_gii=mask))
    # initialize the data set object
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    for s in T.participant_id:
        atlas_maps = []
        data = []
        for i,hem in enumerate(['L','R']):
            adir = mdtb_dataset.anatomical_dir.format(s)
            pial = adir + f'/{s}_space-32k_hemi-{hem}_pial.surf.gii'
            white = adir + f'/{s}_space-32k_hemi-{hem}_white.surf.gii'
            mask = adir + f'/{s}_desc-brain_mask.nii'
            atlas_maps.append(am.AtlasMapSurf(atlas[i].vertices,
                            white,pial, mask))
            atlas_maps[i].build()
            # data = mdtb_dataset.get_data(s,[A])
            # data_files=mdtb_dataset.get_data_fnames(s,'ses-s1')
            data.append(np.random.normal(0,1,(100,atlas_maps[i].P))) # am.get_data(data_files,atlas_maps)
        im = am.data_to_cifti(data,atlas_maps)
        nb.save(im,atlas_dir + '/tpl-fs32k/tpl_gs32k_func.dscalar.nii')
        pass

def make_hcp_suit():
    # Make the atlas object for cerebellum suit
    # mask = atlas_dir + '/tpl-MNI152NLin6AsymC' + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('SUIT',mask_img=mask)
    hcp_dataset = DataSetHcpResting(base_dir + '/HCP')

    # create and calculate the atlas map for each participant
    T = hcp_dataset.get_participants()

    # hcp resting state is in MNI space, so we can use the deformation from MNI to suit space
    deform = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
    mask = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152AsymC_res-2_gmcmask2.nii'

    for s in T.participant_id.values:
        # create a mapping based on atlas deformation
        atlas_map = am.AtlasMapDeform(suit_atlas.world,deform, mask)
        # print(f"building atlas map")
        atlas_map.build(smooth=2.0) # smoothing level?
        # print(f"building atlas map is done")

        # get the data based on atlas map
        data = hcp_dataset.get_data_vol(s,[atlas_map],'ses-01')
    pass

def make_hcp_tessel0042():
    """
    Steps to extract the time series for a cortical tesselation in hcp resting state
    1. Make sure that you have a mask for each hemisphere.
        - these masks are created based on non-zero vertices in left and right hemi
        in hcp resting state
        - if the masks are not created,  run test_cifti.get_cortex
    2. load in the label file for the tesselation
    3. apply the mask to the loaded file (only keep the non-zero elements)
    4.extract the time series for each parcel?
    """

    # get the tessellation file
    tessel_dir = atlas_dir + '/tpl-fs32k'
    # get the gifti file for the label
    gii_label = nb.load(tessel_dir + '/Icosahedron-42.32k.L.label.gii')

    # get the gifti of the mask
    gii_mask = nb.load(atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii')

    # get the time series for cortex left
    hcp_dir = base_dir + '/HCP'
    dir = hcp_dir + '/derivatives/100307/func'
    ts_cifti = nb.load(dir+'/ses-01'+'/sub-100307_ses-01_space-fsLR32k_run-01.dtseries.nii')
    bmf = ts_cifti.header.get_axis(1)
    # get the data array with all the time points, all the structures
    ts_array = ts_cifti.get_fdata()
    ts_list = []
    for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
            # just get the cortical surfaces
            if (idx == 0): # just for cortex left (idx = 0)
                # get the values corresponding to the brain model
                bm_vals = ts_array[:, slc]
                # get the voxels/vertices corresponding to the current brain model
                bm_indices = bm.vertex  # excludes medial wall vertices???????
                surf_data = np.zeros((bm_vals.shape[1], bm_vals.shape[0]), dtype=bm_vals.dtype)
                surf_data[:, :] = bm_vals.T

                ts_list.append(surf_data)

    ts_mean = am.get_average_data(data = bm_vals, labels=gii_label, mask = gii_mask)

    # initialize the data set object
    # hcp_dataset = DataSetHcpResting(base_dir + '/HCP')

    return ts_mean

def make_tessels0042_atlas():

    # get the path to the label file
    tessel_file = atlas_dir + '/tpl-fs32k'+'/Icosahedron-42.32k.L.label.gii'

    # get the mask
    mask = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii'
    # create an instance of AtlasSurface for the mask
    myAtlas = am.AtlasSurfaceParcel('cortex_left',tessel_file)
    p_axis = myAtlas.get_parcel_axis()
    return myAtlas,p_axis

def test_pdscalar_file():
    atlas, p_axis = make_tessels0042_atlas()
    D= np.random.normal(0,1,(5,len(p_axis)))
    names = [f'row {r:02}' for r in range(D.shape[0])]
    row_axis = nb.cifti2.ScalarAxis(names)
    header = nb.Cifti2Header.from_axes((row_axis,p_axis))
    cifti_img = nb.Cifti2Image(dataobj=D,header=header)
    nb.save(cifti_img,'test.pscalar.nii')

def hcp_fc_fp(tessel=162):
    """
    Example to get the functional connectivity finger print for hcp resting state
    """
    # Make the atlas object for cerebellum suit
    # mask = atlas_dir + '/tpl-MNI152NLin6AsymC' + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('SUIT',mask_img=mask)
    # hcp resting state is in MNI space, so we can use the deformation from MNI to suit space
    deform = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
    mask = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152AsymC_res-2_gmcmask2.nii'



    # get the tessellation file
    tessel_dir = atlas_dir + '/tpl-fs32k'
    # get the gifti file for the label
    gii_labels = [nb.load(tessel_dir + f'/Icosahedron-{tessel}.32k.L.label.gii'),
                  nb.load(tessel_dir + f'/Icosahedron-{tessel}.32k.R.label.gii')]

    # get the gifti of the mask
    gii_mask = [nb.load(atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii'),
                nb.load(atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii')]

    # create and calculate the atlas map for each participant
    # hcp_dataset = DataSetHcpResting(base_dir + '/HCP')

    # T = hcp_dataset.get_participants()
    # con_fp = [] # connectivity finger print list
    # for s in T.participant_id.values:
    #     # create a mapping based on atlas deformation
    #     atlas_map = am.AtlasMapDeform(hcp_dataset, suit_atlas, s,deform, mask)
    #     # print(f"building atlas map")
    #     atlas_map.build(smooth=2.0) # smoothing level?

    #     # get the connectivity matrix and put it in a list
    #     con_fp.append(hcp_dataset.get_data(s,[atlas_map], gii_labels, gii_mask,'ses-01'))

    return con_fp

def test_atlas_sym():
    mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-2_gmcmask.nii'
    sym_atlas = am.AtlasVolumeSymmetric('MNISymC2',mask_img=mask)

    # First test: Does indx_full work?
    Left = sym_atlas.world[:,sym_atlas.indx_full[0]]
    Right = sym_atlas.world[:,sym_atlas.indx_full[1]]
    Left[0,:] = -Left[0,:]
    assert(np.all(Left==Right))
    # Second test: Does indx_reduced work?
    New = Right[:,sym_atlas.indx_reduced]
    New[0,sym_atlas.indx_full[0]]*=-1
    assert(np.all(New == sym_atlas.world))
    # Third test: flipping
    Flipped = sym_atlas.world[:,sym_atlas.indx_flip]
    Flipped[0,:] = - Flipped[0,:]
    assert(np.all(Flipped == sym_atlas.world))
    pass

def test_atlas_parcelVol():

    data_file = os.path.join(base_dir, 'WMFS', 'derivatives', 'sub-01', 'data', 'sub-01_space-SUIT3_ses-02_CondHalf.dscalar.nii')
    mask_img = os.path.join(atlas_dir, 'tpl-SUIT', 'tpl-SUIT_res-3_gmcmask.nii')
    label_img = os.path.join(atlas_dir, 'tpl-SUIT', 'atl-MDTB10_space-SUIT_dseg.nii')
    MDTB_parcel = am.AtlasVolumeParcel('cerebellum',label_img,mask_img)

    data_img = nb.load(data_file)
    data = data_img.get_fdata()

    data_parcel = MDTB_parcel.agg_data(data)


    return data_parcel

def test_atlas_parcelSurf():

    data_file = os.path.join(base_dir, 'WMFS', 'derivatives', 'sub-01', 'data', 'sub-01_space-fs32k_ses-02_CondHalf.dscalar.nii')
    
    # load left and right data in one file
    data_img = nb.cifti2.load(data_file)
    hemi = ['L', 'R']
    data = data_img.get_fdata()
    # get brain models
    bmf = data_img.header.get_axis(1)
    data_list = []
    for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
        # get the data corresponding to the brain structure
        data_hemi = data[:, slc]

        # get name to be passed on to the AtlasSurfaceParcel object
        name = nam[16:].lower()

        # create atlas parcel object
        # label_img = os.path.join(atlas_dir, 'tpl-fs32k', f'Icosahedron-642_Sym.32k.{hemi[idx]}.label.gii')
        label_img = os.path.join(atlas_dir, 'tpl-fs32k', f'ROI.32k.{hemi[idx]}.label.gii')
        mask_img = os.path.join(atlas_dir, 'tpl-fs32k', f'tpl-fs32k_hemi-{hemi[idx]}_mask.label.gii')
        MDTB_parcel = am.AtlasSurfaceParcel(name,label_img,mask_img)
        data_list.append(MDTB_parcel.agg_data(data_hemi))

    # concatenate into a single array
    data_parcel = np.concatenate(data_list, axis = 1)
    print(data_parcel.shape)

    return data_parcel

if __name__ == "__main__":
    # make_mdtb_suit()
    # test_atlas_sym()
    test_atlas_parcelSurf()
    pass


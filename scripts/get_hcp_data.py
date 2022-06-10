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


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/HCP'
atlas_dir = base_dir + '/Atlases'
hem_name = ['cortex_left','cortex_right']

def get_hcp_data(res=162):    
    # Make the atlas object
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    
    # initialize the data set object
    hcp_dataset = DataSetHcpResting(data_dir)
    
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
    for s in T.participant_id[0:2]:
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
        nb.save(cifti_img, dest_dir + f'/sub-{s}.dpconn.nii')


if __name__ == "__main__":
    get_hcp_data()
    pass


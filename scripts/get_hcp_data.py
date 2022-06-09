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
    atlas_map.build(smooth=None)

    # Get the parcelation 
    surf_parcel =[] 
    for i,h in enumerate(['L','R']): 
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i],gifti))

    T = hcp_dataset.get_participants() 
    for s in T.participant_id[0:2]:
        print(f'Extract {s}')
        data,info,names = hcp_dataset.get_cereb_connectivity(s,atlas_map,   surf_parcel)
        pass


        C=am.data_to_cifti(data,[atlas_map],names)
        dest_dir = mdtb_dataset.data_dir.format(s)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        nb.save(C, dest_dir + f'/{s}_space-SUIT3_{ses_id}_CondSes.dscalar.nii')
        info.to_csv(dest_dir + f'/{s}_{ses_id}_info-CondSes.tsv',sep='\t')



if __name__ == "__main__":
    get_hcp_data()
    pass


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

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
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
        atlas_map = am.AtlasMapDeform(mdtb_dataset, suit_atlas, s,deform, mask)
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
            atlas_maps.append(am.AtlasMapSurf(mdtb_dataset, atlas[i],
                            s,white,pial, mask))
            atlas_maps[i].build()
            # data = mdtb_dataset.get_data(s,[A])
            # data_files=mdtb_dataset.get_data_fnames(s,'ses-s1')
            data.append(np.random.normal(0,1,(100,atlas_maps[i].P))) # am.get_data(data_files,atlas_maps)
        im = am.data_to_cifti(data,atlas_maps)
        nb.save(im,atlas_dir + '/tpl-fs32k/tpl_gs32k_func.dscalar.nii')
        pass

def make_hcp_suit():
    # Make the atlas object for cerebellum suit
    mask = atlas_dir + '/tpl-MNI152NLin6AsymC' + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('SUIT',mask_img=mask)
    hcp_dataset = DataSetHcpResting(base_dir + '/HCP')
    
    # create and calculate the atlas map for each participant
    T = hcp_dataset.get_participants()

    # hcp resting state is in MNI space, so we can use the deformation from MNI to suit space
    deform = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'

    for s in T.participant_id.values: ## for testing, just one participant is used
        # create a mapping based on atlas deformation 
        atlas_map = am.AtlasMapDeform(hcp_dataset, suit_atlas, s,deform, mask)
        atlas_map.build(smooth=2.0) # smoothing level?

        # get the data based on atlas map
        data = hcp_dataset.get_data(s,[atlas_map],'ses-01')

        return data

if __name__ == "__main__":
    make_mdtb_suit()


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
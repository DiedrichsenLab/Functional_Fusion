# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np 
import atlas_map as am
from dataset import DataSetMDTB

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

def make_mdtb_suit(): 
    # Make the atlas object 
    mask = atlas_dir + '/tpl-SUIT3/tpl-SUIT_res-2_mask.nii'
    suit3_atlas = am.AtlasVolumetric('SUIT3',mask_img=mask)
    # initialize the data set object 
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    atlas_maps = []
    for s in T.participant_id:
        deform = mdtb_dataset.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
        mask = mdtb_dataset.suit_dir.format(s) + f'/{s}_desc-cereb_mask.nii'
        A = am.AtlasMapDeform(mdtb_dataset, suit3_atlas, s,deform, mask)
        A.build(smooth=2.0)
        # data = mdtb_dataset.get_data(s,[A])
        #a=mdtb_dataset.get_data_fnames(s,'ses-s1')
        # data = am.get_data(a,[A])
        am.save_data_to_cifti(data,atlas_maps)
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
    atlas_maps = []
    for s in T.participant_id:
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
        # data = am.get_data(data_files,atlas_maps)
        am.save_data_to_cifti(np.zeros((4,4)),atlas_maps)
        pass


if __name__ == "__main__":
    make_mdtb_fs32k()


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np 
from atlas_map import AtlasVolumetric, AtlasMapDeform
from dataset import DataSetMDTB

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

def make_mdtb_suit3(): 
    # Make the atlas object 
    mask = atlas_dir + '/atl-SUIT3/atl-SUIT3_mask.nii'
    suit3_atlas = AtlasVolumetric('SUIT3',mask_img=mask)
    # initialize the data set object 
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    atlas_maps = []
    for s in T.participant_id:
        deform = mdtb_dataset.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
        mask = mdtb_dataset.data_dir.format(s) + f'/ses-s1/{s}_ses-s1_mask.nii'
        A = AtlasMapDeform(mdtb_dataset, suit3_atlas, s,deform, mask)
        A.build()
        A.save()

if __name__ == "__main__":
    make_mdtb_suit3()


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
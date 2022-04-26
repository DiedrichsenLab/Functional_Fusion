# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np 
from atlas_map import AtlasVolumetric, AtlasMapDeform



base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

if __name__ == "__main__":
    mask = atlas_dir + '/atl-SUIT3/atl-SUIT3_mask.nii'
    suit3_atlas = AtlasVolumetric('SUIT3',mask_img=mask)
    pass
    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
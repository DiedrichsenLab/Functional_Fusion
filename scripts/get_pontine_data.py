# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
# Discuss: This is causing trouble in compatibility
# sys.path.append(
#     '/Users/callithrix/Documents/Projects/Functional_Fusion/code/shared/Functional_Fusion/') # can be removed before push, but currently that is the best way to import atlas_map for me
import atlas_map as am
from dataset import DataSetPontine
import nibabel as nb
import SUITPy as suit


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Pontine7T'
atlas_dir = base_dir + '/Atlases'

def extract_pontine_suit(ses_id='ses-01',type='taskHalf',atlas='SUIT3'):
    p7_dataset = DataSetPontine(data_dir)
    p7_dataset.extract_all_suit(ses_id,type,atlas)

def extract_pontine_fs32k(ses_id='ses-01',type='taskHalf'):
    p7_dataset = DataSetPontine(data_dir)
    p7_dataset.extract_all_fs32k(ses_id,type)

def show_pontine_suit(subj,sess,cond):
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    pontine_dataset = DataSetPontine(data_dir)
    T = pontine_dataset.get_participants()
    s = T.participant_id[subj]
    ses = f'ses-{sess:02d}'
    C = nb.load(pontine_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_taskHalf.dscalar.nii')
    D = pd.read_csv(pontine_dataset.data_dir.format(s) + f'/{s}_{ses}_info-taskHalf.tsv',sep='\t')
    X = C.get_fdata()
    Nifti = suit_atlas.data_to_nifti(X)
    surf_data = suit.flatmap.vol_to_surf(Nifti)
    fig = suit.flatmap.plot(surf_data[:,cond],render='plotly')
    fig.show()
    print(f'Showing {D.cond_name[cond]}')
    pass

if __name__ == "__main__":
    extract_pontine_suit(ses_id='ses-01',type='taskHalf')
    extract_pontine_fs32k(ses_id='ses-01',type='taskHalf')

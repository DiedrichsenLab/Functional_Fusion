# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetPontine
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Pontine'
atlas_dir = base_dir + '/Atlases'


def show_pontine_group(type='TaskHalf', atlas='SUIT3', cond=0, info_column='task_name', savefig=False):
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum', mask_img=mask)
    p7_dataset = DataSetPontine(data_dir)
    C = nb.load(p7_dataset.data_dir.split('/{0}')[0] +
                f'/group/group_space-{atlas}_{type}.dscalar.nii')
    D = pd.read_csv(p7_dataset.data_dir.split('/{0}')[0] +
                    f'/group/group_info-{type}.tsv', sep='\t')
    X = C.get_fdata()


def extract_pontine_group(type='TaskHalf', atlas='SUIT3', info_column='task_name'):
    p7_dataset = DataSetPontine(data_dir)
    p7_dataset.group_average_data(type, atlas, info_column)

def extract_pontine_suit(ses_id='ses-01',type='TaskHalf',atlas='SUIT3'):
    p7_dataset = DataSetPontine(data_dir)
    p7_dataset.extract_all_suit(ses_id,type,atlas)

def extract_pontine_fs32k(ses_id='ses-01',type='TaskHalf'):
    p7_dataset = DataSetPontine(data_dir)
    p7_dataset.extract_all_fs32k(ses_id,type)

def show_pontine_suit(subj,sess,cond):
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    p7_dataset = DataSetPontine(data_dir)
    T = p7_dataset.get_participants()
    s = T.participant_id[subj]
    ses = f'ses-{sess:02d}'
    C = nb.load(p7_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_TaskHalf.dscalar.nii')
    D = pd.read_csv(p7_dataset.data_dir.format(s) + f'/{s}_{ses}_info-TaskHalf.tsv',sep='\t')
    X = C.get_fdata()
    Nifti = suit_atlas.data_to_nifti(X)
    surf_data = suit.flatmap.vol_to_surf(Nifti)
    fig = suit.flatmap.plot(surf_data[:,cond],render='plotly')
    fig.show()
    print(f'Showing {D.cond_name[cond]}')
    pass

if __name__ == "__main__":
    # extract_pontine_group(type='TaskHalf', atlas='MNISymC3')
    #  extract_pontine_fs32k(ses_id='ses-01',type='TaskHalf')
    # extract_pontine_suit(ses_id='ses-01', type='TaskHalf', atlas='MNISymC2')
    # show_pontine_group(type='TaskHalf', atlas='SUIT3',
    #                    cond='all', savefig=True)

    dataset = DataSetPontine(data_dir)
    dataset.group_average_data(atlas='MNISymC3')
    dataset.plot_cerebellum(savefig=True, atlas='MNISymC3', colorbar=True)
    pass

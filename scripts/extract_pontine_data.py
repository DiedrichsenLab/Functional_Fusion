# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
import atlas_map as am
from dataset import DataSetPontine
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

    if cond == 'all':
        conditions = D[info_column]
        # -- as subplot --
        # dim = int(np.ceil(np.sqrt(len(conditions))))
        # fig, axs = plt.subplots(dim,dim)
        # for i, c in enumerate(conditions):
        #     Nifti = suit_atlas.data_to_nifti(X[i, :])
        #     surf_data = suit.flatmap.vol_to_surf(Nifti)
        #     axs[int(i % dim), int(i / dim)]=suit.flatmap.plot(
        #         surf_data, render='matplotlib', new_figure=False)
        #     axs[int(i % dim), int(i / dim)].set_title(c)
        # fig.show()
        # -- each in seperate figures --
        dest_dir = p7_dataset.data_dir.split('/{0}')[0] + f'/group/figures/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        for i, c in enumerate(conditions):
            Nifti = suit_atlas.data_to_nifti(X[i, :])
            surf_data = suit.flatmap.vol_to_surf(Nifti)
            fig = suit.flatmap.plot(
                surf_data, render='matplotlib', new_figure=True)
            fig.set_title(c)
            # save figure
            if savefig:
                plt.savefig(dest_dir + f'group_{c}.png')
            plt.clf()
            pass

    else:
        Nifti = suit_atlas.data_to_nifti(X[cond, :])
        surf_data = suit.flatmap.vol_to_surf(Nifti)
        fig = suit.flatmap.plot(surf_data, render='plotly')
        fig.show()
        print(f'Showing {D[info_column][cond]}')
        pass

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
    extract_pontine_suit(ses_id='ses-01', type='TaskHalf', atlas='MNISymC2')
    # show_pontine_group(type='TaskHalf', atlas='SUIT3',
    #                    cond='all', savefig=True)
    pass

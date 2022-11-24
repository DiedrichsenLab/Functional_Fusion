# Script for importing the IBC data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import *
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt



base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

data_dir = base_dir + '/IBC'
atlas_dir = base_dir + '/Atlases'


def show_ibc_group(ses_id='ses-hcp1', type='CondHalf', atlas='MNISymC3', cond=0, info_column='names', savefig=False):
    if (atlas == 'MNISymC3'):
        mask = atlas_dir + '/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum', mask_img=mask)
    
    ibc_dataset = DataSetIBC(data_dir)
    C = nb.load(ibc_dataset.data_dir.split('/{0}')[0] +
                f'/group/group_{ses_id}_space-{atlas}_{type}.dscalar.nii')
    D = pd.read_csv(ibc_dataset.data_dir.split('/{0}')[0] +
                    f'/group/group_{ses_id}_info-{type}.tsv', sep='\t')
    X = C.get_fdata()

    if cond == 'all':
        conditions = D[info_column]
        # -- each in seperate figures --
        dest_dir = ibc_dataset.data_dir.split('/{0}')[0] + f'/group/figures/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        for i, c in enumerate(conditions):
            Nifti = suit_atlas.data_to_nifti(X[i, :])
            surf_data = suit.flatmap.vol_to_surf(Nifti, atlas[:-1])
            fig = suit.flatmap.plot(
                surf_data, render='matplotlib', new_figure=True)
            fig.set_title(c)
            # save figure
            if savefig:
                plt.savefig(dest_dir + f'group_{ses_id}_{c}.png')
            plt.clf()
            pass

    else:
        Nifti = suit_atlas.data_to_nifti(X[cond, :])
        surf_data = suit.flatmap.vol_to_surf(Nifti)
        fig = suit.flatmap.plot(surf_data, render='plotly')
        fig.show()
        print(f'Showing {D.cond_name[cond]}')
        pass

def extract_all(atlas='MNISym3'):
    ibc_dataset = DataSetIBC(data_dir)
    info = ibc_dataset.get_participants()
    for ses in ibc_dataset.sessions:
        print(f'extracting {ses}')
        if atlas == 'fs32k':
            ibc_dataset.extract_all_fs32k(ses,type='CondHalf')
        else:
            ibc_dataset.extract_all_suit(ses,type='CondHalf',atlas=atlas)


def group_average(atlas='MNISymC3'):
    type = 'CondHalf'
    ibc_dataset = DataSetIBC(data_dir)
    # ---- Extract all data 
    # info = ibc_dataset.get_participants(
    #     # 
    # --- Get group average ---
    for ses in ibc_dataset.sessions:
        ibc_dataset.group_average_data(
            ses_id=ses, type=type, atlas=atlas)
        # write session tsv file for group average
        s = ibc_dataset.get_participants().participant_id[0]
        D = pd.read_csv(Path(ibc_dataset.data_dir.format(s)) /
                        f'{s}_{ses}_info-{type}.tsv', sep='\t')
        D = D.drop(columns=['sn', 'sess', 'run']).drop_duplicates(keep='first')
        D.to_csv(ibc_dataset.data_dir.split('/{0}')[0] +
                        f'/group/group_{ses}_info-{type}.tsv', sep='\t')


def show_group_average(atlas='MNISymC3'):
    ibc_dataset = DataSetIBC(data_dir)
    for ses in ibc_dataset.sessions:
        show_ibc_group(ses_id=ses, type='CondHalf',
                atlas=atlas, cond='all', savefig=True)

    pass

if __name__ == "__main__":
    # extract_all('fs32k')
    group_average(atlas='fs32k')

    # parcel_mdtb_fs32k()
    # 
 
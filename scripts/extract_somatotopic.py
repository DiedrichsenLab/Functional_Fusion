# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
import atlas_map as am
from Functional_Fusion.dataset import DataSetSomatotopic
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Somatotopic'
atlas_dir = base_dir + '/Atlases'

def extract_somatotopic(ses_id='ses-01',type='CondHalf',atlas='MNISymC3'):
    dataset = DataSetSomatotopic(data_dir)
    dataset.extract_all(ses_id,type,atlas)


def show_group(type='CondHalf', atlas='SUIT3', session='motor', cond='all', info_column='cond_name', savefig=False):
    
    dataset = DataSetSomatotopic(data_dir)

    myatlas = am.get_atlas(atlas, dataset.atlas_dir)

    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum', mask_img=mask)

    # Load group average
    C = nb.load(dataset.data_dir.split('/{0}')[0] +
                f'/group/data/group_ses-{session}_space-{atlas}_{type}.dscalar.nii')
    D = pd.read_csv(dataset.data_dir.split('/{0}')[0] +
                    f'/group/data/group_ses-{session}_info-{type}.tsv', sep='\t')
    X = C.get_fdata()

    if cond == 'all':
        conditions = D[info_column]
        # -- each in seperate figures --
        dest_dir = dataset.data_dir.split('/{0}')[0] + f'/group/figures/'
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


if __name__ == "__main__":
    # --- Extracting Estimates ---
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')


    # --- Group Average ---
    dataset = DataSetSomatotopic(data_dir)
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')

    
    # --- Show group average ---
    dataset.plot_group_cerebellum(savefig=True, colorbar=True)
    pass

# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
# import mat73
import numpy as np
import sys
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetLanguage, DataSetPontine
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import Functional_Fusion.connectivity as conn


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Language'
atlas_dir = base_dir + '/Atlases'


# def extract_langloc_group(atlas='SUIT3'):
#     LL_dataset = DataSetLangloc(data_dir)
#     LL_dataset.group_average_data(type, atlas)

def extract_language(atlas='MNISymC2'):
    LL_dataset = DataSetLanguage(data_dir)
    LL_dataset.extract_all(ses_id = 'ses-02',type = 'CondHalf', atlas = atlas)

# def extract_langloc_fs32k(ses_id='ses-01',type='TaskHalf'):
#     LL_dataset = DataSetPontine(data_dir)
#     LL_dataset.extract_all_fs32k(ses_id,type)

# def show_langloc_suit(subj,sess,cond):
#     mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
#     suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
#     LL_dataset = DataSetPontine(data_dir)
#     T = LL_dataset.get_participants()
#     s = T.participant_id[subj]
#     ses = f'ses-{sess:02d}'
#     C = nb.load(LL_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_TaskHalf.dscalar.nii')
#     D = pd.read_csv(LL_dataset.data_dir.format(s) + f'/{s}_{ses}_info-TaskHalf.tsv',sep='\t')
#     X = C.get_fdata()
#     Nifti = suit_atlas.data_to_nifti(X)
#     surf_data = suit.flatmap.vol_to_surf(Nifti)
#     fig = suit.flatmap.plot(surf_data[:,cond],render='plotly')
#     fig.show()
#     print(f'Showing {D.cond_name[cond]}')
#     pass

if __name__ == "__main__":
    # extract_langloc_group( atlas='MNISymC3')
    # extract_langloc_fs32k(ses_id='ses-01',type='TaskHalf')
    # extract_langloc_suit(ses_id='ses-01', type='TaskHalf', atlas='MNISymC2')
    # show_langloc_group(type='TaskHalf', atlas='SUIT3',
    #                    cond='all', savefig=True)

    # dataset = DataSetLanguageLocalizer(data_dir)
    # dataset.group_average_data(atlas='MNISymC3')
    # dataset.plot_cerebellum(savefig=True, atlas='MNISymC3', colorbar=True)

    # extract_language()

    # Exctract Rest timeseries & connectivity fingerprint
    dname = 'Language'
    session = 'rest'
    # session = 'localizer_cond'
    T = pd.read_csv(
            data_dir + '/participants.tsv', delimiter='\t')
    subject_subset_indices = T.participant_id[T['ses-rest'] == 1].index.tolist()
    subject_subset = T.participant_id[T['ses-rest'] == 1].tolist()
    lang_dataset = DataSetLanguage(data_dir)
    # Extract non-fix Tseries
    # lang_dataset.extract_all(ses_id=f'ses-{session}', type='Tseries', atlas='MNISymC3', subj=subject_subset_indices)
    # lang_dataset.extract_all(ses_id=f'ses-{session}', type='Tseries', atlas='fs32k', subj=subject_subset_indices)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06Run', space='MNISymC3', ses_id=f'ses-{session}', subj=subject_subset)
    conn.get_connectivity_fingerprint(dname,
                                      type='Fus06Run', space='fs32k', ses_id=f'ses-{session}', subj=subject_subset)
    
    # Exctract fix-cleaned Tseries
    # lang_dataset.extract_all(ses_id=f'ses-{session}', type='FixTseries', atlas='MNISymC3', subj=subject_subset_indices)
    # lang_dataset.extract_all(ses_id=f'ses-{session}', type='FixTseries', atlas='fs32k', subj=subject_subset_indices)
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Fus06FixRun', space='MNISymC3', ses_id=f'ses-{session}', subj=subject_subset)
    pass

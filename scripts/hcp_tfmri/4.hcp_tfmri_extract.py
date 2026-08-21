"""
Script for importing the Language dataset to general format.

Created Dec 2023
Author: Bassel Arafat
"""

import pandas as pd
from pathlib import Path
# import mat73
import time
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetHcpTask
import Functional_Fusion.util as ut

base_dir = '/home/dzhi/eris_mount'
if not os.path.exists(base_dir):
    base_dir = '/data/tge'

data_dir = base_dir + '/Tian/HCP_img'
atlas_dir = base_dir + '/Tian/UKBB_full/imaging/Atlases'

types = ['ZstatHalf']
atlases  = ['fs32k']
session_list = ['ses-task']

def mask_hcptask_fs32k(ses_id='ses-s1', type='CondHalf', high_percent=0.1, low_percent=0.1,
                    smooth=None, z_transfer=False, binarized=False):
    hcptask_dataset = DataSetHcpTask(data_dir)
    T = hcptask_dataset.get_participants(f'/subj_list/HCP203_test_set_filtered.tsv')

    for s in T.participant_id:
        print(f'Mask data for {s} fs32k {ses_id} in high {high_percent} low {low_percent} ...')

        start = time.perf_counter()
        if smooth is not None:
            file = hcptask_dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        else:
            file = hcptask_dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'

        ut.mask_fs32k_data(file, high_percent=high_percent, low_percent=low_percent,
                           z_transfer=z_transfer, binarized=binarized)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


if __name__ == "__main__":
    types = ['CondAll', 'CondHalf']
    atlases  = ['fs32k']
    session_list = ['ses-task']


    dataset = DataSetHcpTask(data_dir)
    for ses in session_list:
        print(f'extracting session {ses}')
        participants_tsv = pd.read_csv(f'{data_dir}/subj_list/HCP200_test.tsv',sep = '\t')
        subj_list = participants_tsv['participant_id'].tolist()
        hcp_subj_ind = dataset.get_participants().index[dataset.get_participants()['participant_id'].isin(subj_list)].tolist()

        for type in types:
            print(f'extracting type {type}')
            for atlas in atlases:
                print(f'extracting atlas: {atlas}')
                dataset.extract_all(ses_id=ses,type=type, atlas=atlas, smooth=None, subj=hcp_subj_ind)

        # for participant_id in subj_list:
        #     # Make info
        #     dest_dir = dataset.base_dir + f'/derivatives/{participant_id}/data/'
        #     info = pd.read_csv(dest_dir + f'{participant_id}_{ses}_CondAll.tsv', sep='\t')
        #
        #     if info.get('half') is None:
        #         info.insert(0, 'half', 1)
        #         info.to_csv(dest_dir + f'{participant_id}_{ses}_CondAll.tsv', sep='\t')
        #         print(f'Added half column for subject {participant_id}')
        #     else:
        #         print('already has half')
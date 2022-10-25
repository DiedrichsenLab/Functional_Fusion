# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import *
import nibabel as nb
import SUITPy as suit


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

data_dir = base_dir + '/IBC'
atlas_dir = base_dir + '/Atlases'


if __name__ == "__main__":
    # parcel_mdtb_fs32k()
    ibc_dataset = DataSetIBC(data_dir)
    info = ibc_dataset.get_participants()
    for ses in ibc_dataset.sessions:
         ibc_dataset.extract_all_suit(ses,type='CondHalf',atlas='MNISymC3')
    # 
    type = 'CondHalf'
    ibc_dataset = DataSetIBC(data_dir)
    # ---- Extract all data 
    info = ibc_dataset.get_participants()
    for ses in ibc_dataset.sessions:
        ibc_dataset.extract_all_suit(ses, type='CondHalf', atlas='SUIT3')
    # 
    # --- Get group average ---
    for ses in ibc_dataset.sessions:
        ibc_dataset.group_average_data(
            ses_id=ses, type=type, atlas='MNISymC3')
        # write session tsv file for group average
        s = ibc_dataset.get_participants().participant_id[0]
        D = pd.read_csv(Path(ibc_dataset.data_dir.format(s)) /
                        f'{s}_{ses}_info-{type}.tsv', sep='\t')
        D = D.drop(columns=['sn', 'sess', 'run']).drop_duplicates(keep='first')
        D.to_csv(ibc_dataset.data_dir.split('/{0}')[0] +
                        f'/group/group_{ses}_info-{type}.tsv', sep='\t')


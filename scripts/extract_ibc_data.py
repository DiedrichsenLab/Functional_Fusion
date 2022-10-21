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

    pass


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
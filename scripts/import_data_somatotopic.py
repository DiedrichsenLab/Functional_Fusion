"""
Script for importing the Pontine dataset to general format.

Created Sep 2022
Author: caro nettekoven
"""

import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
from import_data import *
from Functional_Fusion.dataset import DataSetSomatotopic

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'


src_base_dir = base_dir + '/Cerebellum//Somatotopic'
dest_base_dir = base_dir + '/FunctionalFusion/Somatotopic'


def import_data():
    dataset = DataSetSomatotopic(dest_base_dir)
    T = dataset.get_participants()
    for p, participant_id in enumerate(T.participant_id):
        fnames, reginfo = dataset.get_data_fnames(
            participant_id, session_id='ses-01')
        pass

if __name__ == '__main__':
    
    # --- Importing Estimates ---
    import_data()

    # --- Importing Estimates ---
    dataset = DataSetSomatotopic(dest_base_dir)
    T = dataset.get_participants()
    dataset.extract_all(type='CondAll')



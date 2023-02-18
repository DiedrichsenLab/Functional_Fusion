# Script for getting all the HCP data for cerebellar-cortical connectivity
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetHcpResting
import nibabel as nb
import SUITPy as suit
import os
import sys
import matplotlib.pyplot as plt
from ProbabilisticParcellation.util import plot_multi_flat, plot_data_flat


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

hcp_dir = base_dir + '/HCP'
atlas_dir = base_dir + '/Atlases'
hem_name = ['cortex_left', 'cortex_right']


def extract_hcp(ses_id='ses-rest1', type='Run', atlas='MNISymC3'):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    T = hcp_dataset.get_participants()
    hcp_dataset.get_data_fnames(T.iloc[0].participant_id, ses_id)
    hcp_dataset.extract_all(ses_id, type, atlas)


if __name__ == "__main__":
    extract_hcp()

    pass

# Script for getting all the HCP data for cerebellar-cortical connectivity
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
import atlas_map as am
from dataset import DataSetLanguage
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Language'
atlas_dir = base_dir + '/Atlases'


def extract_language_suit(ses_id='ses-01', type='TaskHalf', atlas='SUIT3'):
    lang_dataset = DataSetLanguage(data_dir)
    lang_dataset.extract_all_suit(ses_id, type, atlas)

if __name__ == "__main__":
    # extract_language_group(type='TaskAll', atlas='MNISymC3')
    extract_language_suit(ses_id='ses-01', type='TaskAll')

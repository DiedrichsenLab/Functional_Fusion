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
    dataset.plot_cerebellum(subject='group', savefig=True, colorbar=True)
    pass

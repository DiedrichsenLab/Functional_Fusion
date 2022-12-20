# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
import atlas_map as am
from dataset import DataSetDemand
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/Demand'
atlas_dir = base_dir + '/Atlases'

def extract_demand(ses_id='ses-01',type='CondHalf',atlas='MNISymC3'):
    de_dataset = DataSetDemand(data_dir)
    de_dataset.extract_all(ses_id,type,atlas)

if __name__ == "__main__":
    # extract_demand(ses_id='ses-01', type='CondHalf', atlas='MNISymC2')
    de_dataset = DataSetDemand(data_dir)
    # de_dataset.group_average_data(ses_id='ses-01',type='CondHalf',atlas='MNISymC2')
    # de_dataset.group_average_data(ses_id='ses-01',type='CondHalf',atlas='MNISymC3')
    # de_dataset.group_average_data(ses_id='ses-01',type='CondHalf',atlas='fs32k')

    de_dataset.plot_group_cerebellum(savefig=True, colorbar=True)
    pass

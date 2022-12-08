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



def extract_demand_suit(ses_id='ses-01',type='CondHalf',atlas='MNISymC3'):
    de_dataset = DataSetDemand(data_dir)
    de_dataset.extract_all_suit(ses_id,type,atlas)

def extract_pontine_fs32k(ses_id='ses-01',type='TaskHalf'):
    de_dataset = DataSetDemand(data_dir)
    de_dataset.extract_all_fs32k(ses_id,type)

if __name__ == "__main__":
    # extract_pontine_group(type='TaskHalf', atlas='MNISymC3')
    #  extract_pontine_fs32k(ses_id='ses-01',type='TaskHalf')
    extract_demand_suit(ses_id='ses-01', type='CondHalf', atlas='MNISymC3')
    # show_pontine_group(type='TaskHalf', atlas='SUIT3',
    #                    cond='all', savefig=True)
    pass

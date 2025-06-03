"""
Script for creating group average data for HCP tfMRI dataset.

Created May 2025
Author: Bassel Arafat
"""

from pathlib import Path
from Functional_Fusion.dataset import DataSetHcpTask


base_dir = 'Y:/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'


data_dir = base_dir + '/HCP_tfMRI'
atlas_dir = base_dir + '/Atlases'


type='CondAll'
atlas  ='fs32k'


dataset = DataSetHcpTask(data_dir)
dataset.group_average_data(ses_id='ses-task', type=type, atlas=atlas)
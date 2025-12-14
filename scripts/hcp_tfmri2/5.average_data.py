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


data_dir = base_dir + '/HCPur100'
atlas_dir = base_dir + '/Atlases'


type='CondAll'
atlas  ='fs32k'


dataset = DataSetHcpTask(data_dir)
subj_list = ['sub-101309','sub-103111','sub-103414','sub-103818','sub-105014','sub-110411','sub-117122','sub-118730','sub-123117','sub-124422','sub-127630','sub-128632','sub-130013']
dataset.group_average_data(ses_id='ses-task2', type=type, atlas=atlas, subj=subj_list)
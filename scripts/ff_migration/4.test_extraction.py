"""
Script for testing the extraction of data from the MDTB dataset with a universal condense_data in the dataset class.
"""

from pathlib import Path
from Functional_Fusion.dataset import DataSetMDTB


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion_new'


data_dir = base_dir + '/MDTB'

dataset_obj = DataSetMDTB(data_dir)
dataset_obj.extract_all(ses_id = 'ses-s1',type = 'CondHalf', atlas = 'MNISymC3', smooth=None, subj=[0])
 

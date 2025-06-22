# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import sys
import atlas_map as am
from Functional_Fusion.dataset import DataSetDmcc
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Users/lshahsha/Documents/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/DMCC'
atlas_dir = base_dir + '/Atlases'

def extract_dmcc(ses_id='ses-stroop_bas',type='CondHalf',atlas='MNISymC3', smooth = 0):
    dataset = DataSetDmcc(data_dir)
    dataset.extract_all(ses_id=ses_id,
                        type=type,
                        atlas=atlas,
                        smooth=smooth)




if __name__ == "__main__":
    # --- Extracting Estimates ---
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # extract_somatotopic(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    
    extract_dmcc(ses_id='ses-axcpt_bas', type='CondAll', atlas='SUIT3', smooth=2)
    extract_dmcc(ses_id='ses-axcpt_bas', type='CondAll', atlas='fs32k', smooth=2)
    extract_dmcc(ses_id='ses-axcpt_bas', type='CondAll', atlas='MNI2009cAsymBg2', smooth=2)

    extract_dmcc(ses_id='ses-cuedts_bas', type='CondAll', atlas='SUIT3', smooth=2)
    extract_dmcc(ses_id='ses-cuedts_bas', type='CondAll', atlas='fs32k', smooth=2)
    extract_dmcc(ses_id='ses-cuedts_bas', type='CondAll', atlas='MNI2009cAsymBg2', smooth=2)

    extract_dmcc(ses_id='ses-stern_bas', type='CondAll', atlas='SUIT3', smooth=2)
    extract_dmcc(ses_id='ses-stern_bas', type='CondAll', atlas='fs32k', smooth=2)
    extract_dmcc(ses_id='ses-stern_bas', type='CondAll', atlas='MNI2009cAsymBg2', smooth=2)

    extract_dmcc(ses_id='ses-stroop_bas', type='CondAll', atlas='SUIT3', smooth=2)
    extract_dmcc(ses_id='ses-stroop_bas', type='CondAll', atlas='fs32k', smooth=2)
    extract_dmcc(ses_id='ses-stroop_bas', type='CondAll', atlas='MNI2009cAsymBg2', smooth=2)


    # --- Group Average ---
    # dataset = DataSetSomatotopic(data_dir)
    # # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')

    
    # # --- Show group average ---
    # dataset.plot_cerebellum(subject='group', savefig=True, colorbar=True)
    pass

import nibabel as nib
from pathlib import Path
import os.path as op
import preprocessing.paths as paths
import Functional_Fusion.util as util
import numpy as np
import pandas as pd


# set paths
base_dir = paths.set_base_dir()
rest_dir = f'{base_dir}/../Cerebellum/super_cerebellum/resting_state/'
data_dir = f'{rest_dir}/imaging_data'




def calculate_variances(subject, run, clean=False, roi='cerebellum', save=False):
    """Calculate the variance removed by FIX for a given subject and run.
        
        Args:
            subject (str): subject ID
            run (int): run number
            clean (bool): whether to use clean data
            roi (str): ROI to use as seed
            save (bool): whether to save the correlation map
            
        Returns:
            correlation_map (np.ndarray): correlation map
    """
    # Load data
    data_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/filtered_func_data.nii.gz'
    data = nib.load(data_path).get_fdata()

    data_clean = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/filtered_func_data_clean.nii.gz'
    data_clean = nib.load(data_clean).get_fdata()

    # Load mask that excludes non-brain and non-skull background voxels
    mask_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/mask.nii.gz'
    mask = nib.load(mask_path).get_fdata()
    data = data[mask > 0]
    data_clean = data_clean[mask > 0]

    # Calculate variance
    variance = np.var(data)
    variance_clean = np.var(data_clean)

    return variance, variance_clean

def calculate_all_variances(save=True):
    """Calculate the variance removed by FIX for all subjects and runs.
        

    """
    
    variances = []
    folders = Path(data_dir).glob('s*')
    for s, subject in enumerate(folders):
        subject = str(subject).split('/')[-1]
        for run in [1, 2]:
            variance, variance_cleaned = calculate_variances(subject, run)
            #  Calculate percentage of variance removed
            variance_removed = (variance - variance_cleaned) / variance * 100
            variances.append((variance, variance_cleaned, variance_removed, subject, run))
            #
            
    Var = pd.DataFrame(variances, columns=['prefix', 'postfix', 'perc_removed', 'subject', 'run'])
    if save:
        Var.to_csv(f'{rest_dir}/fix_ica/variance.csv')
    
    return Var


if __name__ == "__main__":
    calculate_all_variances()
    
    
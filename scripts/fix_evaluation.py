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



def calculate_variance(subject, run, clean=False, roi='cerebellum', save=False):
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
    if clean:
        data_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/filtered_func_data_clean.nii.gz'
    data = nib.load(data_path).get_fdata()

    # Load mask that excludes non-brain and non-skull background voxels
    mask_path = f'{rest_dir}/imaging_data/{subject}/run{run:02d}.feat/mask.nii.gz'
    mask = nib.load(mask_path).get_fdata()
    data = data[mask > 0]

    # Calculate variance
    variance = np.var(data)

    return variance

def calculate_all_removed_variance():
    """Calculate the variance removed by FIX for all subjects and runs.
        

    """
    
    variances = []
    for s, subject in enumerate(rest_dir.glob('s*')):
        for run in [1, 2]:
            variance = calculate_variance(subject, run)
            variances.append((variance, subject, run))
            

    Var = pd.DataFrame(variances, columns=['variance', 'subject', 'run'])


if __name__ == "__main__":
    calculate_all_removed_variance()
    
    
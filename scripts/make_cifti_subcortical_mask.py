# Script for importing the Pontine data set to general format.
import pandas as pd
import shutil
from pathlib import Path
import nitools as nt 
import sys, os, time
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetDemand
import Functional_Fusion.util as fut
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import numpy as np

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion_new'

data_dir = base_dir + '/Demand'
atlas_dir = base_dir + '/Atlases'

def make_mask(ses_id='ses-01',type='CondHalf',atlas='MNISymC3'):
    de_dataset = DataSetDemand(data_dir)
    fname,_ = de_dataset.get_data_fnames(participant_id='sub-01',session_id='ses-01')
    cifti = nb.load(fname[1])
    V = nt.volume_from_cifti(cifti)
    X = V.get_fdata().squeeze()
    X=X>0
    X = X.astype(np.uint8)
    V= nb.Nifti1Image(X, V.affine)
    
    V.to_filename(f'{atlas_dir}/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_desc_subcortexmask.nii.gz')

    pass 

if __name__ == "__main__":
    make_mask()
    pass

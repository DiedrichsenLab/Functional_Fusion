# Working script for all extractions and types

import pandas as pd
import shutil
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as util
from Functional_Fusion.matrix import indicator
import nibabel as nb

base_dir = util.get_base_dir()

def extract_dataset(dataset,space,type,sessions='all'):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    if not isinstance(space, list):
        space = [space]
    if not isinstance(type, list):
        type = [type]

    for sp in space:
        if sp == 'MNISymC3':
            smooth=2.0
            interpolation=2
        else:
            smooth = None
            interpolation = 1 
        for t in type:
            if sessions=='all':
                sessions = mydataset.sessions
            for sess in sessions:
                print(f'extracting {dataset} type {t} space {sp}')
                mydataset.extract_all(ses_id=sess, type=t, atlas=sp,smooth=smooth,interpolation=interpolation)


if __name__ == "__main__":
    
    extract_dataset('Social', ['fs32k','MNISymC3'], ['Tseries'], ['ses-rest', 'ses-social'])
    extract_dataset('Language', ['fs32k','MNISymC3'], ['Tseries', 'FixTseries'], ['ses-rest', 'ses-localizer'])
    extract_dataset('MDTB', ['fs32k','MNISymC3'], ['Tseries', 'FixTseries'], ['ses-rest', 'ses-s1'])

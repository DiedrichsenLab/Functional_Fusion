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

def extract_dataset(dataset,space,type):
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
            for sess in mydataset.sessions:
                if sess !='ses-rest':
                    print(f'extracting {dataset} type {t} space {sp}')
                    mydataset.extract_all(ses_id=sess, type=t, atlas=sp,smooth=smooth,interpolation=interpolation)


if __name__ == "__main__":
    # datasets = ['WMFS','MDTB','Nishimoto','Somatotopic','IBC']
    # extract_dataset('MDTB', ['fs32k','MNISymC3'], ['CondAll','CondRun'])
    # extract_dataset('MDTB', ['fs32k','MNISymC3'], ['CondHalf'])
    #extract_dataset('WMFS', ['fs32k','MNISymC3'], ['CondHalf'])
    extract_dataset('HCPur100', ['fs32k','MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('Nishimoto', ['fs32k','MNISymC3'], ['CondHalf'])
    # extract_dataset('Somatotopic', ['fs32k','MNISymC3'], ['CondHalf'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # ['Social','Language','WMFS','MDTB','Demand','Nishimoto','Somatotopic','IBC']
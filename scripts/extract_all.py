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

def extract_dataset(dataset,atlas,type):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    if not isinstance(atlas, list):
        atlas = [atlas]
    if not isinstance(type, list):
        type = [type]

    for at in atlas:
        smooth = None
        interpolation = 1 
        for t in type:
            for sess in mydataset.sessions:
                if sess !='ses-rest':
                    print(f'extracting {dataset} type {t} space {at}')
                    mydataset.extract_all(ses_id=sess, type=t, atlas=at,smooth=smooth,interpolation=interpolation)

def extract_dataset_multiatlas(dataset,atlases,type,cifti_atlas_name=None):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    if not isinstance(type, list):
        type = [type]

    for t in type:
        smooth = None
        interpolation = 1 
        for sess in mydataset.sessions:
            if sess !='ses-rest':
                print(f'extracting {dataset} type {t} ')
                mydataset.extract_all(ses_id=sess, type=t, atlas=atlases,smooth=smooth,interpolation=interpolation,cifti_atlas_name=cifti_atlas_name)

def group_average(dataset,atlas,type):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    if not isinstance(atlas, list):
        atlas = [atlas]
    if not isinstance(type, list):
        type = [type]

    for at in atlas:
        for t in type:
            for sess in mydataset.sessions:
                print(f'group averaging {dataset} type {t} space {at}')
                mydataset.group_average_data(ses_id=sess, type=t, atlas=at)


if __name__ == "__main__":
    # datasets = ['WMFS','MDTB','Nishimoto','Somatotopic','IBC']
    # extract_dataset('MDTB', ['fs32k','MNISymC3'], ['CondAll','CondRun'])
    # extract_dataset('MDTB', ['fs32k','MNISymC3'], ['CondHalf'])
    #extract_dataset('WMFS', ['fs32k','MNISymC3'], ['CondHalf'])
    # extract_dataset('HCPur100', ['fs32k','MNISymC3'], ['CondAll','CondRun'])
    # extract_dataset('Nishimoto', ['fs32k','MNISymC3'], ['CondHalf'])
    # extract_dataset('Somatotopic', ['fs32k','MNISymC3'], ['CondHalf'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # extract_dataset('WMFS', ['MNISymC3'], ['CondHalf','CondAll','CondRun'])
    # ['Social','Language','WMFS','MDTB','Demand','Nishimoto','Somatotopic','IBC']
    # extract_dataset('MTLearn', ['fs32k'], ['CondAll'])
    # group_average('MTLearn', ['fs32k','MNISymC3'], ['CondAll'])
    for dataset in ['MDTB','Social','Language','WMFS','Demand','Nishimoto','Somatotopic','IBC']:
        extract_dataset_multiatlas('MDTB', ['MNIAsymHippocampus_L','MNIAsymHippocampus_R'], ['CondHalf'], cifti_atlas_name='MNIAsymHippocampus')
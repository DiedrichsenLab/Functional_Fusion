# Script to test multi-atlas extraction for different datasets

import pandas as pd
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as util

base_dir = util.get_base_dir()

def extract_dataset(dataset,atlas,multiatlas_name,type):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    for t in type:
        for sess in mydataset.sessions:
            if sess !='ses-rest':
                print(f'extracting {dataset} type {t} atlases {atlas} will be combined into into a single cifti called {multiatlas_name}')
                mydataset.extract_all(ses_id=sess, type=t, atlas=atlas,multiatlas_name=multiatlas_name,smooth=None,interpolation=1)


if __name__ == "__main__":
    hcp_atlases = [
    'fs32k_Asym_L',
    'fs32k_Asym_R',
    'MNIAsymAccumbens_L',
    'MNIAsymAccumbens_R',
    'MNIAsymAmygdala_L',
    'MNIAsymAmygdala_R',
    'MNIAsymBrainstem',
    'MNIAsymCaudate_L',
    'MNIAsymCaudate_R',
    'MNIAsymCerebellum_L',
    'MNIAsymCerebellum_R',
    'MNIAsymDiencephalonventral_L',
    'MNIAsymDiencephalonventral_R',
    'MNIAsymHippocampus_L',
    'MNIAsymHippocampus_R',
    'MNIAsymPallidum_L',
    'MNIAsymPallidum_R',
    'MNIAsymPutamen_L',
    'MNIAsymPutamen_R',
    'MNIAsymThalamus_L',
    'MNIAsymThalamus_R'
]
    extract_dataset('Language', hcp_atlases, 'multiatlasHCP', ['CondAll'])
    extract_dataset('MDTB', hcp_atlases, 'multiatlasHCP', ['CondAll'])
    extract_dataset('HCPur100', hcp_atlases, 'multiatlasHCP', ['CondAll'])
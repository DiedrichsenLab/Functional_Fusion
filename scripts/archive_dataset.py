'''Script to archive a dataset directory into a zip file for upload on Zenodo.
For convenience, the function can exclude time series data from the archive.
'''
import Functional_Fusion.util as ut 
import os 
import sys
import pandas as pd
zip_base = '/Users/jdiedrichsen/Data/Datasets'
base_dir = ut.get_base_dir()


def archive_subj(dataset, subj_name,sess,zip_dir):
    '''Archive a subject directory into a zip file for upload on Zenodo.
    Args:
        dataset: the dataset name
        subj_name: the subject name
        exclude_timeseries: whether to exclude time series data from the archive
    '''

    dataset_dir = os.path.join(base_dir, dataset,'derivatives','ffimport')

    if not os.path.exists(dataset_dir):
        raise ValueError('Subject directory does not exist: {}'.format(dataset_dir))

    os.system(f'')
    zip_file = os.path.join(zip_dir, f'{dataset}_{subj_name}.zip')

    patterns = []
    patterns.append(os.path.join(subj_name, 'anat',f'{subj_name}*'))
    for s in sess:
        patterns.append(os.path.join(subj_name, 'func',s, f'{subj_name}_{s}_run-*_reg*'))
        patterns.append(os.path.join(subj_name, 'func',s, f'{subj_name}_{s}_designmatrix*'))
        patterns.append(os.path.join(subj_name, 'func',s, f'{subj_name}_{s}_mask.nii'))
        patterns.append(os.path.join(subj_name, 'func',s, f'{subj_name}_{s}_reginfo.tsv'))
        patterns.append(os.path.join(subj_name, 'func',s, f'{subj_name}_{s}_resms.nii'))

    os.system(f'cd {dataset_dir};zip -r {zip_file} {" ".join(patterns)}')

def archive_dataset(dataset,sess):
    dataset_dir = os.path.join(base_dir, dataset)
    T = pd.read_csv(os.path.join(dataset_dir, 'participants.tsv'), sep='\t')
    zip_dir = os.path.join(zip_base, dataset,'derivatives','ffimport')
    if not os.path.exists(zip_dir):
        os.mkdir(zip_dir)

    for subj in T['participant_id']:
        archive_subj(dataset, subj, sess, zip_dir)



if __name__== '__main__':
    # archive_dataset('WMFS', ['ses-01', 'ses-02'])
    archive_dataset('MDTB',  ['ses-s1', 'ses-s2'])
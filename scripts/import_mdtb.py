# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np 

import import_data as id


base_dir = '/Volumes/diedrichsen_data$/data'
orig_dir = base_dir + '/Cerebellum/super_cerebellum'
target_dir = base_dir + '/FunctionalFusion/MDTB'





if __name__ == "__main__":
    T= pd.read_csv(target_dir + '/participants.tsv',delimiter='\t')
    for s in T.participant_id:
        old_id = s.replace('sub-','s',1)
        # dir1 = orig_dir + f'/sc1/suit/anatomicals/{old_id}'
        # dir2 = target_dir + f'/derivatives/{s}/suit'
        # id.import_suit(dir1,dir2,'anatomical',s)
        # dir1 = orig_dir + f'/sc1/anatomicals/{old_id}'
        # dir2 = target_dir + f'/derivatives/{s}/anat'
        # id.import_anat(dir1,dir2,'anatomical',s)
        # dir1 = orig_dir + f'/sc1/surfaceWB/{old_id}'
        # dir2 = target_dir + f'/derivatives/{s}/anat'
        # id.import_freesurfer(dir1,dir2,old_id,s)
        '''info_dict={'run':'run',
                   'inst':'instruction',
                   'TN':'task_name',
                   'CN':'cond_name',
                   'task':'task_num',
                   'cond':'cond_num'}
        '''
        dir1 = orig_dir + f'/sc1/GLM_firstlevel_7/{old_id}'
        dir2 = target_dir + f'/derivatives/{s}/estimates/ses-s1'
        id.import_spm_designmatrix(dir1,dir2,s,'ses-s1')
        pass
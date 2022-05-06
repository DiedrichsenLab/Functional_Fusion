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


def fix_sc2_reginfo(): 
    D = pd.read_csv(orig_dir + '/sc1_sc2_taskConds_conn.txt',sep='\t')
    T= pd.read_csv(target_dir + '/participants.tsv',delimiter='\t')
    for s in T.participant_id:
        for ses in [1,2]: 
            info = target_dir + f'/derivatives/{s}/estimates/ses-s{ses}/{s}_ses-s{ses}_reginfo.tsv'
            T1 = pd.read_csv(info,delimiter='\t')
            T1['task_num']=T1['task_num'].astype(int)
            T1['cond_num']=T1['cond_num'].astype(int)
            T1['instruction']=T1['instruction'].astype(int)
            T1['study'] = np.ones((T1.shape[0],1),dtype=int)*ses
            T1['cond_num_uni']=T1.cond_num
            T1['task_num_uni']=T1.task_num
            T1['common']=np.zeros((T1.shape[0],),dtype=int)

            D1 = D[D.StudyNum==ses]
            for i in range(T1.cond_num.max()): 
                indx = np.nonzero((D1.condNum==i+1).to_numpy())[0][0]
                cnu = D1.condNumUni.to_numpy()
                T1.cond_num_uni[T1.cond_num==i+1]=cnu[indx]
            for i in range(T1.task_num.max()): 
                indx = np.nonzero((D1.taskNum==i+1).to_numpy())[0][0]
                tnu = D1.taskNumUni.to_numpy()
                shc = D1.overlap.to_numpy()
                T1.task_num_uni[T1.task_num==i+1]=tnu[indx]
                T1.common[T1.task_num==i+1]=shc[indx]
            T1.to_csv(info,sep='\t')



if __name__ == "__main__":
    # fix_sc2_reginfo()
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
        print(s)
        info_dict={'run':'run',
                   'inst':'instruction',
                   'TN':'task_name',
                   'CN':'cond_name',
                   'task':'task_num',
                   'cond':'cond_num'}
        dir1 = orig_dir + f'/sc2/GLM_firstlevel_7/{old_id}'
        dir2 = target_dir + f'/derivatives/{s}/estimates/ses-s2'
        id.import_spm_glm(dir1,dir2,s,'ses-s2',info_dict)

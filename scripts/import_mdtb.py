# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import numpy as np
import Functional_Fusion.import_data as id
import scripts.paths as pt


base_dir = pt.set_base_dir()
atlas_dir = pt.set_atlas_dir(base_dir)

orig_dir = base_dir + '/../Cerebellum/super_cerebellum'
target_dir = base_dir + '/MDTB'


def fix_sc2_reginfo():
    D = pd.read_csv(orig_dir + '/sc1_sc2_taskConds_conn.txt', sep='\t')
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        for ses in [1, 2]:
            info = target_dir + \
                f'/derivatives/{s}/estimates/ses-s{ses}/{s}_ses-s{ses}_reginfo.tsv'
            T1 = pd.read_csv(info, delimiter='\t')
            T1 = T1[['run', 'instruction', 'task_name', 'cond_name',
                     'task_num', 'cond_num', 'reg_num', 'reg_id', 'study']]
            T1['task_num'] = T1['task_num'].astype(int)
            T1['cond_num'] = T1['cond_num'].astype(int)
            T1['instruction'] = T1['instruction'].astype(int)
            T1['study'] = np.ones((T1.shape[0], 1), dtype=int) * ses
            T1['cond_num_uni'] = T1.cond_num
            T1['task_num_uni'] = T1.task_num
            T1['common'] = np.zeros((T1.shape[0],), dtype=int)

            D1 = D[D.StudyNum == ses]
            for i in range(T1.cond_num.max()):
                indx = np.nonzero((D1.condNum == i + 1).to_numpy())[0][0]
                cnu = D1.condNumUni.to_numpy()
                T1.cond_num_uni[T1.cond_num == i + 1] = cnu[indx]
            for i in range(T1.task_num.max()):
                indx = np.nonzero((D1.taskNum == i + 1).to_numpy())[0][0]
                tnu = D1.taskNumUni.to_numpy()
                shc = D1.overlap.to_numpy()
                T1.task_num_uni[T1.task_num == i + 1] = tnu[indx]
                T1.common[T1.task_num == i + 1] = shc[indx]
            T1.to_csv(info, sep='\t', index=False)


if __name__ == "__main__":
    # fix_sc2_reginfo()
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    participants = T.participant_id

    # for s in participants:
    # old_id = s.replace('sub-','s',1)
    # dir1 = orig_dir + f'/sc1/suit/anatomicals/{old_id}'
    # dir2 = target_dir + f'/derivatives/{s}/suit'
    # id.import_suit(dir1,dir2,'anatomical',s)
    # dir1 = orig_dir + f'/sc1/anatomicals/{old_id}'
    # dir2 = target_dir + f'/derivatives/{s}/anat'
    # id.import_anat(dir1,dir2,'anatomical',s)
    # dir1 = orig_dir + f'/sc1/surfaceWB/{old_id}'
    # dir2 = target_dir + f'/derivatives/{s}/anat'
    # id.import_freesurfer(dir1,dir2,old_id,s)
    # print(s)
    # info_dict={'run':'run',
    #           'inst':'instruction',
    #           'TN':'task_name',
    #           'CN':'cond_name',
    #           'task':'task_num',
    #           'cond':'cond_num'}
    # dir1 = orig_dir + f'/sc1/GLM_firstlevel_7/{old_id}'
    # dir2 = target_dir + f'/derivatives/{s}/estimates/ses-s1'
    # id.import_spm_glm(dir1,dir2,s,'ses-s2',info_dict)
    # id.import_spm_designmatrix(dir1,dir2,s,'ses-s1')

    # Import resting-state session
    # (only take participants who have rest data)
    participants = participants[T['ses-rest'] == 1]
    for s in participants:
        old_id = s.replace('sub-', 's', 1)
        dir1 = orig_dir + '/resting_state/imaging_data_fix/'
        dir2 = target_dir + f'/derivatives/{s}/estimates/ses-rest'
        info_dict = {
            'runs': ['01', '02'],
            'reginfo_general': 'sub-02',
        }
        id.import_rest(dir1, dir2, s, 'ses-rest', info_dict)

    # T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    # for s in T.participant_id:
    #     print(f"-Start importing subject {s}")
    #     # old_id = s.replace('sub-','s',1)
    #     dir1 = os.path.join(orig_dir, str(s))
    #     dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
    #     import_func_resting(dir1, dir2, str(s))
    #     print(f"-Done subject {s}")

# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import numpy as np
import Functional_Fusion.import_data as id
import scripts.fusion_paths as paths
import shutil


base_dir = paths.set_base_dir()
atlas_dir = paths.set_atlas_dir(base_dir)

orig_dir = base_dir + '/Cerebellum/super_cerebellum'
target_dir = base_dir + 'FunctionalFusion/MDTB'


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

def import_fix_restingstate():
    # Import session s1 and session s2 fix-cleaned timeseries
    T = pd.read_csv(target_dir + '/FunctionalFusion/MDTB/participants.tsv', delimiter='\t')
    participants = T.participant_id
    participants_rest = participants[T['ses-rest'] == 1]

    runs = [f'{run:02d}' for run in np.arange(1, 17)]
    
    fix=False
    # for session in ['1','2']:
    for session in ['rest']:
        session_name = f'ses-s{session}' if session != 'rest' else 'ses-rest'
        dest_dir = base_dir + '/FunctionalFusion/MDTB/derivatives/{sub}/estimates/' + session_name + '/{sub}_' + session_name
        if fix:
            src_stem = base_dir + '/Cerebellum/super_cerebellum/sc1/' + 'imaging_data_fix/{sub}_' + session_name if session != 'rest' else base_dir + '/Cerebellum/super_cerebellum/resting_state/imaging_data_fix/'
            file_ending = '_run-{run}_fix.nii'
        else:
            src_stem = base_dir + '/Cerebellum/super_cerebellum/sc1//imaging_data/{sub}/' if session != 'rest' else base_dir + '/Cerebellum/super_cerebellum/resting_state/imaging_data/{sub}/'
        
        for s in participants:
            subject_orig = s.replace('sub-', 's', 1)
            src = src_stem.format(sub=subject_orig) + 'rrun_{run}.nii'
            dest = dest_dir.format(sub=s) + '_run-{run}.nii'
            mask_file = base_dir + '/Cerebellum/super_cerebellum/sc1/imaging_data_fix/{sub}_ses-s1_mask.nii'.format(sub=s) 
            id.import_tseries(src, dest, s, session_name, runs, mask_file=mask_file)


def import_mean_bold():
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    participants = T.participant_id
    for s in participants:
        subject_orig = s.replace('sub-', 's', 1)
        src = base_dir + f'/Cerebellum/super_cerebellum/sc1/imaging_data/{subject_orig}/rmeanepi.nii'
        dest = base_dir + f'/FunctionalFusion/MDTB/derivatives/{s}/anat/{s}_meanbold.nii'
        shutil.copyfile(src,dest)

if __name__ == "__main__":
    # fix_sc2_reginfo()
    # import_mean_bold()

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
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    participants = T[T['ses-rest'] == 1].participant_id
    for s in participants:
        old_id = s.replace('sub-', 's', 1)
        dir1 = orig_dir + '/resting_state/imaging_data_fix/'
        dir2 = target_dir + f'/derivatives/{s}/estimates/ses-rest'
        info_dict = {
            'runs': ['01', '02'],
            'reginfo_general': 'sub-02',
        }
        id.import_tseries(dir1, dir2, s, 'ses-rest', info_dict)

    # 
    # for s in T.participant_id:
    #     print(f"-Start importing subject {s}")
    #     # old_id = s.replace('sub-','s',1)
    #     dir1 = os.path.join(orig_dir, str(s))
    #     dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
    #     import_func_resting(dir1, dir2, str(s))
    #     print(f"-Done subject {s}")



# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import numpy as np
import Functional_Fusion.import_data as id
import scripts.fusion_paths as paths
import shutil
from pathlib import Path 
import nibabel as nb

base_dir = paths.set_base_dir()

orig_dir = base_dir + 'Cerebellum/LearningMultiTask_Maedbh'
target_dir = base_dir + 'FunctionalFusion_new/MTLearn'

def safe_copy(src,dest):
    # Extract the destination directory path
    dst_dir = Path(dest).parent

    # Create the directory structure if it doesn't exist
    dst_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src,dest)

def import_anatomical():
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for sub in T.participant_id:
        src = orig_dir + f'/anat/{sub}_desc-preproc_T1w_defaced.nii.gz'
        dest = target_dir + f'/derivatives/ffimport/{sub}/anat/{sub}_T1w.nii'
        safe_copy(src,dest)

def import_surface(): 
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    for sub in T.participant_id:
        src = orig_dir + f'/surfaces/{sub}'
        dest = target_dir + f'/derivatives/ffimport/{sub}/anat'
        id.import_freesurfer(src, dest, sub, sub)

def import_bold():
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    D = pd.read_csv(target_dir + '/reginfo.tsv',delimiter='\t')
    participants = T.participant_id
    for sub in participants:
        for ses in ['ses-1','ses-2','ses-3']:
            subj_dir = target_dir + f'/derivatives/ffimport/{sub}/func/{ses}'
            Path(subj_dir).mkdir(parents=True, exist_ok=True)

            # Average all _rmse.nii files across runs
            X = []
            for run in range(1, 7):
                rmse_file = nb.load(orig_dir + f'/func/{sub}/{sub}_glm-02_{ses}_run-{run:02d}_rmse.nii')
                X.append(rmse_file.get_fdata())
            mean_rmse = np.stack(X, axis=3).mean(axis=3)
            mean_file = nb.Nifti1Image(mean_rmse, rmse_file.affine, rmse_file.header)
            mean_file.to_filename(subj_dir + f'/{sub}_{ses}_resms.nii')

            # Save mask 
            X = mean_rmse >0 
            mask_file = nb.Nifti1Image(X.astype(np.uint8), rmse_file.affine)
            mask_file.to_filename(subj_dir + f'/{sub}_{ses}_mask.nii')

            # Break up the 4d niftis into 3d niftis for each run and save them as _run-xx_reg-xx.nii files
            reginfo=[]
            for run in range(1, 7):
                beta_file = nb.load(orig_dir + f'/func/{sub}/{sub}_glm-02_{ses}_run-{run:02d}_betas.nii')
                X = beta_file.get_fdata()
                if X.shape[3] != 16:
                    print(f'Warning: expected 16 regressors but found {X.shape[3]} for {sub} {ses} run {run}')
                for reg in range(16):
                    reg_file = nb.Nifti1Image(X[:,:,:,reg], beta_file.affine, beta_file.header)
                    reg_file.to_filename(subj_dir + f'/{sub}_{ses}_run-{run:02d}_reg-{reg:02d}.nii')
                d = pd.DataFrame({'run': [run]*16, 'task_name': D['task_name'], 
                                  'cond_name': D['cond_name'], 
                                  'reg_id': D['reg_id'],    
                                  'instruction': D['instruction'], 
                                  'task_code': D['task_code'], 
                                  'cond_code': D['cond_code']
                                  })
                reginfo.append(d)
            reginfo = pd.concat(reginfo, ignore_index=True)
            reginfo.to_csv(subj_dir + f'/{sub}_{ses}_reginfo.tsv', sep='\t', index=False)
            

if __name__ == "__main__":
    # fix_sc2_reginfo()
    # import_anatomical()
    # import_surface() 
    import_bold()   
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
    # participants = participants[T['ses-rest'] == 1]
    # for s in participants:
    #     old_id = s.replace('sub-', 's', 1)
    #     dir1 = orig_dir + '/resting_state/imaging_data_fix/'
    #     dir2 = target_dir + f'/derivatives/{s}/estimates/ses-rest'
    #     info_dict = {
    #         'runs': ['01', '02'],
    #         'reginfo_general': 'sub-02',
    #     }
    #     id.import_tseries(dir1, dir2, s, 'ses-rest', info_dict)

    # T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    # for s in T.participant_id:
    #     print(f"-Start importing subject {s}")
    #     # old_id = s.replace('sub-','s',1)
    #     dir1 = os.path.join(orig_dir, str(s))
    #     dir2 = os.path.join(target_dir, 'derivatives/%s/func' % str(s))
    #     import_func_resting(dir1, dir2, str(s))
    #     print(f"-Done subject {s}")



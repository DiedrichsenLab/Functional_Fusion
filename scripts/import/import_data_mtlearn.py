# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import numpy as np
import Functional_Fusion.import_data as id
import Functional_Fusion.util as ut
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
        dest = target_dir + f'/derivatives/ffimport/{sub}/anat/{sub}_T1w.nii.gz'
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
        print(f'Importing {sub}')
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
                    reg_file.to_filename(subj_dir + f'/{sub}_{ses}_run-{run:02d}_reg-{reg:02d}_beta.nii')
                d = pd.DataFrame({'run': [run]*16, 'task_name': D['task_name'], 
                                  'cond_name': D['cond_name'], 
                                  'reg_id': D['reg_id'],    
                                  'instruction': D['instruction'], 
                                  'task_code': D['task_code'], 
                                  'cond_code': D['cond_code'],
                                  'half': [np.mod(run-1,2)+1]*16
                                  })
                reginfo.append(d)
            reginfo = pd.concat(reginfo, ignore_index=True)
            reginfo.to_csv(subj_dir + f'/{sub}_{ses}_reginfo.tsv', sep='\t', index=False)
            
def run_suit(): 
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    participants = T.participant_id
    for sub in participants:
        print(f'SUITing {sub}')
        id.run_suit(target_dir + f'/derivatives/ffimport/{sub}/anat', sub,space='MNI200')

def rename_bold():
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    participants = T.participant_id
    for sub in participants:
        print(f'Importing {sub}')
        for ses in ['ses-1','ses-2','ses-3']:
            subj_dir = target_dir + f'/derivatives/ffimport/{sub}/func/{ses}'

            for run in range(1, 7):
                    for reg in range(16):
                        src = subj_dir + f'/{sub}_{ses}_run-{run:02d}_reg-{reg:02d}.nii'
                        trg = subj_dir + f'/{sub}_{ses}_run-{run:02d}_reg-{reg:02d}_beta.nii'
                        shutil.move(src, trg)

def rename_xfm():
    T = pd.read_csv(target_dir + '/participants.tsv', delimiter='\t')
    participants = T.participant_id
    for sub in participants:
        print(f'changing {sub}')
        subj_dir = target_dir + f'/derivatives/ffimport/{sub}/anat'
        src = subj_dir + f'/{sub}_space-MNISymC_xfm.nii.gz'
        trg = subj_dir + f'/{sub}_space-MNI152NLin2009cSymC_xfm.nii.gz'
        shutil.move(src, trg)





if __name__ == "__main__":
    # fix_sc2_reginfo()
    # import_anatomical()
    # import_surface() 
    # import_bold()   
    # run_suit()
    rename_xfm()
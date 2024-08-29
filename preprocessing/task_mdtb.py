import ProbabilisticParcellation.util as ut
import nibabel as nib
from pathlib import Path
import os.path as op
import subprocess
import random
from itertools import product
import pandas as pd
from datetime import datetime
import os
import shutil
import numpy as np
import scripts.fix as fx


data_dir = f'{ut.base_dir}/../Cerebellum/super_cerebellum/'
imaging_dir = Path(f'{data_dir}/sc1/imaging_data/')
fusion_dir = Path(f'{ut.base_dir}/MDTB/')
design_dir = Path(
    '~/code/Python/Functional_Fusion/preprocessing/design_files/').expanduser()
runs = np.arange(1, 17)
runs = [f'{run:02d}' for run in runs]
# For the MDTB dataset, task session 1 and 2 (sc1 & sc2) were concatenated, so the runs are from 1 to 32 when dealing with the original data files
runs_sessionscat = np.arange(1, 33)
runs_sessionscat = [f'{run:02d}' for run in runs_sessionscat]

sessions = ["1", "2"]
mdtb_subjects = pd.read_csv(f'{fusion_dir}/participants.tsv', delimiter='\t')
subjects = mdtb_subjects['participant_id']

# def copy_runs(fix=True):
#     """
#         Copies the raw runs into the estimates folder for each subject. Run this function instead of the import_rest function in import_data.py, because this function needs to additionally handle the issue of session 1 and 2 being concatenated in the MDTB dataset.
#         If fix is True, it will copy the FIX-cleaned runs, otherwise it will copy the raw (uncleaned) data.
        
#     """
#     T = pd.read_csv(f'{fusion_dir}/participants.tsv', delimiter='\t')
#     # --- Copy the raw runs into estimates ---
#     for subject in T.iterrows():
#         subject = subject[1].participant_id
#         for session in sessions:
#             imaging_folder = 'imaging_data_fix' if fix else 'imaging_data'
#             task_dir = Path(f'{data_dir}/sc1/{imaging_folder}/') if fix else Path(f'{data_dir}/sc1/{imaging_folder}/s{subject[-2:]}')
#             # Remove 'c' from session
#             for run in runs:
#                 file_name = f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_run-{run}_fix.nii" if fix else f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_run-{run}.nii"
#                 # if session is 1, then the run is the same as the original run, otherwise it is the original run + 16
#                 orig_run = int(run) if session == "1" else int(run) + 16
#                 task_file = f"{str(task_dir)}/{subject}_run-{orig_run:02d}.nii" if fix else f"{str(task_dir)}/rrun_{orig_run}.nii" 
#                 if op.exists(task_file) and not op.exists(file_name):
#                     subprocess.run(
#                         ['cp', task_file, file_name])
#                     print(f'Copied {run} for {subject} in {session}')

def get_labelled_folders():
    """
        Returns a list of the folders containing the labelled components for FIX training.
    """
    classified = pd.read_csv(f'{design_dir}/classifier_short.tsv', delimiter='\t')
    checked = classified[classified['checked'] == 'X']
    folder_pattern = '/{imaging_dir}/{subject}/run{run:02d}.feat/'
    labelled_folders = [folder_pattern.format(imaging_dir=imaging_dir, subject=line['subject'], run=line['run']) for _, line in checked.iterrows()]
    return labelled_folders


def copy_mean_func():
    """
    Copies the mean_func.nii.gz file from the filtered_func_data.ica folder to the parent folder for each subject and run.
    FIX needs the mean_func file in the main (feat) directory, not in the ica directory.
    """
    for i in imaging_dir.glob('s*/run*.feat/filtered_func_data.ica/mean.nii.gz'):
        print(" ")
        if not (i.parent / 'mean_func.nii.gz').exists():
            shutil.copy(i, i.parent / 'mean_func.nii.gz')


# def rename_run(imaging_dir, subject, run, new_run, session):
#     """
#     Renames the run folder from run to new_run for the given subject and session.
#     """
#     old_file = f"{imaging_dir}_fix/sc1/{subject}_run-{run}.nii"
#     new_file = f"{imaging_dir}_fix/sc1/{subject}_ses-s2_run-{new_run}.nii"
#     if op.exists(old_file) and not op.exists(new_file):
#         os.rename(old_file, new_file)


if __name__ == "__main__":
    # copy_runs()
    # make_tinfo_file()

    
    # ================================ FIX CLEANING ================================
    
    # --- Correct the header of the image files by inserting TR ---
    # for subject_path in imaging_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]  # remove the 's' prefix
    #     for run in runs_sessionscat:
    #         img_file = f"{str(subject_path)}/rrun_{run}.nii"
    #         fx.correct_header(img_file)

    # --- Create the design files for each subject and run single-subject ICA ---
    # for subject_path in imaging_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]
    #     for run in runs_sessionscat:
    #         # fx.make_design(subject, run, imaging_dir, design_dir, template_filestem='ssica_task')
    #         fx.run_ica(subject, run, imaging_dir, design_dir, template_filestem='ssica_task')

    # --- Copy motion parameter files to ica folders for feature extraction ---
    # for subject_path in imaging_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]
    #     for run in runs_sessionscat:
    #         fx.copy_motionparams(subject_path, run)

    # # --- Create a balanced subset of subjects and runs to classify into signal or noise ---
    # percent_data = 20
    # df = fx.make_classifier_sample(percent_data, imaging_dir, runs_sessionscat, outfile=f'{design_dir}/task_classifier_sample')

    # # --- Name the checked classification files 'hand_labels_noise.txt' ---
    # labelled_folders = get_labelled_folders()
    # for folder in labelled_folders:
    #     # Make copy of hand_labels_caro and save it as hand_labels_noise
    #     shutil.copy(f'{folder}/filtered_func_data.ica/hand_classification_caro', f'{folder}/hand_labels_noise.txt')

    # # --- Copy motion parameter files to ica folders for feature extraction ---
    # for subject in subjects:
    # for subject in ['sub-12', 'sub-24']:
    #     subject_orig = subject.replace('sub-', 's')
    #     subject_path = f'{imaging_dir}/{subject_orig}'
    #     for run in runs_sessionscat:
    #         fx.copy_motionparams(subject_path, run)
    
    # # --- Copy the mean_func.nii.gz file from the ica folder to the parent folder ---
    # copy_mean_func()

    # --- After classification, run fix training and leave-one-out testing ---
    # labelled_folders = get_labelled_folders()
    # os.chdir(f'{imaging_dir}/../fix_ica/')
    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-t', 'mdtb_task', '-l'] + labelled_folders)

    # # # --- Run leave-one-out testing using HCP training data and standard training data to compare acccuracy ---
    # labelled_folders = get_labelled_folders()
    
    # # change working directory to output directory (this is where the fix results will be saved)
    # os.chdir(f'{imaging_dir}/../fix_ica/')
    
    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-C', f'{imaging_dir}/../../resting_state/fix_ica/mdtb_rest.RData', 'mdtb_rest'] + labelled_folders)
    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-C', '/srv/software/fix/1.06.15/training_files/HCP_hp2000.RData', 'hcp3t'] + labelled_folders)
    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-C', '/srv/software/fix/1.06.15/training_files/Standard.RData', 'standard'] + labelled_folders)


    # --- Extract features for all ICAs that have been run so far - saves time for later fix-cleaning ---
    # for subject_path in imaging_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]
    #     for run in runs_sessionscat:
    #         ica_path = f"{str(subject_path)}/run{run}.feat/"
    #         if op.exists(ica_path) and not op.exists(f"{ica_path}/fix/features.csv"):
    #             subprocess.run(
    #                 ['/srv/software/fix/1.06.15/fix', '-f', ica_path])

    # --- Run FIX cleanup---
    # chosen_threshold = 20
    # # # For those scans that have hand-labelled components, clean noise components from the data
    # # labelled_folders = get_labelled_folders()
    # # for folder in labelled_folders:
    # #     subprocess.run(
    # #         ['/srv/software/fix/1.06.15/fix', '-a', f'{folder}/hand_labels_noise.txt'])
    
    # # For the rest, automatically classify labelled components using mdtb task training set, then clean noise components from the data
    # automatic_folders = [f"{folder}/run{run}.feat" for folder in imaging_dir.glob('s[0-9][0-9]') for run in runs_sessionscat if not op.exists(
    #     f'{folder}/run{run}.feat/filtered_func_data.ica/hand_labels_noise.txt')]
    # for subject in subjects:
    #     subject_orig = subject.replace('sub-', 's')
    #     folder = f'{imaging_dir}/{subject_orig}' 
    #     for run in runs_sessionscat:
    #         feat_path = f"{str(folder)}/run{run}.feat"
    #         if not op.exists(f"{feat_path}/filtered_func_data_clean.nii.gz"):
    #             print(f"Cleaning {subject} run {run}")
    #             subprocess.run(
    #                 ['/srv/software/fix/1.06.15/fix', '-c', feat_path, f'{imaging_dir}/../fix_ica/mdtb_task.RData', str(chosen_threshold)])
    #             subprocess.run(
    #                 ['/srv/software/fix/1.06.15/fix', '-a', f'{feat_path}/fix4melview_mdtb_task_thr{chosen_threshold}.txt'])
        
    # --- Move files ---
    # Move files to imaging_data_fix

    session = 's1'
    for subject in subjects:
        folder = f'{imaging_dir}/{subject}' 
        subject_orig = subject.replace('sub-', 's')
        # fx.move_mask(imaging_dir, subject_orig, session)
        for run in runs_sessionscat:
            
            # Rename runs 17-32 to 1-16 of session 2
            run = int(run)
            if run > 16:
                fx.move_cleaned(imaging_dir, subject_orig, run)
                # rename_run(imaging_dir, subject_orig, run, run-16, session='s2')

    


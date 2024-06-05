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
# For the MDTB dataset, task session 1 and 2 (sc1 & sc2) were concatenated, so the runs are from 1 to 32 when dealing with the original data files
runs_sessionscat = np.arange(1, 33)
# get zeropadded numbers from 1 to 32
runs = [f'{run:02d}' for run in runs]
sessions = ["1", "2"]

def copy_runs(fix=False):
    """
        Copies the raw runs into the estimates folder for each subject.
        If fix is True, it will copy the FIX-cleaned runs, otherwise it will copy the raw (uncleaned) data.
    """
    T = pd.read_csv(f'{fusion_dir}/participants.tsv', delimiter='\t')
    # --- Copy the raw runs into estimates ---
    for subject in T.iterrows():
        subject = subject[1].participant_id
        for session in sessions:
            imaging_folder = f'imaging_data_fix' if fix else f'imaging_data'
            task_dir = Path(f'{data_dir}/sc1/imaging_data/s{subject[-2:]}')
            # Remove 'c' from session
            for run in runs:
                # if session is 1, then the run is the same as the original run, otherwise it is the original run + 16
                orig_run = int(run) if session == "1" else int(run) + 16
                task_file = f"{str(task_dir)}/rrun_{orig_run}.nii"
                if op.exists(task_file):
                    subprocess.run(
                        ['cp', task_file, f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_run-{run}.nii"])
                    print(f'Copied {run} for {subject} in {session}')

def make_tinfo_file():
    """
        Creates info file for time series data in the following format:

        run	timepoint	task	time_id
        1	T0001	rest	1
        1	T0002	rest	2
        1	T0003	rest	3
        1	T0004	rest	4
        1	T0005	rest	5
        1	T0006	rest	6
        1	T0007	rest	7
    """

    T = pd.read_csv(f'{fusion_dir}/participants.tsv', delimiter='\t')
    for subject in T.iterrows():
        subject = subject[1].participant_id
        for session in sessions:
            tinfo_file = f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_tinfo.tsv"
            # Load last run
            run = runs[-1]
            img = nib.load(f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_run-{run}.nii")
            # Get the number of timepoints
            timepoints = img.shape[-1]

            with open(tinfo_file, 'w') as f:
                f.write('run\ttimepoint\ttask\ttime_id\n')
                for run in runs:
                    for timepoint in range(1, timepoints+1):
                        # make time_id continuous across runs
                        time_id = (int(run) - 1) * timepoints + timepoint
                        f.write(f"{int(run)}\tT{time_id:04d}\ttask\t{time_id}\n")
                
            print(f'Created {tinfo_file}')

if __name__ == "__main__":
    # copy_runs()
    # make_tinfo_file()

    # ================================ FIX CLEANING ================================
    
    # --- Create the design files for each subject and run single-subject ICA ---
    for subject_path in imaging_dir.glob('s[0-9][0-9]'):
        subject = subject_path.name[1:]
        for run in runs_sessionscat:
            fx.make_design(subject, run, imaging_dir, design_dir, template_filstem='ssica_task')
            # fx.run_ica(subject, run, imaging_dir, design_dir, template_filstem='ssica_task')

    # --- Copy motion parameter files to ica folders for feature extraction ---
    for subject_path in imaging_dir.glob('s[0-9][0-9]'):
        subject = subject_path.name[1:]
        for run in runs_sessionscat:
            fx.copy_motionparams(subject_path, run)

    # # --- After classification, run fix training and leave-one-out testing ---
    # labelled_folders = get_labelled_folders()
    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-t', 'mdtb_rest', '-l'] + labelled_folders)

    # # --- Run leave-one-out testing using HCP training data and standard training data to compare acccuracy ---
    # labelled_folders = get_labelled_folders()
    # # change working directory to output directory (this is where the fix results will be saved)
    # os.chdir(f'{rest_dir}/../fix_ica/')

    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-C', '/srv/software/fix/1.06.15/training_files/HCP_hp2000.RData', 'hcp3t'] + labelled_folders)
    # subprocess.run(
    #     ['/srv/software/fix/1.06.15/fix', '-C', '/srv/software/fix/1.06.15/training_files/Standard.RData', 'standard'] + labelled_folders)


    # --- Copy the FIX-cleaned runs into estimates ---
    # copy_runs(fix=True)

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
from Functional_Fusion.scripts.fix import copy_motionparams

data_dir = f'{ut.base_dir}/../Cerebellum/super_cerebellum/'
rest_dir = Path(f'{data_dir}/resting_state/imaging_data/')
anat_dir = Path(f'{data_dir}/sc1/anatomicals/')
design_dir = Path(
    '~/code/Python/Functional_Fusion/preprocessing/design_files/').expanduser()
runs = ["01", "02"]



def bet_anatomical(img_file, use='opti'):
    """Brain extract the anatomical image using FSL's BET.

    N.B. Using BBR registration for functional to anatomical (required for FIX classification) requires that the anatomical image be brain extracted.

    Args:
        img_file (string): path to the anatomical image file to be brain extracted
    """
    if use == 'opti':
        out_file = Path(f"{img_file.strip('.nii')}_optiBET_brain.nii.gz")
    elif use == 'bet':
        out_file = Path(f"{img_file.strip('.nii')}_brain.nii.gz")
    img_file = Path(img_file)

    if not out_file.exists() and img_file.exists():
        if use == 'opti':
            print(f"Brain extracting {img_file} with optiBET")

            # Brain extract the anatomical image
            subprocess.run(
                ['/srv/diedrichsen/shell/optiBET.sh', '-i', str(img_file)])

        elif use == 'bet':
            print(f"Brain extracting {img_file} with fsl")

            # Brain extract the anatomical image
            subprocess.run(['bet', str(img_file), str(out_file), '-R'])

    else:
        print(f"{img_file} already extracted")


def rename_anatomical(img_file):
    """Rename the anatomical image file to have a _brain.nii.gz suffix.

    Args:
        img_file (string): path to the anatomical image file
    """

    out_file = f"{img_file.strip('.nii')}_brain.nii.gz"
    opti_file = f"{img_file.strip('.nii')}_optiBET_brain.nii.gz"
    subprocess.run(['mv', str(opti_file), str(out_file)])

def classify_components():
    runs = ["01", "02"]
    # Load the list of subjects to classify
    classified_file = op.join(design_dir, 'classified_subjects.tsv')
    classified_df = pd.read_csv(classified_file, sep='\t', header=0)

    for subject_number, run in classified_df.values:
        subject_path = op.join(
            rest_dir, subject_number)

        # Check if components are not classified
        ica_path = op.join(
            subject_path, f'run{run:02d}.feat', 'filtered_func_data.ica')
        labels_file = op.join(ica_path, 'hand_labels_noise.txt')
        if op.exists(ica_path) and not op.exists(labels_file):
            print(f"Classifying {subject_number} run{run}")

            # Open motion parameter plots
            rp_file = op.join(subject_path, f'rp_run_{run}.txt')
            if op.exists(rp_file):
                if not op.exists(op.join(subject_path, f'rot_{run}.png')):
                    subprocess.run(['fsl_tsplot', '-i', rp_file, '-t', 'SPM estimated rotations (radians)',
                                    '-u', '1', '--start=1', '--finish=3', '-a', 'x,y,z', '-w', '800', '-h', '300',
                                    '-o', op.join(subject_path, f'rot_{run}.png')])
                    subprocess.run(['fsl_tsplot', '-i', rp_file, '-t', 'SPM estimated translations (mm)',
                                    '-u', '1', '--start=4', '--finish=6', '-a', 'x,y,z', '-w', '800', '-h', '300',
                                    '-o', op.join(subject_path, f'trans_{run}.png')])
                    subprocess.run(['fsl_tsplot', '-i', rp_file, '-t', 'SPM estimated rotations and translations (mm)',
                                    '-u', '1', '--start=1', '--finish=6', '-a', 'x(rot),y(rot),z(rot),x(trans),y(trans),z(trans)',
                                    '-w', '800', '-h', '300', '-o', op.join(subject_path, f'rot_trans_{run}.png')])

            # Open plots or report_prestats.html
            rot_trans_plot = op.join(subject_path, f'rot_trans_{run}.png')
            if op.exists(rot_trans_plot):
                subprocess.run(['open', op.join(subject_path, f'*_{run}.png')])
            else:
                subprocess.run(
                    ['open', op.join(subject_path, f'run{run}.feat', 'report_prestats.html')])

            # Open fsleyes
            subprocess.run(['fsleyes', '--scene', 'melodic', '-ad',
                           op.join(subject_path, f'run{run}.feat', 'filtered_func_data.ica')])

        else:
            print(f"Already classified {subject_number} run{run}")


def get_labelled_folders():
    labelled_folders = [f"{folder}/run{run}_smoothed.ica" for folder in rest_dir.glob('s[0-9][0-9]') for run in runs if op.exists(
        f'{folder}/run{run}_smoothed.ica/filtered_func_data.ica/hand_labels_noise.txt')]
    labelled_folders = labelled_folders + [f"{folder}/run{run}.feat" for folder in rest_dir.glob(
        's[0-9][0-9]') for run in runs if op.exists(f"{folder}/run{run}.feat/filtered_func_data.ica/hand_labels_noise.txt")]

    return labelled_folders



def move_cleaned(subject, run):
    """Move the cleaned data into the imaging_data_fix folder."""
    
    # move data into the corresponding session folder
    src = 'filtered_func_data_clean.nii.gz'
    dest = (f'/sub-{subject[1:]}_run-{run}.nii.gz')

    source_dir = f"{rest_dir}/{subject}/run{run}.feat/"
    dest_dir = f"{rest_dir}/../imaging_data_fix/"
    try:
        shutil.copyfile(source_dir + src,
                        dest_dir + dest)
        gunzip_cmd = f"gunzip {dest_dir + dest}"
        subprocess.run(gunzip_cmd, shell=True)
    except:
        print('skipping ' + src)


def move_mask(subject):    
    """Move the functional space grey matter mask into the imaging_data_fix folder."""
    src = 'rmask_noskull.nii'
    dest = (f'/sub-{subject[1:]}_ses-rest_mask.nii')

    source_dir = f"{rest_dir}/{subject}/"
    dest_dir = f"{rest_dir}/../imaging_data_fix/"
    try:
        shutil.copyfile(source_dir + src,
                        dest_dir + dest)
    except:
        print('skipping ' + src)

if __name__ == "__main__":

    # --- Brain-extract anatomical to be used in registration ---
    # for subject_path in anat_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]  # remove the 's' prefix
    #     for run in runs:
    #         img_file = f"{str(subject_path)}/anatomical.nii"
    #         bet_anatomical(img_file)
    #         rename_anatomical(img_file)

    # --- Correct the header of the image files by inserting TR ---
    # for subject_path in rest_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]  # remove the 's' prefix
    #     for run in runs:
    #         img_file = f"{str(subject_path)}/rrun_{run}.nii"
    #         correct_header(img_file)

    # --- Create the design files for each subject and run single-subject ICA ---
    # for subject_path in rest_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]
    #     for run in runs:
    #         make_design(subject, run)
    #         run_ica(subject, run)

    # # --- Create a balanced subset of subjects and runs to classify into signal or noise ---
    # make_classifier_sample()

    # --- Classify components ---
    # classify_components()

    # # --- Copy motion parameter files to ica folders for feature extraction ---
    # for subject_path in rest_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]
    #     for run in runs:
    #         copy_motionparams(subject_path, run)

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

    # --- Run FIX cleanup---
    # chosen_threshold = 20
    # # # For those scans that have hand-labelled components, clean noise components from the data
    # labelled_folders = [f"{folder}/run{run}.feat" for folder in rest_dir.glob('s[0-9][0-9]') for run in runs if op.exists(
    #     f'{folder}/run{run}.feat/filtered_func_data.ica/hand_labels_noise.txt')]
    # for folder in labelled_folders:
    #     subprocess.run(
    #         ['/srv/software/fix/1.06.15/fix', '-a', f'{folder}/hand_labels_noise.txt'])

    # For the rest, automatically classify labelled components using mdtb training set, then clean noise components from the data
    # automatic_folders = [f"{folder}/run{run}.feat" for folder in rest_dir.glob('s[0-9][0-9]') for run in runs if not op.exists(
    #     f'{folder}/run{run}.feat/filtered_func_data.ica/hand_labels_noise.txt')]
    # for folder in automatic_folders:
    #     # subprocess.run(
    #     #     ['/srv/software/fix/1.06.15/fix', '-c', folder, f'{str(rest_dir)}/../fix_ica/mdtb_rest.RData', str(chosen_threshold)])
    #     subprocess.run(
    #         ['/srv/software/fix/1.06.15/fix', '-a', f'{folder}/fix4melview_mdtb_rest_thr{chosen_threshold}.txt'])
        
    # --- Move files ---
    # Move files to imaging_data_fix
    folders = rest_dir.glob('s[0-9][0-9]')
    for folder in folders:
        folder = str(folder)
        subject = folder.split('/')[-1]
        move_mask(subject)
        # for run in runs:
        #     move_cleaned(subject, run)
        
    
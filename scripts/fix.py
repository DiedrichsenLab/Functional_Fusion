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


def correct_header(img_file):
    """Correct the header of the image file to have a TR of 1. Saves the file as a _hdr.nii.gz file for use with fsl.
    Args:
        img_file (string): path to the image file to be corrected
    """
    out_file = Path(f"{img_file.strip('.nii')}_hdr.nii.gz")
    img_file = Path(img_file)

    if not out_file.exists() and img_file.exists():
        print(f"Adding TR to header of {img_file}")

        img = nib.load(img_file)

        # Modify the TR field in the header
        img.header['pixdim'][4] = 1

        # Save the modified image as image file ending in '_hdr.nii.gz'
        nib.save(img, out_file)

    else:
        print(f"{img_file} already processed")


def make_design(subject, run, imaging_dir, design_dir, template_filstem='ssica'):
    """Create the design file for the single-subject ICA for the subject and run.

    Args:
        subject (string): subject ID
        run (string): run ID
    """

    img_file = Path(f"{imaging_dir}/s{subject}/rrun_{run}_hdr.nii.gz")
    design_template = Path(f"{design_dir}/{template_filstem}_template.fsf")
    design_output = Path(f"{design_dir}/{template_filstem}_{subject}_run-{run}.fsf")

    if img_file.is_file() and not design_output.is_file():
        # Read the contents of the template file
        with open(design_template, 'r') as template_file:
            design = template_file.read()

        # Replace placeholders in fsf content
        design = design.replace('XX', str(subject))
        design = design.replace('YY', str(run))

        # Write the modified content to the output file
        with open(design_output, 'w') as output_file:
            output_file.write(design)

    elif not img_file.is_file():
        print(f"{subject} {run}: missing image file")
    else:
        print(f"{subject} {run}: design file already created")


def run_ica(subject, run, imaging_dir, design_dir, template_filestem='ssica'):
    """Run the single-subject ICA on the resting state data for the subject and run.

    """

    img_file = Path(f"{imaging_dir}/s{subject}/rrun_{run}_hdr.nii.gz")
    feat_dir = Path(f"{imaging_dir}/s{subject}/run{run}.feat")
    ica_dir = Path(f"{imaging_dir}/s{subject}/run{run}.feat/filtered_func_data.ica")
    design_output = Path(f"{design_dir}/{template_filestem}_{subject}_run-{run}.fsf")

    if img_file.is_file() and not feat_dir.is_dir():
        print(f"Running ssica for subject {subject} run {run}")
        subprocess.Popen(['feat', str(design_output)])
    elif img_file.is_file() and feat_dir.is_dir() and not ica_dir.is_dir():
        print(f"Feat already run, but not completed. Re-running ssica for subject {subject} run {run}")
        subprocess.Popen(['mv', str(feat_dir), f"{feat_dir}_bak"])
        subprocess.Popen(['feat', str(design_output)])
    elif not img_file.is_file():
        print(f"{subject} {run}: missing image file")
    elif feat_dir.is_dir():
        print(f"{subject} {run}: ica already run")
        # use firefox if ut.base_dir.startswith('/Volumes') else 'open'
        # if ut.base_dir.startswith('/Volumes'):
        #     subprocess.run(['open', str(ica_dir / 'report.html')])
        # else:
        #     subprocess.run(['firefox', str(ica_dir / 'report.html')])

def copy_motionparams(subject_path, run):
    ica_path = f"{str(subject_path)}/run{run}.feat"
    rp_file = f"{str(subject_path)}/rp_run_{run}.txt"
    if op.exists(ica_path) and op.exists(rp_file):
        subprocess.run(['mkdir', f"{ica_path}/mc"])
        subprocess.run(
            ['cp', rp_file, f"{ica_path}/mc/prefiltered_func_data_mcf.par"])

def move_cleaned(imaging_dir, subject, run):
    """Move the cleaned data into the imaging_data_fix folder."""
    
    # move data into the corresponding session folder
    src = 'filtered_func_data_clean.nii.gz'
    dest = (f'/sub-{subject[1:]}_run-{run}.nii.gz')

    source_dir = f"{imaging_dir}/{subject}/run{run}.feat/"
    dest_dir = f"{imaging_dir}/../imaging_data_fix/"
    try:
        shutil.copyfile(source_dir + src,
                        dest_dir + dest)
        gunzip_cmd = f"gunzip {dest_dir + dest}"
        subprocess.run(gunzip_cmd, shell=True)
    except:
        print('skipping ' + src)



def move_mask(imaging_dir, subject):    
    """Move the functional space grey matter mask into the imaging_data_fix folder."""
    src = 'rmask_noskull.nii'
    dest = (f'/sub-{subject[1:]}_ses-rest_mask.nii')

    source_dir = f"{imaging_dir}/{subject}/"
    dest_dir = f"{imaging_dir}/../imaging_data_fix/"
    try:
        shutil.copyfile(source_dir + src,
                        dest_dir + dest)
    except:
        print('skipping ' + src)

# def balanced_subset(subjects, runs, percent_data):
#     # Calculate the number of subjects to select
#     num_scans = len(subjects) * len(runs)
#     num_scans_to_select = int(num_scans * percent_data * 0.01)
#     print(f"Selecting {num_scans_to_select} scans")

#     # Generate all combinations of subjects and runs
#     all_combinations = list(product(subjects, runs))

#     # Shuffle the combinations to ensure randomness
#     random.shuffle(all_combinations)

#     # Initialize counters
#     selected_subjects = set()
#     subset = []

#     # Select subjects while maintaining balance across runs
#     while len(selected_subjects) < num_scans_to_select:
#         for subject, run in all_combinations:
#             if num_scans_to_select < len(subjects):
#                 selected_subjects.add(subject)
#                 subset.append((subject, run))
#                 if subject not in selected_subjects:
#                     selected_subjects.add(subject)
#                     subset.append((subject, run))
#             else:
#                 selected_subjects.add(subject)
#                 subset.append((subject, run)) 
#             if len(selected_subjects) == num_scans_to_select:
#                 break   

#     # Separate the subset into lists of subjects and runs
#     subset_subjects, subset_runs = zip(*subset)

#     return list(subset_subjects), list(subset_runs)

def balanced_subset(subjects, runs, percent_data):
    """Create a balanced subset of subjects and runs to classify into signal or noise
    
    Args:
        subjects (list): list of subjects
        runs (list): list of runs
        percent_data (int): percentage of data to select
        
    Returns:
        list: list of subjects
        list: list of runs
    """
    # Calculate the number of subjects to select
    num_scans = len(subjects) * len(runs)
    num_scans_to_select = int(num_scans * percent_data * 0.01)
    print(f"Selecting {num_scans_to_select} scans")

    # Generate all combinations of subjects and runs
    all_combinations = list(product(subjects, runs))

    # Shuffle the combinations to ensure randomness
    random.shuffle(all_combinations)

    # Initialize the subset
    subset = []

    # Select scans while maintaining balance
    while len(subset) < num_scans_to_select:
        for subject, run in all_combinations:
            if len(subset) < num_scans_to_select:
                subset.append((subject, run))
            else:
                break

    # Separate the subset into lists of subjects and runs
    subset_subjects, subset_runs = zip(*subset) if subset else ([], [])

    return list(subset_subjects), list(subset_runs)

def make_classifier_sample(percent_data, imaging_dir, runs, outfile, add_new_subjects=False):
    """Create a balanced subset of subjects and runs to classify into signal or noise

    Args:
        percent_data (int): percentage of data to select
        imaging_dir (Path): path to the imaging data
        runs (list): list of runs
        outfile (Path): path to the output file
        add_new_subjects (bool): whether to add new subjects to the existing dataframe
    """
    # get first element of subject folders
    subject_list = [subject.name for subject in imaging_dir.glob('s[0-9][0-9]')]
    subset_subjects, subset_runs = balanced_subset(
        subject_list, runs, percent_data)
    # Save
    df = pd.DataFrame({'subject': subset_subjects, 'run': subset_runs})
    df = df.sort_values(by=['subject', 'run'])
    # if file already exists, add a timestamp to the filename
    outfile = Path(outfile)
    if outfile.is_file():
        outfile = Path(
            f"{outfile}_{datetime.now().strftime('%Y%m%d')}.tsv")

    if outfile.is_file() and add_new_subjects:
        # Add the new subjects to the existing dataframe
        akready_classified = pd.read_csv(
            Path(f"{outfile}.tsv"), sep='\t')
        akready_classified.run = akready_classified.run.astype(int)
        df = pd.concat([akready_classified, df], ignore_index=True)
        # Remove duplicates
        df = df.drop_duplicates()
        df = df.drop_duplicates(subset=['subject'], keep='first').sort_values(
            by=['subject', 'run'])
        # Check if  the subset is still balanced across runs
        print(f"Runs : {df.groupby(['run']).size()}")
        # Make run into zero-padded string
        df.run = df.run.apply(lambda x: f"{x:02d}")

    # Save
    df.to_csv(outfile, sep='\t', index=False)

    return df
import ProbabilisticParcellation.util as ut
import nibabel as nib
from pathlib import Path
import subprocess
import random
from itertools import product
import pandas as pd
from datetime import datetime


data_dir = f'{ut.base_dir}/../Cerebellum/super_cerebellum/'
rest_dir = Path(f'{data_dir}/resting_state/imaging_data/')
anat_dir = Path(f'{data_dir}/sc1/anatomicals/')
design_dir = Path(
    '~/code/Python/Functional_Fusion/preprocessing/design_files/').expanduser()
runs = ["01", "02"]


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
            subprocess.run(['/srv/diedrichsen/shell/optiBET.sh', '-i', str(img_file)])

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
    new_file = f"{img_file.strip('.nii')}_optiBET_brain.nii.gz"
    subprocess.run(['mv', str(out_file), str(new_file)])


def make_design(subject, run):
    """Create the design file for the single-subject ICA for the subject and run.

    Args:
        subject (string): subject ID
        run (string): run ID
    """

    img_file = Path(f"{rest_dir}/s{subject}/rrun_{run}_hdr.nii.gz")
    design_template = Path(f"{design_dir}/ssice_template.fsf")
    design_output = Path(f"{design_dir}/ssice_{subject}_run-{run}.fsf")

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


def run_ica(subject, run):
    """Run the single-subject ICA on the resting state data for the subject and run.

    """

    img_file = Path(f"{rest_dir}/s{subject}/rrun_{run}_hdr.nii.gz")
    ica_dir = Path(f"{rest_dir}/s{subject}/run{run}.ica")
    design_output = Path(f"{design_dir}/ssica_{subject}_run-{run}.fsf")

    if img_file.is_file():
        print(f"Running ssica for subject {subject} run {run}")
        subprocess.Popen(['feat', str(design_output)])

    elif not img_file.is_file():
        print(f"{subject} {run}: missing image file")
    else:
        print(f"{subject} {run}: ica already run")
        # use firefox if ut.base_dir.startswith('/Volumes') else 'open'
        if ut.base_dir.startswith('/Volumes'):
            subprocess.run(['open', str(ica_dir / 'report.html')])
        else:
            subprocess.run(['firefox', str(ica_dir / 'report.html')])


def balanced_subset(subjects, runs, percent_data):
    # Calculate the number of subjects to select
    num_scans = len(subjects) * len(runs)
    num_scans_to_select = int(num_scans * percent_data * 0.01)

    # Generate all combinations of subjects and runs
    all_combinations = list(product(subjects, runs))

    # Shuffle the combinations to ensure randomness
    random.shuffle(all_combinations)

    # Initialize counters
    selected_subjects = set()
    subset = []

    # Select subjects while maintaining balance across runs
    for subject, run in all_combinations:
        if subject not in selected_subjects:
            selected_subjects.add(subject)
            subset.append((subject, run))

        if len(selected_subjects) == num_scans_to_select:
            break

    # Separate the subset into lists of subjects and runs
    subset_subjects, subset_runs = zip(*subset)

    return list(subset_subjects), list(subset_runs)


def make_classifier_sample(add_new_subjects=False):
    # Create a balanced subset of subjects and runs to classify into signal or noise
    # get first element of subject folders
    subject_list = [subject.name for subject in rest_dir.glob('s[0-9][0-9]')]
    percent_data = 30
    subset_subjects, subset_runs = balanced_subset(
        subject_list, runs, percent_data)
    # Save
    df = pd.DataFrame({'subject': subset_subjects, 'run': subset_runs})
    df = df.sort_values(by=['subject', 'run'])
    # if file already exists, add a timestamp to the filename
    subset_file = Path(
        f"{design_dir}/classified_subjects_{datetime.now().strftime('%Y%m%d')}.tsv")

    if add_new_subjects:
        # Add the new subjects to the existing dataframe
        akready_classified = pd.read_csv(
            Path(f"{design_dir}/classified_subjects.tsv"), sep='\t')
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
    df.to_csv(subset_file, index=False, sep='\t')


if __name__ == "__main__":

    # --- Brain-extract anatomical to be used in registration ---
    for subject_path in anat_dir.glob('s[0-9][0-9]'):
        subject = subject_path.name[1:]  # remove the 's' prefix
        for run in runs:
            img_file = f"{str(subject_path)}/anatomical.nii"
            bet_anatomical(img_file)
            rename_anatomical(img_file)
    

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

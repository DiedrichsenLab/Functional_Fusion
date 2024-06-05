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


def make_design(subject, run, imaging_dir, design_dir, template_filstem='ssica'):
    """Create the design file for the single-subject ICA for the subject and run.

    Args:
        subject (string): subject ID
        run (string): run ID
    """

    img_file = Path(f"{imaging_dir}/s{subject}/rrun_{run}_hdr.nii.gz")
    design_template = Path(f"{design_dir}/{template_file}_template.fsf")
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
    ica_dir = Path(f"{imaging_dir}/s{subject}/run{run}.ica")
    design_output = Path(f"{design_dir}/{template_filestem}_{subject}_run-{run}.fsf")

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

def copy_motionparams(subject_path, run):
    ica_path = f"{str(subject_path)}/run{run}.feat"
    rp_file = f"{str(subject_path)}/rp_run_{run}.txt"
    if op.exists(ica_path) and op.exists(rp_file):
        subprocess.run(['mkdir', f"{ica_path}/mc"])
        subprocess.run(
            ['cp', rp_file, f"{ica_path}/mc/prefiltered_func_data_mcf.par"])

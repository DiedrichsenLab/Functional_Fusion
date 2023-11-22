import ProbabilisticParcellation.util as ut
import nibabel as nib
from pathlib import Path
import subprocess

data_dir = Path(f'{ut.base_dir}/../Cerebellum/super_cerebellum/resting_state/imaging_data/')
design_dir = Path('~/code/Python/Functional_Fusion/preprocessing/design_files/').expanduser()
runs = ["01", "02"]

def correct_header(img_file):
    """Correct the header of the image file to have a TR of 1. Saves the file as a _hdr.nii.gz file for use with fsl.
    Args:
        img_file (string): path to the image file to be corrected
    """
    out_file = Path(f"{img_file}.gz")
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


def make_design(subject, run):
    img_file = Path(f"{str(subject_path)}/rrun_{run}_hdr.nii.gz")
    design_template = Path(f"{design_dir}/ssica_template.fsf")
    design_output = Path(f"{design_dir}/rest_{subject}_run-{run}.fsf")

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

    img_file = Path(f"{str(subject_path)}/rrun_{run}_hdr.nii.gz")
    ica_dir = Path(f"{subject_path}/run{run}.ica")
    design_template = Path(f"{design_dir}/ssica_template.fsf")
    design_output = Path(f"{design_dir}/rest_{subject}_run-{run}.fsf")

    if img_file.is_file() and not ica_dir.is_dir():
        print(f"Running SS melodic for subject {subject} run {run}")
        subprocess.run(['feat', str(design_output)])

    elif not img_file.is_file():
        print(f"{subject} {run}: missing image file")
    else:
        print(f"{subject} {run}: ica already run")
        # use firefox if ut.base_dir.startswith('/Volumes') else 'open'
        if ut.base_dir.startswith('/Volumes'):
            subprocess.run(['open', str(ica_dir / 'report.html')])
        else:
            subprocess.run(['firefox', str(ica_dir / 'report.html')])


if __name__ == "__main__":

    # for subject_path in data_dir.glob('s[0-9][0-9]'):
    #     subject = subject_path.name[1:]  # remove the 's' prefix
    #     for run in runs:
    #         img_file = f"{str(subject_path)}/rrun_{run}.nii"
    #         correct_header(img_file)

    for subject_path in data_dir.glob('s[0-9][0-9]'):
        subject = subject_path.name[1:]
        for run in runs:
            make_design(subject, run)
            run_ica(subject, run)

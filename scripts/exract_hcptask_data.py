import os, subprocess, re
import zipfile
import shutil
import nibabel as nb
import pandas as pd
from pathlib import Path
from Functional_Fusion.import_data import *
import numpy as np
import gzip
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetHcpTask

# directory = '/data/tge/dzhi/projects/HCP_tfMRI'
directory = '/data/tge/Tian/BANDA/BANDA_U01'
if not os.path.exists(directory):
    directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_new'

functional_fusion_dir = f'/data/tge/Tian/BANDA'
HCP_dir = f'/data/tge/Tian/BANDA/BANDA_U01'

# step 3: Reorganize anat stuff
def process_anat(directory):
    for participant_id in os.listdir(directory):
        # Original folder path
        participant_path = os.path.join(directory, participant_id)
        
        if os.path.isdir(participant_path):
            participant_number = re.search(r'\d+', participant_id).group()
            output_path = os.path.join(directory, f'sub-{participant_number}')
            anat_folder = os.path.join(output_path, 'anat')
            os.makedirs(anat_folder, exist_ok=True)

            mni_folder = os.path.join(participant_path, 'MNINonLinear')
            
            if os.path.isdir(mni_folder):
                # Copy the 'xfms' folder if it exists and isn't already in 'anat'
                xfms_folder = os.path.join(mni_folder, 'xfms')
                if os.path.isdir(xfms_folder) and not os.path.exists(os.path.join(anat_folder, 'xfms')):
                    shutil.copytree(xfms_folder, os.path.join(anat_folder, 'xfms'))

                # Decompress and copy 'T1w.nii.gz' and 'BiasField.nii.gz' into the 'anat' folder
                t1w_file = os.path.join(mni_folder, 'T1w.nii.gz')
                bias_field_file = os.path.join(mni_folder, 'BiasField.nii.gz')
                shutil.copy(t1w_file, os.path.join(anat_folder, 'T1w.nii.gz'))
                shutil.copy(bias_field_file, os.path.join(anat_folder, 'BiasField.nii.gz'))

                # Move and rename 'fsaverage_LR32k' to 'SurfaceWB' change this name?
                fsaverage_folder = os.path.join(mni_folder, 'fsaverage_LR32k')
                surface_wb_folder = os.path.join(output_path, 'SurfaceWB')

                if os.path.isdir(fsaverage_folder):
                    shutil.copytree(fsaverage_folder, surface_wb_folder)

            # Remove the structural preproc folder
            # shutil.rmtree(mni_folder, ignore_errors=True)

        
#  step 4: Reorganize func stuff
def process_func(directory):
    for subject_id in os.listdir(directory):  # Iterate over subject directories
        subject_path = os.path.join(directory, subject_id)
        subject_number = re.search(r'\d+', subject_id).group()
        
        if os.path.isdir(subject_path):
            output_path = os.path.join(directory, f'sub-{subject_number}')
            func_directory = os.path.join(output_path, 'func')
            os.makedirs(func_directory, exist_ok=True)

            session_path = os.path.join(subject_path, 'MNINonLinear', 'Results')
            # Iterate over session directories that need preprocessing
            for run in os.listdir(session_path):
                # Process directories with '3T_tfMRI' in the name
                if ('_AP' in run) or ('_PA' in run):
                    session_name = run.split('_')[-2]
                    this_run_path = os.path.join(session_path, run)
                                        
                    # Create a session directory inside 'func'
                    session_func_dir = os.path.join(func_directory, f'ses-{session_name}') 
                    os.makedirs(session_func_dir, exist_ok=True)
                    run_func_dir = os.path.join(session_func_dir, run)

                    # Move the run folder to the session directory in 'func'
                    shutil.copytree(this_run_path, run_func_dir, dirs_exist_ok=True)
                    
                    if os.path.isdir(run_func_dir):
                        # Clean unwanted files
                        for file in os.listdir(run_func_dir):
                            file_path = os.path.join(run_func_dir, file)
                            if file.endswith((
                                '.dtseries.nii', '.func.gii', 'SBRef_dc.nii.gz',
                                'gdc_dc.nii.gz', 'MSMAll.dtseries.nii'
                            )) or 'RibbonVolumeToSurfaceMapping' in file_path:
                                if os.path.isfile(file_path):
                                    os.remove(file_path)
                                elif os.path.isdir(file_path):
                                    shutil.rmtree(file_path, ignore_errors=True)

                    print(f'Done cleaning subject {subject_number} {run}')

                # shutil.rmtree(os.path.join(subject_path, session), ignore_errors=True)


# step 1: change the smoothing kernel in the fsf files
def update_fsf_smooth(file_path):
    """
    Update the `set fmri(smooth)` value in the specified .fsf file to 0.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        updated_lines = []
        for line in lines:
            if line.startswith("set fmri(smooth)"):
                # Change the value to 0
                updated_lines.append("set fmri(smooth) 0\n")
            elif line.startswith("set fmri(analysis)"):
                updated_lines.append("set fmri(analysis) 2\n")
            else:
                updated_lines.append(line)
        
        # Overwrite the file with updated content
        with open(file_path, 'w') as file:
            file.writelines(updated_lines)
        print(f"Updated: {file_path}")
    except Exception as e:
        print(f"Error processing {file_path}: {e}")

def find_and_update_fsf_files(base_dir):
    """
    Walk through the directory structure starting from base_dir, look for .fsf files,
    and update the `set fmri(smooth)` value in each file.
    """
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".fsf"):
                fsf_path = os.path.join(root, file)
                update_fsf_smooth(fsf_path)


def update_outputdir_and_run_feat(base_dir):
    """
    Updates the output directory in each .fsf file to match its location
    and runs FEAT on the updated file.
    """
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith("level1.fsf"):
                fsf_path = os.path.join(root, file)
                output_dir = os.path.join(root, f"{os.path.splitext(file)[0]}.feat")
                
                # Update the `set fmri(outputdir)` in the .fsf file
                updated_lines = []
                with open(fsf_path, 'r') as f:
                    for line in f:
                        if line.startswith("set fmri(outputdir)"):
                            updated_lines.append(f'set fmri(outputdir) "{output_dir}"\n')
                        else:
                            updated_lines.append(line)
                
                with open(fsf_path, 'w') as f:
                    f.writelines(updated_lines)
                
                print(f"Updated output directory in: {fsf_path}")
                
                # Run FEAT
                try:
                    print(f"Running FEAT on: {fsf_path}")
                    subprocess.run(["feat", fsf_path], check=True)
                    print(f"FEAT successfully completed for: {fsf_path}")
                except subprocess.CalledProcessError as e:
                    print(f"Error running FEAT on {fsf_path}: {e}")
                except Exception as e:
                    print(f"Unexpected error: {e}")


def make_participant_tsv(source_dir, dest_dir):
    if not Path(dest_dir).exists():
        os.makedirs(dest_dir)

    dest_dir = Path(dest_dir)
    source_dir = Path(source_dir)
    
    subj_list = []
    for subj in Path(source_dir).iterdir():
        if subj.is_dir() and subj.name.startswith("sub-"):
            subj_list.append(subj.name)
    pd.DataFrame(subj_list, columns=["participant_id"]).to_csv(dest_dir / "participants.tsv", sep="\t", index=False)

    return


def import_anat_data(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        anat_file = f'{source_dir}/{participant}/anat/T1w.nii.gz'
        dest_folder = f'{dest_dir}/derivatives/{participant}/anat'
        dest_file = f'{dest_folder}/{participant}_T1w.nii.gz'

        if not Path(dest_folder).exists():
            os.makedirs(dest_folder, exist_ok=True)
        
        # copy anat file
        shutil.copyfile(anat_file, dest_file)
        print(f'Copied {anat_file} to {dest_file} for participant {participant}')

    return


def import_freesurfer(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        participant_number = 'BANDA' + participant.split("-")[1] + '_MR'
        pial_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.L.pial.32k_fs_LR.surf.gii'
        pial_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.R.pial.32k_fs_LR.surf.gii'
        white_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.L.white.32k_fs_LR.surf.gii'
        white_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.R.white.32k_fs_LR.surf.gii'
        sulc_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.L.sulc.32k_fs_LR.shape.gii'
        sulc_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.R.sulc.32k_fs_LR.shape.gii'

        pial_L_dest = f'{dest_dir}/derivatives/{participant}/anat/{participant}_space-32k_hemi-L_pial.surf.gii'
        pial_R_dest = f'{dest_dir}/derivatives/{participant}/anat/{participant}_space-32k_hemi-R_pial.surf.gii'
        white_L_dest = f'{dest_dir}/derivatives/{participant}/anat/{participant}_space-32k_hemi-L_white.surf.gii'
        white_R_dest = f'{dest_dir}/derivatives/{participant}/anat/{participant}_space-32k_hemi-R_white.surf.gii'
        sulc_L_dest = f'{dest_dir}/derivatives/{participant}/anat/{participant}_space-32k_hemi-L_sulc.shape.gii'
        sulc_R_dest = f'{dest_dir}/derivatives/{participant}/anat/{participant}_space-32k_hemi-R_sulc.shape.gii'

        # copy all files
        shutil.copyfile(pial_L_source, pial_L_dest)
        shutil.copyfile(pial_R_source, pial_R_dest)
        shutil.copyfile(white_L_source, white_L_dest)
        shutil.copyfile(white_R_source, white_R_dest)
        shutil.copyfile(sulc_L_source, sulc_L_dest)
        shutil.copyfile(sulc_R_source, sulc_R_dest)

        print(f'Copied freesurfer files for participant {participant}')

    return


def import_resms(source_dir,dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        resms_files = list(Path(f'{source_dir}/{participant}').rglob('sigmasquareds.nii.gz'))
        dest_folder = f'{dest_dir}/derivatives/{participant}/estimates/ses-task'
        dest_file = f'{dest_folder}/{participant}_ses-task_resms.nii'
        if not Path(dest_folder).exists():
            os.makedirs(dest_folder, exist_ok=True)
        resms_data = []
        for resms_file in resms_files:
            img = nb.load(resms_file)
            data = img.get_fdata()
            resms_data.append(data)
        resms_data = np.mean(resms_data, axis=0)
        resms_img = nb.Nifti1Image(resms_data, img.affine)
        nb.save(resms_img, dest_file)
        print(f'Copied resms file to {dest_file} for participant {participant}')

    return

def import_betas(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        # Gather all task session directories
        participant_dir = Path(source_dir) / participant / 'func'
        session_dirs = [d for d in participant_dir.iterdir() if d.is_dir() 
                        and d.name.startswith("ses-") and not 'REST' in d.name]

        session_run_mapping = {}

        for session_index, session_dir in enumerate(session_dirs):
            run_dirs = [d for d in session_dir.iterdir() if d.is_dir()]
            session_run_mapping[session_index] = run_dirs

        max_runs = max(len(runs) for runs in session_run_mapping.values())
        global_run_counter = 1  # Global counter for runs across sessions

        for run_index in range(max_runs):
            reg_num = 1
            for session_index, run_dirs in session_run_mapping.items():
                if run_index < len(run_dirs):
                    run_dir = run_dirs[run_index]

                    # Path to the stats folder
                    feat_dir = list(run_dir.glob("*.feat"))
                    if not feat_dir:
                        continue
                    stats_dir = feat_dir[1] / "stats"

                    if stats_dir.exists():
                        pe_files = list(stats_dir.glob("pe*.nii.gz"))

                        # Filter for odd-numbered PE files (derivatives)
                        odd_pe_files = []
                        for pe in pe_files:
                            beta_number = int(pe.name.replace(".nii.gz", "").replace("pe", ""))
                            if beta_number % 2 == 1:
                                odd_pe_files.append(pe)

                        # Enumerate odd PE files
                        for  pe_file in odd_pe_files:
                            reg_label = f"reg-{reg_num:02d}"

                            # Construct output filename
                            dest_folder = Path(dest_dir) / "derivatives" / participant / "estimates" / "ses-task"
                            dest_file = dest_folder / f"{participant}_ses-task_run-{global_run_counter:02d}_{reg_label}_beta.nii.gz"

                            if not dest_folder.exists():
                                os.makedirs(dest_folder, exist_ok=True)

                            # Copy beta file
                            shutil.copyfile(pe_file, dest_file)
                            decompressed_file = dest_file.with_suffix('')
                            with gzip.open(dest_file, 'rb') as f_in:
                                with open(decompressed_file, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)

                            # Remove the original .gz file
                            os.remove(dest_file)

                            reg_num += 1
                global_run_counter += 1
        print(f"Copied beta files for participant {participant}")

    return

def import_masks(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        mask_files = list(Path(f'{source_dir}/{participant}').rglob('mask.nii.gz'))

        dest_folder = f'{dest_dir}/derivatives/{participant}/estimates/ses-task'
        dest_file = f'{dest_folder}/{participant}_ses-task_mask.nii'
        if not Path(dest_folder).exists():
            os.makedirs(dest_folder, exist_ok=True)
        
        mask_data = []
        for mask_file in mask_files:
            img = nb.load(mask_file)
            data = img.get_fdata()
            mask_data.append(data)
        
        # Average the mask data
        mask_data = np.mean(mask_data, axis=0)
        
        mask_img = nb.Nifti1Image(mask_data, img.affine)
        nb.save(mask_img, dest_file)
        print(f"Copied mask file to {dest_file} for participant {participant}")

    return


def make_reginfo(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        reginfo_data = []  # To store reginfo entries
        participant_dir = Path(source_dir) / participant / "func"

        # Get session directories
        session_dirs = [d for d in participant_dir.iterdir() if d.is_dir() and d.name.startswith("ses-")]

        session_run_mapping = {}
        for session_index, session_dir in enumerate(session_dirs):
            run_dirs = [d for d in session_dir.iterdir() if d.is_dir()]
            session_run_mapping[session_index] = run_dirs

        max_runs = max(len(runs) for runs in session_run_mapping.values())
        global_run_counter = 1  # Global counter for runs

        for run_index in range(max_runs):
            reg_id = 1  # Global regressor ID for each run
            for session_index, run_dirs in session_run_mapping.items():
                if run_index < len(run_dirs):
                    run_dir = run_dirs[run_index]

                    # Path to the stats folder and design.fsf
                    feat_dir = list(run_dir.glob("*.feat"))
                    if not feat_dir:
                        continue
                    design_fsf = feat_dir[0] / "design.fsf"

                    # Parse design.fsf to extract condition names
                    with open(design_fsf, "r") as fsf_file:
                        lines = fsf_file.readlines()

                    conditions = []
                    for line in lines:
                        if line.startswith("set fmri(evtitle"):
                            cond_name = line.split('"')[1]  # Extract condition name
                            conditions.append(cond_name)

                    # Populate reginfo entries for this run
                    reg_num = 1  # Local regressor ID for each run
                    for cond_name in conditions:
                        reginfo_data.append({
                            "run": global_run_counter,
                            "task_name": feat_dir[0].stem.split('_')[1],
                            "cond_name": cond_name,
                            "reg_id": reg_id,
                            "reg_num": reg_num,  # Increment regressor number for conditions within run
                            "half": 1 if run_index % 2 == 0 else 2
                        })
                        reg_num += 1
                        reg_id += 1

                    # Increment the global run counter after processing each run
                    global_run_counter += 1

        # Save reginfo.tsv
        dest_folder = Path(dest_dir) / "derivatives" / participant / "estimates" / "ses-task"
        dest_file = dest_folder / f"{participant}_ses-task_reginfo.tsv"
        if not dest_folder.exists():
            os.makedirs(dest_folder, exist_ok=True)

        reginfo_df = pd.DataFrame(reginfo_data, columns=[ "run", "task_name", "cond_name", "reg_id","reg_num", "half"])
        reginfo_df.to_csv(dest_file, sep="\t", index=False)
        print(f"Saved reginfo file: {dest_file}")

    return


if __name__ == "__main__":
    # Gathering all necessary files
    # process_anat(directory)
    # process_func(directory)

    #step 1; change the smoothing kernel in the fsf files
    # find_and_update_fsf_files(directory)

    # # step 2 run feat
    # update_outputdir_and_run_feat(directory)

    # make_participant_tsv(HCP_dir, functional_fusion_dir)
    # import_anat_data(HCP_dir, functional_fusion_dir)
    # import_freesurfer(HCP_dir, functional_fusion_dir)

    # import_resms(HCP_dir, functional_fusion_dir)
    # import_masks(HCP_dir, functional_fusion_dir)
    import_betas(HCP_dir, functional_fusion_dir)
    # make_reginfo(HCP_dir, functional_fusion_dir)


    data_dir = f'/data/tge/Tian/BANDA'
    atlas_dir = '/data/tge/Tian/UKBB_full/imaging/Atlases'


    types = ['CondAll']
    atlases  = ['fs32k']
    session_list = ['ses-task']


    dataset = DataSetHcpTask(data_dir)
    for ses in session_list:
        print(f'extracting session {ses}')
        participants_tsv = pd.read_csv(f'{data_dir}/participants.tsv',sep = '\t')
        subj_list = participants_tsv['participant_id'].tolist()

        for type in types:
            print(f'extracting type {type}')
            for atlas in atlases:
                print(f'extracting atlas: {atlas}')
                dataset.extract_all(ses_id = ses,type = type, atlas = atlas, smooth=None, subj='all')
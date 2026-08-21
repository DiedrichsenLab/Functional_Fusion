import pandas as pd
from pathlib import Path
from Functional_Fusion.import_data import *
import Functional_Fusion.dataset as ds
import os
import nibabel as nb
import numpy as np
import gzip

base_dir = 'Y:/data'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data'
if not Path(base_dir).exists():
    base_dir = '/Volumes/diedrichsen_data$/data'

ERIS_DIR = '/home/dzhi/eris_mount'
if not Path(ERIS_DIR).exists():
    ERIS_DIR = '/data/tge'
if not Path(ERIS_DIR).exists():
    raise (NameError('Could not find hcp_dir'))

functional_fusion_dir = ERIS_DIR + f'/Tian/HCP_img'
HCP_dir = '/mnt/sda/HCP_tfMRI'

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


def import_anat_data(source_dir, dest_dir, subj_list):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
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


def import_freesurfer(source_dir, dest_dir, subj_list):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:

        pial_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant}.L.pial.32k_fs_LR.surf.gii'
        pial_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant}.R.pial.32k_fs_LR.surf.gii'
        white_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant}.L.white.32k_fs_LR.surf.gii'
        white_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant}.R.white.32k_fs_LR.surf.gii'
        sulc_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant}.L.sulc.32k_fs_LR.shape.gii'
        sulc_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant}.R.sulc.32k_fs_LR.shape.gii'

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


def import_resms(source_dir,dest_dir, subj_list):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
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

def import_betas(source_dir, dest_dir, subj_list):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        # Gather all session directories
        participant_dir = Path(source_dir) / str(participant) / 'func'
        session_dirs = sorted([d for d in participant_dir.iterdir()
                               if d.is_dir() and d.name.startswith("ses-")],
                               key=lambda d: d.name)

        session_run_mapping = {}

        for session_index, session_dir in enumerate(session_dirs):
            run_dirs = sorted([d for d in session_dir.iterdir()
                               if d.is_dir() and (d.name.endswith('LR') or d.name.endswith('RL'))],
                               key=lambda d: d.name)
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
                    stats_dir = feat_dir[0] / "stats"

                    if stats_dir.exists():
                        pe_files = sorted(list(stats_dir.glob("pe*.nii.gz")), key=lambda d: d.name)

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
                            dest_folder = Path(dest_dir) / "derivatives" / str(participant) / "estimates" / "ses-task"
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

def import_masks(source_dir, dest_dir, subj_list):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
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

def download_2lvl_glm_from_s3_server(subject_id, directory):
    # AWS S3 Bucket URL
    s3_base_url = "s3://hcp-openaccess/HCP_1200"

    """Download specific data folder for each subject."""
    # Construct the full S3 path for the subject and folder
    s3_mni_result_path = f"{s3_base_url}/{subject_id}/MNINonLinear/Results"
    local_folder = f"{directory}/{subject_id}"
    commands = []

    ## functional
    func_folder = os.path.join(local_folder, 'func')
    os.makedirs(func_folder, exist_ok=True)
    session_names = ['EMOTION', 'GAMBLING', 'LANGUAGE', 'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM']

    # Iterate over session directories that need preprocessing
    for session_name in session_names:
        ses_dir = os.path.join(func_folder, f'ses-{session_name}')
        os.makedirs(ses_dir, exist_ok=True)

        source_folder = os.path.join(s3_mni_result_path, f'tfMRI_{session_name}')
        commands.append(["aws", "s3", "sync", source_folder,
                         os.path.join(ses_dir, f'tfMRI_{session_name}'),
                         "--region", "us-east-1", "--exclude", "*", "--include", "*.dscalar.nii"])

    # Run the command to download the data
    for cmd in commands:
        try:
            subprocess.run(cmd, check=True)
            print(f"Successfully downloaded {cmd}")
        except subprocess.CalledProcessError as e:
            print(f"Error {cmd}: {e}")

def import_2lvl_task_contrast(source_dir, dest_dir, subj_list):
    session_names = ['EMOTION', 'GAMBLING', 'LANGUAGE', 'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM']
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        # Iterate over session directories that need preprocessing
        for session_name in session_names:
            print(f'Copying {participant}/ses-{session_name}')
            ses_dir = f'{source_dir}/{participant}/func/ses-{session_name}/tfMRI_{session_name}'
            contrast_dir = list(Path(ses_dir).rglob('*.dscalar.nii'))

            if len(contrast_dir) != 0:
                dest_folder = f'{dest_dir}/derivatives/{participant}/func/ses-{session_name}'
                if not Path(dest_folder).exists():
                    os.makedirs(dest_folder, exist_ok=True)

                # Copy files
                for contrast_file in contrast_dir:
                    d_file = shutil.copy(contrast_file, dest_folder)
                    os.chmod(d_file, 0o755)
            else:
                print(f'No contrasts for subject {participant}, skipping')

        # Average the mask data
        print(f"Copied all task contrasts for participant {participant}!")


def import_contrast(source_dir, dest_dir, subj_list, import_type="zstat"):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        # Gather all session directories
        participant_dir = Path(source_dir) / str(participant) / 'func'
        session_dirs = sorted([d for d in participant_dir.iterdir()
                               if d.is_dir() and d.name.startswith("ses-")],
                               key=lambda d: d.name)

        session_run_mapping = {}

        for session_index, session_dir in enumerate(session_dirs):
            run_dirs = sorted([d for d in session_dir.iterdir()
                               if d.is_dir() and (d.name.endswith('LR') or d.name.endswith('RL'))],
                               key=lambda d: d.name)
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
                    stats_dir = feat_dir[0] / "stats"

                    if stats_dir.exists():
                        pe_files = sorted(list(stats_dir.glob(f"{import_type}*.nii.gz")), key=lambda d: d.name)
                        for pe_file in pe_files:
                            reg_label = f"reg-{reg_num:02d}"

                            # Construct output filename
                            dest_folder = Path(dest_dir) / "derivatives" / str(participant) / "estimates" / "ses-task"
                            dest_file = dest_folder / f"{participant}_ses-task_run-{global_run_counter:02d}_{reg_label}_{import_type}.nii.gz"
                            if not dest_folder.exists():
                                os.makedirs(dest_folder, exist_ok=True)

                            # Copy contrast file
                            shutil.copyfile(pe_file, dest_file)
                            decompressed_file = dest_file.with_suffix('')
                            with gzip.open(dest_file, 'rb') as f_in:
                                with open(decompressed_file, 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)

                            # Remove the original .gz file
                            os.remove(dest_file)

                            reg_num += 1
                global_run_counter += 1
        print(f"Copied {import_type} files for participant {participant}")

    return

def make_task_contrasts(dataset_dir, subj_list="/subj_list/HCP203_test_set.tsv", smooth='2_MSMAll'):
    session_names = ['EMOTION', 'GAMBLING', 'LANGUAGE', 'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM']
    hcp_ds = ds.DataSetHcpTask(dataset_dir)
    T = hcp_ds.get_participants(subj_list)
    positive_ind = [1, 1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    beta_ind = [1, 1, 0, 0, 0, 0,
                    1, 1, 0, 0, 0, 0,
                    1, 1, 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 0, 0, 0, 0,
                    1, 1, 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    for s in T.participant_id:
        out_file = hcp_ds.func_dir.format(s) + \
                    f'/{s}_tfMRI_contrasts_level2_hp200_s{smooth}.dscalar.nii'
        if os.path.exists(out_file):
            print(f"Subject {s} combined contrasts already exists, skipping")
            continue
        else:
            ses_data, ses_info, ses_domains = [], [], []
            for sess in session_names:
                print(f'Loading {s}/ses-{sess}')
                ses_dir = hcp_ds.func_dir.format(s) + f'/ses-{sess}'
                contrast_file = ses_dir + f'/{s}_tfMRI_{sess}_level2_hp200_s{smooth}.dscalar.nii'

                # Load data / info
                dat = nb.load(contrast_file)
                this_info = dat.header.get_axis(0).name.tolist()
                prefix = os.path.commonprefix(this_info)
                this_info = [s[len(prefix):] for s in this_info]
                dat = dat.get_fdata().astype(np.float32)

                ses_data.append(dat)
                ses_info.append(this_info)
                ses_domains.append([sess] * len(this_info))

            data = np.vstack(ses_data)
            info = np.concatenate(ses_info)
            domains = np.concatenate(ses_domains)

            # Build new header
            C = nb.load(contrast_file)
            new_axis = nb.cifti2.ScalarAxis(info)
            bm = C.header.get_axis(1)  # brain models axis
            new_header = nb.cifti2.Cifti2Header.from_axes((new_axis, bm))

            # Save combined task contrasts
            C = nb.Cifti2Image(dataobj=data, header=new_header)
            nb.save(C, hcp_ds.func_dir.format(s) +
                    f'/{s}_tfMRI_contrasts_level2_hp200_s{smooth}.dscalar.nii')

            info_com = pd.DataFrame({'contrast_name': info,
                                     'task_name': domains,
                                     'positive': positive_ind,
                                     'betas': beta_ind})
            info_com.to_csv(hcp_ds.func_dir.format(s) +
                            f'/{s}_tfMRI_contrasts_level2_hp200.tsv', sep='\t', index=False)

            # Average the mask data
            print(f"Combined all task contrasts for participant {s}!")


def make_reginfo(source_dir, dest_dir, subj_list):
    participants = pd.read_csv(dest_dir + subj_list, sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        reginfo_data = []  # To store reginfo entries
        participant_dir = Path(source_dir) / str(participant) / "func"

        # Get session directories
        session_dirs = sorted([d for d in participant_dir.iterdir()
                               if d.is_dir() and d.name.startswith("ses-")],
                              key=lambda d: d.name)

        session_run_mapping = {}
        for session_index, session_dir in enumerate(session_dirs):
            run_dirs = sorted([d for d in session_dir.iterdir()
                               if d.is_dir() and (d.name.endswith('LR') or d.name.endswith('RL'))],
                              key=lambda d: d.name)
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
        dest_folder = Path(dest_dir) / "derivatives" / str(participant) / "estimates" / "ses-task"
        dest_file = dest_folder / f"{participant}_ses-task_reginfo.tsv"
        if not dest_folder.exists():
            os.makedirs(dest_folder, exist_ok=True)

        reginfo_df = pd.DataFrame(reginfo_data, columns=[ "run", "task_name", "cond_name", "reg_id","reg_num", "half"])
        reginfo_df.to_csv(dest_file, sep="\t", index=False)
        print(f"Saved reginfo file: {dest_file}")

    return

def make_zstat_reginfo(dataset_dir, subj_list):
    session_names = ['EMOTION', 'GAMBLING', 'LANGUAGE', 'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM']
    hcp_ds = ds.DataSetHcpTask(dataset_dir)
    T = hcp_ds.get_participants(subj_list)
    positive_ind = [1, 1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0] * 2
    beta_ind = [1, 1, 0, 0, 0, 0,
                1, 1, 0, 0, 0, 0,
                1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 0, 0, 0, 0,
                1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] * 2

    for s in T.participant_id:
        out_file = hcp_ds.data_dir.format(s) + f'/{s}_ses-task_ZstatHalf.tsv'
        if os.path.exists(out_file):
            df = pd.read_csv(out_file, sep="\t")
            df["positive"] = positive_ind
            df["beta"] = beta_ind
            df.to_csv(out_file, sep="\t", index=False)
            print(f"Saved reginfo file: {out_file}")
        else:
            print(f"Subject {s} zstat tsv not found")

if __name__ == '__main__':
    subj_list = "/subj_list/HCP200_test.tsv"
    # make_participant_tsv(HCP_dir, functional_fusion_dir)
    # import_anat_data(HCP_dir, functional_fusion_dir, subj_list)
    # import_freesurfer(HCP_dir, functional_fusion_dir, subj_list)
    # import_resms(HCP_dir, functional_fusion_dir, subj_list)
    # import_masks(HCP_dir, functional_fusion_dir, subj_list)
    #
    # import_2lvl_task_contrast('/mnt/sda/HCP_tfMRI', '/home/dzhi/eris_mount/Tian/HCP_img', subj_list)

    # make_task_contrasts('/home/dzhi/eris_mount/Tian/HCP_img', smooth='2_MSMAll', subj_list=subj_list)
    # import_contrast(HCP_dir, functional_fusion_dir, subj_list, import_type='zstat')
    # import_betas(HCP_dir, functional_fusion_dir, subj_list)
    # make_reginfo(HCP_dir, functional_fusion_dir, subj_list)
    make_zstat_reginfo(functional_fusion_dir, subj_list)



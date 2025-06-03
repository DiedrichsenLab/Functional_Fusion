import pandas as pd
from pathlib import Path
import Functional_Fusion.import_data as fi
import os
import nibabel as nb
import numpy as np
import gzip
import shutil

base_dir = 'Y:/data'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data'
if not Path(base_dir).exists():
    base_dir = '/Volumes/diedrichsen_data$/data'

functional_fusion_dir = f'{base_dir}/FunctionalFusion_new/HCP_tfMRI'
HCP_dir = f'{base_dir}/ExternalOpenData/HCP_UR100_tfMRI_full'


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
        anat_file = f'{source_dir}/{participant}/anat/T1W.nii'
        dest_folder = f'{dest_dir}/derivatives/ffimport/{participant}/anat/'
        dest_file = f'{dest_folder}/{participant}_T1w.nii'

        if not Path(dest_folder).exists():
            os.makedirs(dest_folder, exist_ok=True)
        
        # copy anat file
        shutil.copyfile(anat_file, dest_file)
        print(f'Copied {anat_file} to {dest_file} for participant {participant}')

    return


def import_suit_data(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"]

    for participant in participants:

        source_folder = f'{source_dir}/{participant}/suit'
        dest_folder = f'{dest_dir}/derivatives/ffimport/{participant}/anat/'
        dest_file = f'{dest_folder}/{participant}_T1w.nii'

        fi.import_suit(source_folder, dest_folder, 'T1w', participant)
        print(f'Copied participant {participant}')

    return


def import_freesurfer(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        participant_number = participant.split("-")[1]
        pial_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.L.pial.32k_fs_LR.surf.gii'
        pial_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.R.pial.32k_fs_LR.surf.gii'
        white_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.L.white.32k_fs_LR.surf.gii'
        white_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.R.white.32k_fs_LR.surf.gii'
        sulc_L_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.L.sulc.32k_fs_LR.shape.gii'
        sulc_R_source = f'{source_dir}/{participant}/SurfaceWB/{participant_number}.R.sulc.32k_fs_LR.shape.gii'

        pial_L_dest = f'{dest_dir}/derivatives/ffimport/{participant}/anat/{participant}_space-32k_hemi-L_pial.surf.gii'
        pial_R_dest = f'{dest_dir}/derivatives/ffimport/{participant}/anat/{participant}_space-32k_hemi-R_pial.surf.gii'
        white_L_dest = f'{dest_dir}/derivatives/ffimport/{participant}/anat/{participant}_space-32k_hemi-L_white.surf.gii'
        white_R_dest = f'{dest_dir}/derivatives/ffimport/{participant}/anat/{participant}_space-32k_hemi-R_white.surf.gii'
        sulc_L_dest = f'{dest_dir}/derivatives/ffimport/{participant}/anat/{participant}_space-32k_hemi-L_sulc.shape.gii'
        sulc_R_dest = f'{dest_dir}/derivatives/ffimport/{participant}/anat/{participant}_space-32k_hemi-R_sulc.shape.gii'

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
        dest_folder = f'{dest_dir}/derivatives/ffimport/{participant}/func/ses-task'
        dest_file = f'{dest_folder}/{participant}_ses-task_resms.nii.'
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
        print(f"Processing participant: {participant}")
        # Gather all session directories
        participant_dir = Path(source_dir) / participant / 'func'
        session_dirs = [d for d in participant_dir.iterdir() if d.is_dir() and d.name.startswith("ses-")]

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
                    stats_dir = feat_dir[0] / "stats"

                    if stats_dir.exists():
                        pe_files = list(stats_dir.glob("pe*.nii.gz"))
                        pe_files.sort(key=lambda f: int(f.name[2:].split(".")[0]))


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
                            dest_folder = Path(dest_dir) / "derivatives" /'ffimport'/participant/ "func" / "ses-task"
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

        dest_folder = f'{dest_dir}/derivatives/ffimport/{participant}/func/ses-task'
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
        dest_folder = Path(dest_dir) / "derivatives" / 'ffimport'/ participant / "func" / "ses-task"
        dest_file = dest_folder / f"{participant}_ses-task_reginfo.tsv"
        if not dest_folder.exists():
            os.makedirs(dest_folder, exist_ok=True)

        reginfo_df = pd.DataFrame(reginfo_data, columns=[ "run", "task_name", "cond_name", "reg_id","reg_num", "half"])
        reginfo_df.to_csv(dest_file, sep="\t", index=False)
        print(f"Saved reginfo file: {dest_file}")

    return



if __name__ == '__main__':
    # make_participant_tsv(HCP_dir, functional_fusion_dir)
    # import_anat_data(HCP_dir, functional_fusion_dir)
    # import_freesurfer(HCP_dir, functional_fusion_dir)
    # import_resms(HCP_dir, functional_fusion_dir)
    # import_masks(HCP_dir, functional_fusion_dir)
    import_betas(HCP_dir, functional_fusion_dir)
    # make_reginfo(HCP_dir, functional_fusion_dir)
    # import_suit_data(HCP_dir, functional_fusion_dir)


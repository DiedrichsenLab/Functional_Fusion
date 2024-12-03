import pandas as pd
from pathlib import Path
from Functional_Fusion.import_data import *
import os
import nibabel as nb
import numpy as np
import gzip

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data'

functional_fusion_dir = f'{base_dir}/FunctionalFusion/HCP_tfMRI'
HCP_dir = f'{base_dir}/ExternalOpenData/HCP_UR100_tfMRI_new'


def make_participant_tsv(source_dir, dest_dir):
    if not Path(dest_dir).exists():
        os.makedirs(dest_dir)
    
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
        dest_folder = f'{dest_dir}/derivatives/{participant}/anat/'
        dest_file = f'{dest_folder}/{participant}_T1w.nii'

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
        participant_number = participant.split("-")[1]
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
        dest_folder = f'{dest_dir}/derivatives/{participant}/estimates/ses-all'
        dest_file = f'{dest_folder}/{participant}_ses-all_resms.nii.'
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
        run_counter = 1

        # Get session dirs
        participant_dir = Path(source_dir) / participant / 'func'
        session_dirs = []

        for d in participant_dir.iterdir():
            if d.is_dir():
                if d.name.startswith("ses-"):
                    session_dirs.append(d)

        for session_dir in session_dirs:
            # Get run dirs
            run_dirs = [d for d in session_dir.iterdir() if d.is_dir()]

            for run_dir in run_dirs:
                # Path to the stats folder
                feat_dir = list(run_dir.glob("*.feat"))
                stats_dir = feat_dir[0] / "stats"

                if stats_dir.exists():
                    pe_files = list(stats_dir.glob("pe*.nii.gz"))

                    # Filter for odd-numbered PE files (derivatives)
                    odd_pe_files = []
                    for pe in pe_files:
                        beta_number = int(pe.name.replace(".nii.gz", "").replace("pe", ""))
                        if beta_number % 2 == 1:
                            odd_pe_files.append(pe)

                    # Enumerate odd PE files for 0-indexed numbering
                    for new_index, pe_file in enumerate(odd_pe_files):
                        beta_number = int(pe_file.name.replace(".nii.gz", "").replace("pe", ""))
                        reg_label = f"reg-{new_index:02d}"

                        # Construct output filename
                        dest_folder = Path(dest_dir) / "derivatives" / participant / "estimates" / "ses-all"
                        dest_file = dest_folder / f"{participant}_ses-all_run-{run_counter:02d}_{reg_label}_beta.nii.gz"

                        if not dest_folder.exists():
                            os.makedirs(dest_folder, exist_ok=True)

                        # copy beta file
                        shutil.copyfile(pe_file, dest_file)
                        decompressed_file = dest_file.with_suffix('') 
                        with gzip.open(dest_file, 'rb') as f_in:
                            with open(decompressed_file, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)

                        # Remove the original .gz file
                        os.remove(dest_file)

                # Add to run counter
                run_counter += 1
        print(f"Copied beta files for participant {participant}")

    return
def import_masks(source_dir, dest_dir):
    participants = pd.read_csv(Path(dest_dir) / "participants.tsv", sep="\t")
    participants = participants["participant_id"].tolist()

    for participant in participants:
        mask_files = list(Path(f'{source_dir}/{participant}').rglob('mask.nii.gz'))

        dest_folder = f'{dest_dir}/derivatives/{participant}/estimates/ses-all'
        dest_file = f'{dest_folder}/{participant}_ses-all_mask.nii'
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
        reginfo_data = []
        participant_dir = Path(source_dir) / participant / "func"
        
        session_dirs = [d for d in participant_dir.iterdir() if d.is_dir() and d.name.startswith("ses-")]

        run_counter = 1
        reg_num = 0

        for session_dir in session_dirs:
            # get the session name after the -
            session_name = session_dir.name.split("-")[1]
            run_dirs = [d for d in session_dir.iterdir() if d.is_dir()]
            
            for run_dir in run_dirs:
                # Locate .feat folder
                feat_dir = list(run_dir.glob("*.feat"))
                design_fsf = feat_dir[0] / "design.fsf"
                
                # Parse design.fsf to extract condition names
                with open(design_fsf, "r") as fsf_file:
                    lines = fsf_file.readlines()
                
                conditions = []
                for line in lines:
                    if line.startswith("set fmri(evtitle"):
                        cond_name = line.split('"')[1] 
                        conditions.append(cond_name)
                
                for cond_name in conditions:
                    reginfo_data.append({
                        "cond_name": cond_name,
                        "run": run_counter,
                        "reg_num": reg_num,
                        "run_type": session_name
                    })
                    reg_num += 1
                
                run_counter += 1

        # Save reginfo.tsv
        dest_folder = Path(dest_dir) / "derivatives" / participant / "estimates" / "ses-all"
        dest_file = dest_folder / f"{participant}_ses-all_reginfo.tsv"
        if not dest_folder.exists():
            os.makedirs(dest_folder, exist_ok=True)
        
        reginfo_df = pd.DataFrame(reginfo_data, columns=["cond_name", "run", "reg_num", "run_type"])
        reginfo_df.to_csv(dest_file, sep="\t", index=False)
        print(f"Saved reginfo file: {dest_file}")

    return


if __name__ == '__main__':
    # make_participant_tsv(HCP_dir, functional_fusion_dir)
    # import_anat_data(HCP_dir, functional_fusion_dir)
    # import_freesurfer(HCP_dir, functional_fusion_dir)
    # import_resms(HCP_dir, functional_fusion_dir)
    # import_masks(HCP_dir, functional_fusion_dir)
    # import_betas(HCP_dir, functional_fusion_dir)
    make_reginfo(HCP_dir, functional_fusion_dir)



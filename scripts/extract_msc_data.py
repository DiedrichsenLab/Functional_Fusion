# Script for importing the MSC data set to general format.
import pandas as pd
import shutil, re, gzip
from pathlib import Path
import mat73, subprocess
import numpy as np
import sys, os, time
import Functional_Fusion.atlas_map as am
import Functional_Fusion.util as ut
import Functional_Fusion.connectivity as conn
from Functional_Fusion.dataset import DataSetMSC
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/home/dzhi/eris_mount/Tian/UKBB_full/imaging'
if not Path(base_dir).exists():
    base_dir = '/data/tge/Tian/UKBB_full/imaging'

data_dir = '/data/tge/Tian/MSC'
atlas_dir = base_dir + '/Atlases'

def extract_MSC(ses_id='ses-01',type='CondHalf',atlas='MNISymC3'):
    dataset = DataSetMSC(data_dir)
    dataset.extract_all(ses_id, type, atlas, smooth=None)


def smooth_MSC_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1, kernel='gaussian'):
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Smoothing data for {s} fs32k {ses_id} in {smooth}mm {kernel} ...')

        start = time.perf_counter()
        file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        ut.smooth_fs32k_data(file, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


def mask_MSC_fs32k(ses_id='ses-s1', type='CondHalf', high_percent=0.1, low_percent=0.1,
                           smooth=None, z_transfer=False, binarized=False):
    myatlas, _ = am.get_atlas('fs32k')
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Mask data for {s} fs32k {ses_id} in high {high_percent} low {low_percent} ...')

        start = time.perf_counter()
        if smooth is not None:
            file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        else:
            file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'

        ut.mask_fs32k_data(file, high_percent=high_percent, low_percent=low_percent,
                           z_transfer=z_transfer, binarized=binarized)

        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


def import_betas(source_dir, dest_dir):
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()

    for i, s in enumerate(T.participant_id):
        orig_subid = T.participant_num[i]
        # Gather all session directories
        betas_dir = Path(source_dir) / 'GLM_level1'
        session_info = ['motor','mem','mixed']

        dest_folder = Path(dataset.estimates_dir.format(s) + "/ses-task")
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder, exist_ok=True)

        reginfo_data = []  # To store reginfo entries
        global_run_counter = 1  # Global counter for runs across sessions
        for run in range(1,11):
            reg_id = 1
            for session_index, session_name in enumerate(session_info):

                beta_files = list(betas_dir.glob(f"{orig_subid}_{session_name}_session{run}_*beta.nii.gz"))

                for reg_num, file in enumerate(beta_files):
                    # 1. reorient the file to MNI standard
                    reoriented_file = file.parent / ('reoriented_' + file.name)
                    cmd1 = f'fslswapdim {file} -x y z {reoriented_file}'
                    cmd2 = f'fslorient -swaporient {reoriented_file}'

                    # 2. upsampling from 3mm to 2mm
                    upsampled_file = file.parent / ('upsampled_' + file.name)
                    cmd3 = f'flirt -in {reoriented_file} -ref ' \
                        f'$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz ' \
                        f'-out {upsampled_file} -applyisoxfm 2 -interp trilinear -usesqform'

                    pattern = f"{re.escape(f'{session_name}_session{run}_')}(.*?){re.escape('_beta.nii.gz')}"
                    cond_name = re.search(pattern, file.name).group(1)

                    # Run commands sequentially
                    try:
                        for cmd in [cmd1, cmd2, cmd3]:
                            process = subprocess.run(cmd, shell=True, check=True)

                        # make reg info file
                        reginfo_data.append({"sn": run,
                                            "run": global_run_counter,
                                            "task_name": session_name,
                                            "task_uni_num": session_index + 1,
                                            "cond_name": cond_name,
                                            "reg_id": reg_id,
                                            "reg_num": reg_num + 1,
                                            "half": 2 if run % 2 == 0 else 1})

                        # Write-in beta files in destination folder
                        dest_file = dest_folder / f"{s}_ses-task_run-{global_run_counter:02d}_reg-{reg_id:02d}_beta.nii.gz"
                        # Copy beta file and decompress as as .nii
                        shutil.copyfile(upsampled_file, dest_file)
                        decompressed_file = dest_file.with_suffix('')
                        with gzip.open(dest_file, 'rb') as f_in:
                            with open(decompressed_file, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)

                        reg_id += 1
                        # remove intermediate files
                        os.remove(reoriented_file)
                        os.remove(upsampled_file)
                        os.remove(dest_file)
                        print(f'Successfully import betas for subject {s}: {file}!')

                    except:
                        print(f'Failed import betas for subject {s}: {file}')

                global_run_counter += 1

        # Write-in reg info file
        reginfo_df = pd.DataFrame(reginfo_data, columns=["sn", "run", "task_name", "cond_name", "reg_id","reg_num", "half"])
        reginfo_df.to_csv(dest_folder / f"{s}_ses-task_reginfo.tsv", sep="\t", index=False)
        print(f"Saved task reginfo file for subject {s}")


def import_freesurfer(source_dir, dest_dir):
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()
    fs32k_dir = source_dir + '/derivatives/surface_pipeline/{0}/fs_LR_Talairach/fsaverage_LR32k'

    for i, s in enumerate(T.participant_id):
        orig_subid = T.participant_num[i]

        if not os.path.exists(dataset.anatomical_dir.format(s)):
            os.makedirs(dataset.anatomical_dir.format(s))

        pial_L_source = fs32k_dir.format(s) + f'/{orig_subid}.L.pial.32k_fs_LR.surf.gii'
        pial_R_source = fs32k_dir.format(s) + f'/{orig_subid}.R.pial.32k_fs_LR.surf.gii'
        white_L_source = fs32k_dir.format(s) + f'/{orig_subid}.L.white.32k_fs_LR.surf.gii'
        white_R_source = fs32k_dir.format(s) + f'/{orig_subid}.R.white.32k_fs_LR.surf.gii'
        sulc_L_source = fs32k_dir.format(s) + f'/{orig_subid}.L.sulc.32k_fs_LR.shape.gii'
        sulc_R_source = fs32k_dir.format(s) + f'/{orig_subid}.R.sulc.32k_fs_LR.shape.gii'

        pial_L_dest = dataset.anatomical_dir.format(s) + f'/{s}_space-32k_hemi-L_pial.surf.gii'
        pial_R_dest = dataset.anatomical_dir.format(s) + f'/{s}_space-32k_hemi-R_pial.surf.gii'
        white_L_dest = dataset.anatomical_dir.format(s) + f'/{s}_space-32k_hemi-L_white.surf.gii'
        white_R_dest = dataset.anatomical_dir.format(s) + f'/{s}_space-32k_hemi-R_white.surf.gii'
        sulc_L_dest = dataset.anatomical_dir.format(s) + f'/{s}_space-32k_hemi-L_sulc.shape.gii'
        sulc_R_dest = dataset.anatomical_dir.format(s) + f'/{s}_space-32k_hemi-R_sulc.shape.gii'

        # copy all files
        shutil.copyfile(pial_L_source, pial_L_dest)
        shutil.copyfile(pial_R_source, pial_R_dest)
        shutil.copyfile(white_L_source, white_L_dest)
        shutil.copyfile(white_R_source, white_R_dest)
        shutil.copyfile(sulc_L_source, sulc_L_dest)
        shutil.copyfile(sulc_R_source, sulc_R_dest)

        print(f'Copied freesurfer files for participant {s}')

    return


def import_masks():
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()

    for participant in T.participant_id:
        task_ses_dir = dataset.estimates_dir.format(participant) + '/ses-task'
        beta_files = list(Path(task_ses_dir).rglob('*beta.nii'))
        unique_runs = sorted(set(f.name.split('_')[2] for f in beta_files))

        masks = []
        for run in unique_runs:
            mask_data = []
            run_beta_files = list(Path(task_ses_dir).rglob(f'*{run}*beta.nii'))
            for f in run_beta_files:
                img = nb.load(f)
                data = img.get_fdata()
                mask = (data != 0).astype(np.uint8)
                mask_data.append(mask)

            # Average the mask data
            mask_data = np.mean(mask_data, axis=0)
            if len(np.unique(mask_data)) != 2:
                print(f'Mismatch in {participant} run {run}!')

            masks.append(mask_data)

        masks = np.mean(masks, axis=0)
        mask_img = nb.Nifti1Image(masks, img.affine)
        nb.save(mask_img, f'{task_ses_dir}/{participant}_ses-task_mask.nii')
        print(f"Mask file successfully generated for participant {participant}")

    return

def import_resms():
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()

    for participant in T.participant_id:
        task_ses_dir = dataset.estimates_dir.format(participant) + '/ses-task'
        mask_files = task_ses_dir + f'/{participant}_ses-task_mask.nii'
        img = nb.load(mask_files)
        mask = img.get_fdata()

        resms = np.where(mask!=0 , 1, 0).astype(np.float64)
        resms_img = nb.Nifti1Image(resms, img.affine)
        nb.save(resms_img, f'{task_ses_dir}/{participant}_ses-task_resms.nii')
        print(f"Residual mean squared file successfully generated for participant {participant}")


def import_preprocessed_task_contrasts(source_dir, target_dir):
    dataset = DataSetMSC(data_dir)
    T = dataset.get_participants()

    for i, s in enumerate(T.participant_num):
        contrasts_folder = source_dir + f'/derivatives/surface_pipeline/sub-{s}/task_contrasts_cifti'
        dest_dir = dataset.base_dir + f'/derivatives/sub-{s}/task_contrasts'
        os.makedirs(dest_dir, exist_ok=True)

        # Copy each folder from the source to the destination
        for folder_name in os.listdir(contrasts_folder):
            folder_path = os.path.join(contrasts_folder, folder_name)
            if os.path.isdir(folder_path):
                shutil.copytree(folder_path, os.path.join(dest_dir, folder_name), dirs_exist_ok=True)
                print(f'Copied {os.path.join(dest_dir, folder_name)}')


def extract_msc_rs_timeseries(data_dir, derivative_dir, ses_id=1):
    myatlas, _ = am.get_atlas('fs32k')
    myatlas.calculate_symmetry()

    T = pd.read_csv(data_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        dest_dir = f'{derivative_dir}/derivatives/{s}/data'
        # if os.path.exists(dest_dir + f'/{s}_space-fs32k_{ses_id}_Tseries.dscalar.nii'):
        #     print(f"Already imported subject {s} {ses_id} Tseries, skipping...")
        # else:
        source_dir = f'{data_dir}/derivatives/surface_pipeline/{s}' \
                     f'/processed_restingstate_timecourses/ses-func{ses_id:02}/cifti'

        print(f"-- Start importing subject {s} {ses_id} Tseries --")
        try:
            info = pd.DataFrame()
            start = time.perf_counter()
            data_file = source_dir + f'/{s}_ses-func{ses_id:02}_task-rest_bold_32k_fsLR.dtseries.nii'
            mask_file = source_dir + f'/{s}_ses-func{ses_id:02}_task-rest_bold_32k_fsLR_tmask.txt'
            tmp_file = source_dir + f'/{s}_ses-func{ses_id:02}_tmp_cifti2.dtseries.nii'

            # Write in smoothed surface data (filled with 0)
            cmd = f"wb_command -file-convert -cifti-version-convert " \
                        f"{data_file} 2 {tmp_file}"
            subprocess.run(cmd, shell=True)

            data = myatlas.cifti_to_data(tmp_file)
            t_mask = np.loadtxt(mask_file)
            data[~t_mask.astype(bool),:] = np.nan
            os.remove(tmp_file)

            num_tpoint = [f'T{i+1:04}' for i in range(data.shape[0])]
            info = pd.DataFrame({'sn': [s] * data.shape[0],
                                      'run': [ses_id] * data.shape[0],
                                      'timepoint': num_tpoint,
                                      'task': 'rest',
                                      'time_id': [i+1 for i in range(data.shape[0])],
                                      'names': num_tpoint,
                                      'mask': t_mask.astype(int)})

            C = myatlas.data_to_cifti(data, num_tpoint)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)

            nb.save(C, dest_dir +
                        f'/{s}_space-fs32k_ses-rest{ses_id}_Tseries.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_ses-rest{ses_id}_Tseries.tsv',
                        sep='\t', index=False)

            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - time {elapse}")
        except:
            print(f"   Failed - {s} ses-rest{ses_id} Tseries!")


if __name__ == "__main__":
    # import_preprocessed_task_contrasts('/data/tge/dzhi/projects/ds000224-download', data_dir)
    # import_freesurfer('/data/tge/dzhi/projects/ds000224-download', data_dir)
    # import_betas('/data/tge/dzhi/projects/ds000224-download', data_dir)
    # import_masks()
    # import_resms()

    # Resting-state preprocessing
    for i in range(1,11):
        # extract_msc_rs_timeseries('/data/tge/dzhi/projects/ds000224-download',
        #                     data_dir, ses_id=i)
        conn.get_connectivity_fingerprint('MSC',
                                    type='Ico642Run', space='fs32k', ses_id=f'ses-rest{i}')
    # smooth_MSC_fs32k(ses_id='ses-motor', type='CondHalf', smooth=4, kernel='fwhm')

    # --- Extracting Estimates ---
    # extract_MSC(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    extract_MSC(ses_id='ses-task', type='CondRun', atlas='fs32k')
    # extract_MSC(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # extract_MSC(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')
    for s in [4,6,8,10]:
        print(f'Doing processing for {s}fwhm ...')
        mask_MSC_fs32k(ses_id='ses-motor', type=f'CondHalf', high_percent=0.1,
                         low_percent=0.1, smooth=f'{s}fwhm', z_transfer=True, binarized=False)

    # --- Group Average ---
    dataset = DataSetMSC(data_dir)
    dataset.extract_all(type='CondAll', ses_id='ses-motor', atlas='MNISymC3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='SUIT3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC3')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='fs32k')
    # dataset.group_average_data(ses_id='ses-motor', type='CondHalf', atlas='MNISymC2')


    # --- Show group average ---
    # dataset.plot_cerebellum(subject='group', savefig=True, colorbar=True)
    pass

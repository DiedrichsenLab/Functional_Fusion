import os
import subprocess
import shutil

new_directory = 'Y:/data/ExternalOpenData/HCP_UR100_new/tasktest'
if not os.path.exists(new_directory):
    new_directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_new/tasktest'

def update_all_fsf_files(new_directory):

    for root, dirs, files in os.walk(new_directory):
        if root.endswith(".feat"):
            continue

        for file in files:
            # Only consider .fsf files
            if not file.endswith(".fsf"):
                continue

            # Only update true HCP Level1 fsf designs
            if not file.startswith("tfMRI_"):
                # debug: print(f"[SKIP] Not tfMRI_* fsf: {file}")
                continue

            fsf_path = os.path.join(root, file)
            fsf_dir = os.path.dirname(fsf_path)

            print(f"\nUpdating FSF file: {fsf_path}")

            # Parse task + run from filename
            parts = file.split("_")
            task = parts[1]
            run  = parts[2]

            # ----- Build output directory path -----
            base = os.path.splitext(file)[0]
            outdir = os.path.join(fsf_dir, f"{base}.feat")

            # normalize + force Linux cluster path
            outdir = outdir.replace("\\", "/").replace("Y:/data", "/cifs/diedrichsen/data")

            # ----- Build cleaned nifti path -----
            new_nii_name = f"tfMRI_{task}_{run}_hp0_clean_rclean_tclean.nii.gz"
            new_nii_path = os.path.join(fsf_dir, new_nii_name)
            new_nii_path = new_nii_path.replace("\\", "/").replace("Y:/data", "/cifs/diedrichsen/data")

            # ----- Rewrite FSF -----
            updated = []
            with open(fsf_path, "r") as f:
                for line in f:
                    if line.startswith("set fmri(outputdir)"):
                        updated.append(f'set fmri(outputdir) "{outdir}"\n')
                    elif line.startswith("set feat_files(1)"):
                        updated.append(f'set feat_files(1) "{new_nii_path}"\n')
                    else:
                        updated.append(line)

            with open(fsf_path, "w") as f:
                f.writelines(updated)


def run_all_glms(new_directory):
    for root, dirs, files in os.walk(new_directory):
        if root.endswith(".feat"):
            continue

        for file in files:
            # must be .fsf
            if not file.endswith(".fsf"):
                continue

            # only HCP level1 designs
            if not file.startswith("tfMRI_"):
                continue

            fsf_path = os.path.join(root, file)
            fsf_dir = os.path.dirname(fsf_path)

            print(f"\nFound FSF: {fsf_path}")

            # Determine expected .feat output dir
            base = os.path.splitext(file)[0]
            feat_output_dir = os.path.join(fsf_dir, f"{base}.feat")

            # Run FEAT
            try:
                subprocess.run(["feat", fsf_path], check=True)
                print(f"[DONE] Completed: {fsf_path}")

            except subprocess.CalledProcessError as e:
                print(f"[ERROR] FEAT FAILED for {fsf_path}\n{e}")

            except Exception as e:
                print(f"[UNEXPECTED ERROR] {e}")

def reorg(new_dir):

    for subj_folder in os.listdir(new_dir):
        subj_path = os.path.join(new_dir, subj_folder)
        if not os.path.isdir(subj_path):
            continue
        
        # Extract subject ID from something like 101309_Task3TRecommended
        subj_id = subj_folder.split("_")[0]

        # Find original Results folder
        results_dir = os.path.join(subj_path, subj_id, "MNINonLinear", "Results")
        if not os.path.exists(results_dir):
            continue
    
        # Create reorganized output inside the SAME root
        out_sub_dir = os.path.join(new_dir, f"sub-{subj_id}", "func")
        os.makedirs(out_sub_dir, exist_ok=True)

        # Loop through all tfMRI folders
        for run_folder in os.listdir(results_dir):

            if not run_folder.startswith("tfMRI_"):
                continue

            # skip concat folders
            if "concat" in run_folder.lower():
                continue

            run_path = os.path.join(results_dir, run_folder)
            if not os.path.isdir(run_path):
                continue

            # Parse tfMRI_<SESSION>_<RUN>
            parts = run_folder.split("_")
            session = parts[1]      # e.g. EMOTION

            # Create session folder
            ses_dir = os.path.join(out_sub_dir, f"ses-{session}")
            os.makedirs(ses_dir, exist_ok=True)

            # Destination folder
            dest_run_dir = os.path.join(ses_dir, run_folder)

            # Copy folder
            shutil.copytree(run_path, dest_run_dir)



if __name__ == "__main__":
    # update_all_fsf_files(new_directory)
    # run_all_glms(new_directory)
    reorg(new_directory)

    
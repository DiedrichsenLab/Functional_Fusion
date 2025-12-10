import os
import zipfile
import shutil
import nibabel as nb
import os
import shutil
import glob

new_directory = 'Y:/data/ExternalOpenData/HCP_UR100_new/tasktest'
old_dir     = "Y:/data/ExternalOpenData/HCP_UR100_tfMRI_full"
if not os.path.exists(new_directory):
    new_directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_new/tasktest'
if not os.path.exists(old_dir): 
    old_dir = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_full'

# step 1: Delete .md5 files and extract .zip files into individual folders
def unzip_clean(directory):
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)

        # Delete all .md5 files
        if item.endswith('.md5'):
            os.remove(item_path)

        # If it's a zip file, extract it into a folder with the same name and delete the zip
        elif item.endswith('.zip'):
            extract_dir = os.path.join(directory, os.path.splitext(item)[0])
            os.makedirs(extract_dir, exist_ok=True)
            
            with zipfile.ZipFile(item_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
            
            os.remove(item_path)


def copy_fsf_files(new_directory, old_dir):
    # Task names used by HCP
    task_list = [
        "EMOTION",
        "LANGUAGE",
        "GAMBLING",
        "MOTOR",
        "RELATIONAL",
        "SOCIAL",
        "WM"
    ]

    runs = ["LR", "RL"]

    # Loop over all subject folders in tasktest (e.g., '101309_Task3TRecommended')
    for subj_folder in os.listdir(new_directory):
        subj_id = subj_folder.split("_")[0]  
        print(f"\n=== Subject {subj_id} ===")

        # Paths
        results_dir = os.path.join(new_directory, subj_folder, subj_id, "MNINonLinear", "Results")
        full_subj_dir = os.path.join(old_dir, f"sub-{subj_id}", "func")

        # For each task × run
        for task in task_list:
            for run in runs:
                res_task_dir = os.path.join(results_dir, f"tfMRI_{task}_{run}")
                src_task_dir = os.path.join(full_subj_dir, f"ses-{task}", f"tfMRI_{task}_{run}")

                # Find fsf file(s)
                fsf_files = glob.glob(os.path.join(src_task_dir, "*.fsf"))
                if len(fsf_files) == 0:
                    print(f"   No .fsf file in {src_task_dir}")
                    continue

                for fsf in fsf_files:
                    dst = os.path.join(res_task_dir, os.path.basename(fsf))
                    print(f"   Copying {fsf} → {dst}")
                    shutil.copy2(fsf, dst)



unzip_clean(new_directory)
copy_fsf_files(new_directory, old_dir)

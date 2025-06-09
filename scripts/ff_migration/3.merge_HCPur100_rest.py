import os
import shutil

def copy_hcp_func_data(source_base, dest_base):
    for sub_folder in os.listdir(dest_base):
        if not sub_folder.startswith("sub-"):
            continue

        subject_num = sub_folder.replace("sub-", "")
        source_func_dir = os.path.join(source_base, subject_num, "func")
        dest_func_dir = os.path.join(dest_base, sub_folder, "func", "ses-rest")

        os.makedirs(dest_func_dir, exist_ok=True)

        # Copy each file
        for file in os.listdir(source_func_dir):
            src_file = os.path.join(source_func_dir, file)
            dest_file = os.path.join(dest_func_dir, file)

            shutil.copy2(src_file, dest_file)
            print(f"Copied: {src_file} → {dest_file}")


if __name__ == '__main__':
    source_hcp = r'Y:\data\FunctionalFusion\HCP\derivatives'
    dest_ffimport = r'Y:\data\FunctionalFusion_new\HCPur100\derivatives\ffimport'

    copy_hcp_func_data(source_hcp, dest_ffimport)

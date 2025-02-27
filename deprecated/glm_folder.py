import os
import shutil

def clean_output_directory(output_dir):
    """Cleans the output directory if it exists, and recreates it."""
    if os.path.exists(output_dir):
        print(f"Cleaning existing directory: {output_dir}")
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

def process_feat_folder(feat_folder, output_dir):
    """Processes a single .feat folder to copy specific files."""
    print(f"Processing feat folder: {feat_folder}")
    
    file_patterns = [
        ("mean_func.nii.gz", "Copying mean_func file..."),
        ("sigmasquared", "Copying sigmasquared files..."),
        ("pe", "Copying files starting with pe...")
    ]
    
    for pattern, message in file_patterns:
        print(message)
        for root, _, files in os.walk(feat_folder):
            for file_name in files:
                if pattern == "mean_func.nii.gz" and file_name == pattern:
                    copy_file(file_name, root, output_dir)
                elif pattern == "sigmasquared" and file_name.startswith(pattern):
                    copy_file(file_name, root, output_dir)
                elif pattern == "pe" and file_name.startswith(pattern):
                    copy_file(file_name, root, output_dir)

def copy_file(file_name, src_dir, dest_base_dir):
    """Copies a file from the source directory to the destination directory, preserving folder structure."""
    src_path = os.path.join(src_dir, file_name)
    relative_path = os.path.relpath(src_dir, start=feat_base_dir)
    dest_dir = os.path.join(dest_base_dir, relative_path)
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, file_name)
    shutil.copy(src_path, dest_path)
    print(f"Copied {file_name} to {dest_path}")

# Main script execution
if __name__ == "__main__":
    feat_base_dir = "/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_new"
    output_base_dir = os.path.join(feat_base_dir, "glm")
    clean_output_directory(output_base_dir)

    for root, dirs, _ in os.walk(feat_base_dir):
        for dir_name in dirs:
            if dir_name.endswith(".feat"):
                feat_folder = os.path.join(root, dir_name)
                process_feat_folder(feat_folder, output_base_dir)

    print("Processing complete!")
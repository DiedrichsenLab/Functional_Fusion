import os

# Base directory to start searching
base_dir = "/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_new"

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

if __name__ == "__main__":
    find_and_update_fsf_files(base_dir)
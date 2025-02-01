import os
import subprocess

def update_outputdir_and_run_feat(base_dir):
    """
    Updates the output directory in each .fsf file to match its location
    and runs FEAT on the updated file.
    """
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".fsf"):
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

# Base directory to search for .fsf files
base_directory = "/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_new"

# Run the function
update_outputdir_and_run_feat(base_directory)
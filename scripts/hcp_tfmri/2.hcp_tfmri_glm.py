import os
import subprocess




directory = 'Y:\\data\\ExternalOpenData\\HCP_UR100_tfMRI_new'
if not os.path.exists(directory):
    directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_new'


# step 1: change the smoothing kernel in the fsf files
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


if __name__ == "__main__":
    #step 1; change the smoothing kernel in the fsf files
    find_and_update_fsf_files(directory)

    # # step 2 run feat
    update_outputdir_and_run_feat(directory)



    

    
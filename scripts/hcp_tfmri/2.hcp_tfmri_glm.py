import os, time, sys
import subprocess
import pandas as pd
from multiprocessing import Pool, cpu_count


directory = '/data/tge/dzhi/projects/HCP_tfMRI'
if not os.path.exists(directory):
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
                if os.path.exists(output_dir):
                    to_remove = f"{os.path.splitext(file)[0]}*.feat"
                    subprocess.run(f"rm -r {os.path.join(root, to_remove)}", 
                                   shell=True, check=True)
                
                # Run FEAT
                try:
                    print(f"Running FEAT on: {fsf_path}")
                    subprocess.run(["feat", fsf_path], check=True)
                    print(f"FEAT successfully completed for: {fsf_path}")
                except subprocess.CalledProcessError as e:
                    print(f"Error running FEAT on {fsf_path}: {e}")
                except Exception as e:
                    print(f"Unexpected error: {e}")


def run_feat_all(set_index):
    T = pd.read_csv('/data/tge/Tian/HCP_img/subj_list/HCP203_test_set.tsv', delimiter='\t')
    subj_idx = [set_index, set_index+100]
    
    for i in subj_idx:
        s = T.participant_id[i]

        if not os.path.exists(f'{directory}/{s}/func/ses-WM/tfMRI_WM_RL/tfMRI_WM_RL_hp200_s4_level1.feat'):
            print(f"-- Start FEAT on subject {s}")

            start = time.perf_counter()
            this_dir = directory + f'/{s}'
            # update_outputdir_and_run_feat(this_dir)

            print(this_dir)
            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f'-- Done {elapse}')
        else:
            print(f'Already processed subject {s}')


if __name__ == "__main__":
    # # step 1; change the smoothing kernel in the fsf files
    # find_and_update_fsf_files(directory)

    # # step 2 run feat
    # update_outputdir_and_run_feat(directory)

    if len(sys.argv) != 2:
        print("Usage: python hcp_tfmri_glm.py <num_cpus>")
        sys.exit(1)

    T = pd.read_csv('/data/tge/Tian/HCP_img/subj_list/HCP203_test_set.tsv', delimiter='\t')

    # num_cpus = 100  # Get CPUs from SLURM
    # set_indices = list(range(num_cpus))  # Modify based on the number of parallel tasks

    # # Set up multiprocessing Pool
    # with Pool(num_cpus) as pool:
    #     pool.map(run_feat_all, set_indices)  # Distribute tasks across CPUs

    # print("Processing complete.")


    s = T.participant_id[int(sys.argv[1])-1]
    print(f"-- Start FEAT on subject {s}")

    start = time.perf_counter()
    update_outputdir_and_run_feat(directory + f'/{s}')
    finish = time.perf_counter()
    elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
    print(f'-- Done {elapse}')

    

    
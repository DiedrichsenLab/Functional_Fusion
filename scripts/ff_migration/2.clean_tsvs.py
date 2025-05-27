# script to cleanup the tsv files for each dataset according to the new naming convention in docs/task_naming.tsv

import os
import pandas as pd



def clean_hcp_tfmri(dir,subject_list,task_map):
    """
    1.Clean up the HCP tfMRI dataset tsv files, maintain task_name,cond_name,run,reg_id and remove everything else.
    2.add matching task_code and conde_code columns using the taks_naming.tsv file.


    Parameters:
    - dir (str): Path to the hcp tfMRI dataset directory (after copying)
    - subject_list (list): List of subject folder names to clean up
    """

    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'HCP-task']

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'HCP_tfMRI', 'derivatives', subject, 'estimates')
        if not os.path.exists(subject_dir):
            print(f"{subject_dir} not found.")
            continue
        for session_name in os.listdir(subject_dir):
            session_path = os.path.join(subject_dir, session_name)
            if not os.path.isdir(session_path):
                continue # skip if not a directoy needed since there is a dsstore file

            for fname in os.listdir(session_path):
                if fname.endswith('.tsv'):
                    tsv_path = os.path.join(session_path, fname)
                    df = pd.read_csv(tsv_path, sep='\t')

                    # keep some columns that are fine
                    cols_to_keep = ['task_name', 'cond_name', 'run', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'cond_name', 'task_code', 'cond_code']],
                                  on=['task_name', 'cond_name'], how='left')

                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")



if __name__=='__main__':
    # dirs
    datashare_dir = 'Y:/data/'
    if not os.path.exists(datashare_dir):
        datashare_dir = '/cifs/diedrichsen/data/'

    ff_dir = f'{datashare_dir}/FunctionalFusion_new/'

    # task_naming.tsv file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    taskmap_file = os.path.join(script_dir, '..', '..', 'docs', 'task_naming.tsv')
    taskmap_file = os.path.abspath(taskmap_file)
    task_map = pd.read_csv(taskmap_file, sep='\t')

    # what to cleanup
    subject_list = ['sub-101309']
    clean_hcp_tfmri(ff_dir, subject_list, task_map)



    pass
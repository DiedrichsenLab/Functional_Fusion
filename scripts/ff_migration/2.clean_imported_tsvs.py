# script to cleanup the tsv files for each dataset according to the new naming convention in docs/task_naming.tsv

import os
import pandas as pd
import numpy as np



def clean_hcp_tfmri(dir,subject_list,task_map):
    """
    1.Clean up the HCP tfMRI dataset tsv files, maintain task_name,cond_name,run,reg_id and remove everything else.
    2.add matching task_code and conde_code columns using the taks_naming.tsv file.


    Parameters:
    - dir (str): Path to the hcp tfMRI dataset directory (after copying)
    - subject_list (list): List of subject folder names to clean up
    """

    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'HCP']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'HCP_tfMRI', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'HCP_tfMRI', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task_name', 'cond_name', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'cond_name', 'task_code', 'cond_code']],
                                  on=['task_name', 'cond_name'], how='left')
                    
                    df['half'] = np.where(df['run'] <= 7, 1, 2)

                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")

def clean_MDTB(dir,subject_list,task_map):
    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'MDTB']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'MDTB', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'MDTB', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task_name', 'cond_name', 'reg_id','instruction','common']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # strip white spaces in task_name and cond_name
                    df['task_name'] = df['task_name'].str.strip()
                    df['cond_name'] = df['cond_name'].str.strip()

                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'cond_name', 'task_code', 'cond_code']],
                                  on=['task_name', 'cond_name'], how='left')
                    
                    df['half'] = np.where(df['run'] < 9, 1, 2)

                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")

def clean_pontine(dir,subject_list,task_map):

    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'Pontine']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'Pontine', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'Pontine', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','taskName', 'inst', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # rename  taskName to task_name and inst to instruction
                    df.rename(columns={'taskName': 'task_name', 'inst': 'instruction'}, inplace=True)

                    # strip white spaces in task_name
                    df['task_name'] = df['task_name'].str.strip()



                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'task_code', 'cond_code']],
                                  on=['task_name'], how='left')
                    
                    df['half'] = np.where(df['run'] < 9, 1, 2)


                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")
def clean_language(dir,subject_list,task_map):
    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'Language']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'Language', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'Language', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','taskName', 'inst', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # rename  taskName to task_name and inst to instruction
                    df.rename(columns={'taskName': 'task_name', 'inst': 'instruction'}, inplace=True)

                    # strip white spaces in task_name
                    df['task_name'] = df['task_name'].str.strip()


                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'task_code', 'cond_code']],
                                on=['task_name'], how='left')
                    
                    df['half'] = np.where(df['run'] % 2 == 1, 1, 2)

                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")



def clean_nishimoto(dir,subject_list,task_map):

    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'Nishimoto']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'Nishimoto', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'Nishimoto', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task_name', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # strip white spaces in task_name
                    df['task_name'] = df['task_name'].str.strip()

                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'task_code', 'cond_code']],
                                  on=['task_name'], how='left')
                    
                    df['half'] = 2 - (df.run < (len(np.unique(df.run)) / 2 + 1))


                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")

def clean_demand(dir,subject_list,task_map):

    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'Demand']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'Demand', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'Demand', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task','cond_name', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # rename task to task_name
                    df.rename(columns={'task': 'task_name'}, inplace=True)

                    # strip white spaces in task_name and cond_name
                    df['task_name'] = df['task_name'].astype(str).str.strip()
                    df['cond_name'] = df['cond_name'].astype(str).str.strip()
                    task_map['task_name'] = task_map['task_name'].astype(str).str.strip()
                    task_map['cond_name'] = task_map['cond_name'].astype(str).str.strip()
                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'cond_name', 'task_code', 'cond_code']],
                                  on=['task_name', 'cond_name'], how='left')
                    
                    
                    df['half'] = (df.run % 2) + 1


                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")

def clean_somatotopic(dir,subject_list,task_map):
    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'Somatotopic']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'Somatotopic', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'Somatotopic', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task','cond_name', 'reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]
                    

                    # rename task to task_name
                    df.rename(columns={'task': 'task_name'}, inplace=True)

                    # strip white spaces in task_name and cond_name
                    df['task_name'] = df['task_name'].str.strip()
                    df['cond_name'] = df['cond_name'].str.strip()




                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'cond_name', 'task_code', 'cond_code']],
                                  on=['task_name', 'cond_name'], how='left')
                    
                    
                    df['half'] = (df.run % 2) + 1


                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")

def clean_wmfs(dir,subject_list,task_map):
    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'WMFS']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'WMFS', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'WMFS', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task_name','cond_name','reg_id']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    #remove white spaces in task_name and cond_name
                    df['task_name'] = df['task_name'].str.strip()
                    df['cond_name'] = df['cond_name'].str.strip()


                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'cond_name', 'task_code', 'cond_code']],
                                  on=['task_name', 'cond_name'], how='left')
                    
                    
                    df['half'] = 2 - (df.run < 3)


                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")

def clean_social(dir,subject_list,task_map):
    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'Social']
    task_map = task_map.drop_duplicates(subset=['task_name', 'cond_name'])

    participants_tsv = os.path.join(dir,'Social', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'Social', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task_name', 'reg_id','Instruction']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # rename Instruction to instruction
                    df.rename(columns={'Instruction': 'instruction'}, inplace=True)

                    # strip white spaces in task_name
                    df['task_name'] = df['task_name'].str.strip()

                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name', 'task_code', 'cond_code']],
                                  on=['task_name'], how='left')
                    
                    df['half'] = 2 - (df.run < 5)


                    # overwrite the tsv file
                    df.to_csv(tsv_path, sep='\t', index=False)
                    print(f"Cleaned: {tsv_path}")


def clean_ibc(dir,subject_list,task_map):
    # Filter to HCP-task dataset only
    task_map = task_map[task_map['Dataset'] == 'IBC']

    participants_tsv = os.path.join(dir,'IBC', 'participants.tsv')

    if not subject_list:
        T = pd.read_csv(participants_tsv, sep='\t')
        subject_list = T['participant_id'].tolist()

    # Loop through each subject
    for subject in subject_list:
        subject_dir = os.path.join(dir, 'IBC', 'derivatives','ffimport', subject, 'func')
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
                    cols_to_keep = ['run','task_name', 'cond_name','reg_id','half']
                    df = df[[col for col in cols_to_keep if col in df.columns]]

                    # add task_code and cond_code
                    df = df.merge(task_map[['task_name','cond_name', 'task_code', 'cond_code']],
                                  on=['task_name','cond_name'], how='left')


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
    subject_list = None
    # clean_hcp_tfmri(ff_dir, subject_list, task_map)
    # clean_nishimoto(ff_dir, subject_list, task_map)
    # clean_pontine(ff_dir, subject_list, task_map)
    # clean_language(ff_dir, subject_list, task_map)
    # clean_MDTB(ff_dir, subject_list, task_map)
    # clean_demand(ff_dir, subject_list, task_map)
    # clean_somatotopic(ff_dir, subject_list, task_map)
    # clean_wmfs(ff_dir, subject_list, task_map)
    # clean_social(ff_dir, subject_list, task_map)
    # clean_ibc(ff_dir, subject_list, task_map)



    pass
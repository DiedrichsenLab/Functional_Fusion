#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for importing the UK Biobank dataset to general format.
Created on 2/26/2024 at 4:34 PM
Author: dzhi
"""
import pandas as pd
import shutil, logging
from pathlib import Path
import os, sys, time, subprocess
from copy import deepcopy
import numpy as np
import nibabel as nb

import Functional_Fusion.dataset as ds
import Functional_Fusion.util as ut

BASE_DIR = '/data/tge/Tian/UKBB_full/imaging'
if not Path(BASE_DIR).exists():
    BASE_DIR = '/Volumes/diedrichsen_data$/data'
if not Path(BASE_DIR).exists():
    BASE_DIR = '/srv/diedrichsen/data'
if not Path(BASE_DIR).exists():
    BASE_DIR = 'Y:/data'
if not Path(BASE_DIR).exists():
    print('Cannot find data folder')

WORK_DIR = os.path.join(BASE_DIR, 'scripts')
IMG_DIR = os.path.join(BASE_DIR, 'rfMRI')
DERIVATIVES_DIR = os.path.join(BASE_DIR, 'derivatives')

# Configure the logging settingsd
logging.basicConfig(filename=WORK_DIR + '/log/download_ukb_rfMRI.log',
                    level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def is_folder_exist_and_empty(folder_path):
    ''' Check if a folder exists and empty

    Args:
        folder_path: give the folder directory

    Returns:
        whether a given folder is exist and empty
    '''
    if not os.path.exists(folder_path):
        return False

    return not any(os.scandir(folder_path))

def add_instance_to_lines(input_file, string_to_add):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    added_lines = [l.strip() + " " + string_to_add for l in lines]

    with open(input_file, 'w') as file:
        file.write('\n'.join(added_lines))

def modify_bulk(input_file, output_file):
    id_occurrences = {}

    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            # Read id-value pair / store occurence
            id_value_pair = line.strip().split()
            if len(id_value_pair) == 2:
                id_value = id_value_pair[0]
                id_occurrences[id_value] = id_occurrences.get(id_value, 0) + 1

    # Write in extracted to new bulk
    with open(output_file, 'w') as outfile:
        for line in lines:
            id_value_pair = line.strip().split()

            # Remove single entry
            if len(id_value_pair) == 2:
                id_value = id_value_pair[0]
                if id_occurrences[id_value] > 1:
                    outfile.write(line)

def split_bulk(input_file, out_dir, lines_per_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Calculate the number of files needed using ceiling division
    num_files = -(-len(lines) // lines_per_file)

    for file_num in range(num_files):
        start_index = file_num * lines_per_file
        end_index = (file_num + 1) * lines_per_file
        output_file = f"{out_dir}/{Path(input_file).stem}_split_{file_num+1}.bulk"

        with open(output_file, 'w') as output_file:
            output_file.writelines(lines[start_index:end_index])

    print(f"{num_files} files created.")

def split_tsv(input_file, out_dir, m):
    df = pd.read_csv(input_file, sep='\t')

    # Split the DataFrame into m pieces
    # chunks = np.array_split(df, m)
    chunks = [df.iloc[i:i+m] for i in range(0, len(df), m)]

    # Save each chunk to a separate TSV file with the same column name
    for i, chunk in enumerate(chunks):
        filename = f"{out_dir}/{Path(input_file).stem}_split_{i+1}.tsv"
        chunk.to_csv(filename, sep='\t', index=False)

    print(f"{len(chunks)} files created.")

def import_rfMRI_timeseries(bulk_file, img_dir=IMG_DIR):
    os.chdir(WORK_DIR)
    fetch_ukb = '/data/tge/dzhi/software/ukbfetch'
    ukb_key = 'k32568r674571.ukbkey'

    with open(bulk_file, 'r') as file:
        for line in file:
            # Strip newline characters and other whitespaces
            subj_id, instance_id = line.strip().split()
            dest_dir = os.path.join(img_dir, subj_id, instance_id)

            if not os.path.exists(dest_dir) or is_folder_exist_and_empty(dest_dir):
                tic = time.perf_counter()
                # Creat folder for current subject / instance
                Path(dest_dir).mkdir(parents=True, exist_ok=True)

                # Download this subject's data if not exist
                cmd = f'{fetch_ukb} -e{subj_id} -d{instance_id} -a{ukb_key}'
                print(f"Start downloading subject {subj_id}, {instance_id} rfMRI...")
                
                try:
                    res = subprocess.run(cmd, check=True, shell=True)
                    subprocess.run(f'mv {subj_id}_{instance_id}.zip {dest_dir}',
                                   shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print("Error:", e)
                    print("Command output:", e.stdout)
                    print("Command error:", e.stderr)

                toc = time.perf_counter()
                min, sec = divmod(toc - tic, 60)
                print(f'Done - time used {int(min):02d}:{int(sec):02d}.')

                try:
                    # Unzip and move to image folder
                    cmd = (f'unzip {dest_dir}/{subj_id}_{instance_id}.zip -d {dest_dir}'
                           f'&& rm {dest_dir}/{subj_id}_{instance_id}.zip')
                    subprocess.run(cmd, shell=True, check=True)

                    # Extract only necessary files we will use
                    subprocess.run(f'mv {dest_dir}/fMRI/rfMRI.ica '
                                   f'{dest_dir}/fMRI/rfMRI_25.dr '
                                   f'{dest_dir}/fMRI/rfMRI_100.dr {dest_dir} '
                                   f'&& rm -r {dest_dir}/fMRI',
                                   shell=True, check=True)
                except:
                    print(f'skipping extract {subj_id}_{instance_id}, the file incomplete!')
                    logging.error(f"Missing fMRI data {subj_id} {instance_id} ")

            else:
                print(f'{subj_id}_{instance_id} already exist!')


def check_rfMRI_timeseries(bulk_file, img_dir=IMG_DIR):
    os.chdir(WORK_DIR)

    with open(bulk_file, 'r') as file:
        for line in file:
            # Strip newline characters and other whitespaces
            subj_id, instance_id = line.strip().split()
            dest_dir = os.path.join(img_dir, subj_id, instance_id)

            if not os.path.exists(dest_dir):
                print(f'{subj_id} {instance_id} data is not downloaded!')
            else:
                # Real checking start here
                print(f"Checking subject {subj_id}, {instance_id} rfMRI...")
                if os.path.exists(dest_dir + f'/{subj_id}_{instance_id}.zip'):
                    # Unzip and move to image folder
                    try:
                        cmd = (f'unzip {dest_dir}/{subj_id}_{instance_id}.zip -d {dest_dir} '
                               f'&& rm {dest_dir}/{subj_id}_{instance_id}.zip')
                        res = subprocess.run(cmd, shell=True, check=True)

                        # Extract only necessary files we will use
                        subprocess.run(f'mv {dest_dir}/fMRI/rfMRI.ica '
                                    f'{dest_dir}/fMRI/rfMRI_25.dr '
                                    f'{dest_dir}/fMRI/rfMRI_100.dr {dest_dir} '
                                    f'&& rm -r {dest_dir}/fMRI',
                                    shell=True, check=True)
                        
                    except subprocess.CalledProcessError as e:
                        print("Error: ", e)
                        print("Command output: ", e.stdout)
                        print("Command error: ", e.stderr)
                elif os.path.exists(dest_dir + '/fMRI/rfMRI.ica'):
                    print(f'No zip file but already unzipped for {subj_id} {instance_id}')
                    
                    try:
                        # Extract only necessary files we will use
                        subprocess.run(f'mv {dest_dir}/fMRI/rfMRI.ica '
                                    f'{dest_dir}/fMRI/rfMRI_25.dr '
                                    f'{dest_dir}/fMRI/rfMRI_100.dr {dest_dir} '
                                    f'&& rm -r {dest_dir}/fMRI',
                                    shell=True, check=True)
                    except subprocess.CalledProcessError as e:
                        print ("Error: ", e)

                else:
                    if os.path.exists(dest_dir + '/fMRI/unusable'):
                        print(f'Real missing - {subj_id}_{instance_id}, the file incomplete!')
                    else:
                        logging.error(f"Missing fMRI data {subj_id} {instance_id} ")


def get_incomplete_list(bulk_file):
    out_file = WORK_DIR + '/participants_missing.bulk'

    with open(out_file, 'w') as outfile:
        with open(bulk_file, 'r') as file:
            for line in file:
                subj_id, instance_id = line.strip().split()
                dest_dir = os.path.join(IMG_DIR, subj_id, 
                                        instance_id, 'rfMRI.ica')

                if not os.path.exists(dest_dir):
                    outfile.write(f'{subj_id} {instance_id}\n')  

def get_complete_list(tsv_file, out_file):
    T = pd.read_csv(tsv_file, sep='\t')
    new_id = []

    for i, s in enumerate(T.participant_id):
        dir_1 = os.path.join(IMG_DIR, str(s), '20227_2_0', 'rfMRI.ica',
                             'filtered_func_data_clean.nii.gz')
        dir_2 = os.path.join(IMG_DIR, str(s), '20227_3_0', 'rfMRI.ica', 
                             'filtered_func_data_clean.nii.gz')
        if os.path.exists(dir_1) and os.path.exists(dir_2):
            new_id.append(s)

    df = pd.DataFrame(new_id, columns=['participant_id'])
    df.to_csv(out_file, sep='\t', index=False)

def makeAndSplit_single_scan_subjlist_bulk(input_file, out_dir, lines_per_file=100):

    with open(input_file, 'r') as file:
        lines = file.readlines()  # Read all lines into a list
        lines = [line.strip() + ' 20227_2_0' for line in lines]

    # Calculate the number of files needed using ceiling division
    num_files = -(-len(lines) // lines_per_file)

    for file_num in range(num_files):
        start_index = file_num * lines_per_file
        end_index = (file_num + 1) * lines_per_file
        output_file = f"{out_dir}/{Path(input_file).stem}_split_{file_num+1}.bulk"

        with open(output_file, 'w') as o_file:
            o_file.writelines('\n'.join(lines[start_index:end_index]))

    print(f"{num_files} files created.")



if __name__ == "__main__":
    # b_file = '/ukb674571_rsfmri-return.bulk'
    # split_bulk(BASE_DIR + b_file, WORK_DIR + '/subj_list', 100)
    b_file = '/subj_list/validation/participants_validation.txt'
    # add_instance_to_lines(WORK_DIR + b_file, '20227_2_0')
    # split_bulk(WORK_DIR + b_file, WORK_DIR + '/subj_list/validation', 50)
    split_tsv(WORK_DIR + b_file, WORK_DIR + '/subj_list/validation', 50)

    # get_incomplete_list(BASE_DIR + '/ukb674571_rsfmri-return.bulk')
    # get_complete_list(BASE_DIR + '/participants_filtered.tsv', 
    #                   BASE_DIR + '/participants_filtered_final.tsv')
    # check_rfMRI_timeseries(WORK_DIR + '/participants_missing.bulk')
    # makeAndSplit_single_scan_subjlist_bulk(WORK_DIR + '/subjlist_eur_unrel_65y_woRep_woF.txt',
    #                                        WORK_DIR + '/subj_list', lines_per_file=100)

    # # Run script from command line
    if len(sys.argv) != 2:
        print("Usage: python import_ukb.py <bulk files>")
        sys.exit(1)

    # Call import function
    import_rfMRI_timeseries(WORK_DIR + '/subj_list/' + sys.argv[1],
                            img_dir=IMG_DIR)
    # import_rfMRI_timeseries(WORK_DIR + '/subj_list/discovery/participants_discovery_split_1.bulk',
    #                         img_dir=IMG_DIR)

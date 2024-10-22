import os
import pandas as pd

base_path = "/cifs/diedrichsen/data/FunctionalFusion/MDTB/derivatives"
task_table_path = "/home/ROBARTS/barafat/Documents/task_table.tsv"

task_table = pd.read_csv(task_table_path, sep="\t")

# Create a dictionary to map each condition to its code
condition_map = {}
for _, row in task_table.iterrows():
    cond = row['cond']
    cond_code = row['code']
    
    # if condition contains dashes split it into sub-conditions
    if "-" in cond:
        sub_conditions = cond.split("-")
    else:
        sub_conditions = [cond]

    # same for codes
    if "-" in cond_code:
        sub_codes = cond_code.split("-")
    else:
        sub_codes = [cond_code]

    
    # Map each sub-condition to the cond_code
    for i, sub_cond in enumerate(sub_conditions):
        condition_map[sub_cond] = sub_codes[i]

# Iterate through each "sub-x" folder
for sub_folder in os.listdir(base_path):
    sub_path = os.path.join(base_path, sub_folder)
    
    # Check if it contains a "data" folder
    data_folder = os.path.join(sub_path, "data")
    if not os.path.isdir(data_folder):
        continue
    
    # Iterate through each TSV file in the "data" folder
    for file in os.listdir(data_folder):
        if file.endswith(".tsv"):
            file_path = os.path.join(data_folder, file)
            
            # Load the TSV file
            tsv_data = pd.read_csv(file_path, sep="\t")
            
            # Check if "cond_name" column exists
            if "cond_name" in tsv_data.columns:
                
                # Add "cond_code" column by mapping the cond_name to condition_map
                tsv_data['cond_code'] = tsv_data['cond_name'].apply(lambda x: condition_map.get(x, ''))

                #save
                tsv_data.to_csv(file_path, sep="\t", index=False)

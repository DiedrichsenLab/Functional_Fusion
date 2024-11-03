import os
import zipfile
import shutil
import nibabel as nb



directory = 'Y:\\data\\ExternalOpenData\\HCP_UR100_tfMRI_new'
if not os.path.exists(directory):
    directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_new'

# step 1: Delete .md5 files and extract .zip files into individual folders
def unzip_clean(directory):
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)

        # Delete all .md5 files
        if item.endswith('.md5'):
            os.remove(item_path)

        # If it's a zip file, extract it into a folder with the same name and delete the zip
        elif item.endswith('.zip'):
            extract_dir = os.path.join(directory, os.path.splitext(item)[0])
            os.makedirs(extract_dir, exist_ok=True)
            
            with zipfile.ZipFile(item_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)
            
            os.remove(item_path)

# step 2: Group all folders by participant ID
def group_participant_folders(directory):
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)

        if os.path.isdir(item_path) and '_' in item:
            participant_id = item.split('_')[0]
            participant_folder = os.path.join(directory, f'sub-{participant_id}')

            # Create a directory for the participant if it doesn't exist
            os.makedirs(participant_folder, exist_ok=True)
            
            # Move the extracted folder into the participant's main folder
            shutil.move(item_path, participant_folder)


# step 3: Reorganize anat stuff
def process_anat(directory):
    for participant_id in os.listdir(directory):
        participant_path = os.path.join(directory, participant_id)
        
        if os.path.isdir(participant_path):
            participant_number = participant_id.split('-')[1]
            anat_folder = os.path.join(participant_path, 'anat')
            os.makedirs(anat_folder, exist_ok=True)

            structural_preproc_folder = os.path.join(participant_path, f'{participant_number}_3T_Structural_preproc', participant_number)
            
            if os.path.isdir(structural_preproc_folder):
                mni_folder = os.path.join(structural_preproc_folder, 'MNINonLinear')
                
                if os.path.isdir(mni_folder):
                    # Copy the 'xfms' folder if it exists and isn't already in 'anat'
                    xfms_folder = os.path.join(mni_folder, 'xfms')
                    if os.path.isdir(xfms_folder) and not os.path.exists(os.path.join(anat_folder, 'xfms')):
                        shutil.copytree(xfms_folder, os.path.join(anat_folder, 'xfms'))

                    # Decompress and copy 'T1w.nii.gz' and 'BiasField.nii.gz' into the 'anat' folder
                    t1w_file = os.path.join(mni_folder, 'T1w.nii.gz')
                    bias_field_file = os.path.join(mni_folder, 'BiasField.nii.gz')
                    
                    if os.path.isfile(t1w_file):
                        t1w_img = nb.load(t1w_file)
                        t1w_output_path = os.path.join(anat_folder, 'T1w.nii')
                        nb.save(t1w_img, t1w_output_path)
                    
                    if os.path.isfile(bias_field_file):
                        bias_field_img = nb.load(bias_field_file)
                        bias_field_output_path = os.path.join(anat_folder, 'BiasField.nii')
                        nb.save(bias_field_img, bias_field_output_path)

                    # Move and rename 'fsaverage_LR32k' to 'SurfaceWB' change this name?
                    fsaverage_folder = os.path.join(mni_folder, 'fsaverage_LR32k')
                    final_surface_wb_folder = os.path.join(participant_path, 'SurfaceWB')

                    if os.path.isdir(fsaverage_folder):
                        shutil.move(fsaverage_folder, final_surface_wb_folder)

                # Remove the structural preproc folder
                shutil.rmtree(structural_preproc_folder, ignore_errors=True)

        
#  step 4: Reorganize func stuff
def process_func(directory):
    for subject_id in os.listdir(directory):  # Iterate over subject directories
        subject_path = os.path.join(directory, subject_id)
        subject_number = subject_id.split('-')[1]
        if os.path.isdir(subject_path):
            func_directory = os.path.join(subject_path, 'func')
            os.makedirs(func_directory, exist_ok=True)

            # Iterate over session directories that need preprocessing
            for session in os.listdir(subject_path):
                # Skip non-functional directories like 'anat', 'func', and 'SurfaceWB'
                if session in ['anat', 'func', 'SurfaceWB']:
                    continue
                
                # Process directories with '3T_tfMRI' in the name
                if '3T_tfMRI' in session:
                    session_name = session.split('_')[-2]
                    session_path = os.path.join(subject_path, session,subject_number, 'MNINonLinear', 'Results')
                    if os.path.isdir(session_path):                        
                        # Create a session directory inside 'func'
                        session_func_dir = os.path.join(func_directory, f'ses-{session_name}')
                        os.makedirs(session_func_dir, exist_ok=True)

                        for run_name in os.listdir(session_path):  
                            if 'LR' in run_name or 'RL' in run_name:  # Only process runs with 'LR' or 'RL' (ignore folder that has 2ndlevel fsf)
                                runs_found = True
                                run_path = os.path.join(session_path, run_name)
                                if os.path.isdir(run_path):
                                    # Clean unwanted files
                                    for file in os.listdir(run_path):
                                        file_path = os.path.join(run_path, file)
                                        if file.endswith((
                                            '.dtseries.nii', '.func.gii', 'SBRef_dc.nii.gz',
                                            'gdc_dc.nii.gz', 'MSMAll.dtseries.nii'
                                        )) or 'RibbonVolumeToSurfaceMapping' in file_path:

                                            if os.path.isfile(file_path):
                                                os.remove(file_path)
                                            elif os.path.isdir(file_path):
                                                shutil.rmtree(file_path, ignore_errors=True)

                                    # Move the run folder to the session directory in 'func'
                                    new_run_path = os.path.join(session_func_dir, run_name)
                                    shutil.move(run_path, new_run_path)

                shutil.rmtree(os.path.join(subject_path, session), ignore_errors=True)

unzip_clean(directory)
group_participant_folders(directory)
process_anat(directory)
process_func(directory)
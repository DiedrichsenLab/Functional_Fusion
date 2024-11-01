import os
import zipfile
import shutil


directory = 'Y:\\data\\ExternalOpenData\\HCP_UR100_tfMRI_new'

# step 1: Delete .md5 files and extract .zip files into individual folders
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
for participant_id in os.listdir(directory):
    participant_path = os.path.join(directory, participant_id)
    participant_number = participant_id.split('-')[1]
    
    if os.path.isdir(participant_path):
        anat_folder = os.path.join(participant_path, 'anat')
        os.makedirs(anat_folder, exist_ok=True)

        structural_preproc_folder = os.path.join(directory, participant_id, f'{participant_number}_3T_Structural_preproc',participant_number)
        
            
        if os.path.isdir(structural_preproc_folder):
            mni_folder = os.path.join(structural_preproc_folder, 'MNINonLinear')
            
            if os.path.isdir(mni_folder):
                xfms_folder = os.path.join(mni_folder, 'xfms')
                if os.path.isdir(xfms_folder) and not os.path.join(anat_folder, 'xfms'):
                    shutil.copytree(xfms_folder, os.path.join(anat_folder, 'xfms'))
                
                t1w_file = os.path.join(mni_folder, 'T1w.nii.gz')
                bias_field_file = os.path.join(mni_folder, 'BiasField.nii.gz')
                #unzip both
                with zipfile.ZipFile(t1w_file , 'r') as zip_ref:
                    zip_ref.extractall(anat_folder)
                with zipfile.ZipFile(bias_field_file , 'r') as zip_ref:
                    zip_ref.extractall(anat_folder)

                # unzipped files paths
                t1w_file = os.path.join(anat_folder, 'T1w.nii')
                bias_field_file = os.path.join(anat_folder, 'BiasField.nii')
                if os.path.isfile(t1w_file):
                    shutil.copy2(t1w_file, anat_folder)
                if os.path.isfile(bias_field_file):
                    shutil.copy2(bias_field_file, anat_folder)

                surface_wb_folder = os.path.join(participant_path, 'SurfaceWB')

                fsaverage_folder = os.path.join(structural_preproc_folder, 'MNINonLinear', 'fsaverage_LR32k')

                if os.path.isdir(fsaverage_folder):
                    shutil.move(fsaverage_folder, surface_wb_folder)

            # delete structural_preproc_folder
            os.removedirs(structural_preproc_folder)
        




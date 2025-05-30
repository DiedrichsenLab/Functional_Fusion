import os
import subprocess
import cerebellum_isolate.isolate as iso
import pandas as pd
import numpy as np



directory = 'Y:\\data\\ExternalOpenData\\HCP_UR100_tfMRI_full'
if not os.path.exists(directory):
    directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_tfMRI_full'
if not os.path.exists(directory):
    directory = '/Volumes/diedrichsen_data$/data/ExternalOpenData/HCP_UR100_tfMRI_full'

pinfo = pd.read_csv(f"{directory}/participants.tsv", sep="\t")
participants = pinfo["participant_id"]

# step 1: change the smoothing kernel in the fsf files
def run_isolate(sn=None):
    """
    Update the `set fmri(smooth)` value in the specified .fsf file to 0.
    """
    if sn is None:
        sn = np.arange(len(participants))
    for part in participants[sn]:
        res_dir = f'{directory}/{part}/suit' 
        anat_file = f'{res_dir}/T1W.nii'
        iso.isolate(anat_file,result_folder=res_dir)
        print(f'Isolated participant {part}')


if __name__ == "__main__":
    run_isolate()



    

    
# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
base_dir = '/Volumes/diedrichsen_data$/data'
orig_dir = base_dir + '/Cerebellum/super_cerebellum'
target_dir = base_dir + '/FunctionalFusion/MDTB'


def import_suit(source_dir,dest_dir,anat_name,participant_id):
    """

    Args:
        source_dir (str): source directory 
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file (submitted to ) 
        participant_id (str): ID of participant
    """
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = f'/c1{anat_name}.nii' 
    dest = f'/{participant_id}_label-GMc_probseg.nii'
    shutil.copyfile(source_dir+src,dest_dir+dest)
    src = f'/c2{anat_name}.nii' 
    dest = f'/{participant_id}_label-WMc_probseg.nii'
    shutil.copyfile(source_dir+src,dest_dir+dest)
    # Masking image 
    src = f'/maskbrainSUITGrey.nii' 
    dest = f'/{participant_id}_desc-cereb_mask.nii'
    shutil.copyfile(source_dir+src,dest_dir+dest)
    # Deformation map file
    src = f'/y_{anat_name}_suitdef.nii' 
    dest = f'/{participant_id}_space-SUIT_xfm.nii'
    shutil.copyfile(source_dir+src,dest_dir+dest)


if __name__ == "__main__":
    T= pd.read_csv(target_dir + '/participants.tsv',delimiter='\t')
    for s in T.participant_id:
        old_id = s.replace('sub-','s',1)
        dir1 = orig_dir + f'/sc1/suit/anatomicals/{old_id}'
        dir2 = target_dir + f'/derivatives/{s}/suit'
        import_suit(dir1,dir2,'anatomical',s)
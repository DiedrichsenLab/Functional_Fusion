# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
base_dir = '/Volumes/diedrichsen_data$/data'
orig_dir = base_dir + '/Cerebellum/super_cerebellum'
target_dir = base_dir + '/FunctionalFusion/MDTB'


def import_suit(source_dir,dest_dir,anat_name,participant_id):
    """ Imports a suit folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory 
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file (submitted to ) 
        participant_id (str): ID of participant
    """
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src=[]
    dest =[] 
    src.append(f'/c1{anat_name}.nii') 
    dest.append(f'/{participant_id}_label-GMc_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/{participant_id}_label-WMc_probseg.nii')
    src.append(f'/maskbrainSUITGrey.nii')
    dest.append(f'/{participant_id}_desc-cereb_mask.nii')
    src.append(f'/y_{anat_name}_suitdef.nii')
    dest.append(f'/{participant_id}_space-SUIT_xfm.nii')
    for i in range(len(src)): 
        try: 
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except: 
            print('skipping ' + src[i])


def import_anat(source_dir,dest_dir,anat_name,participant_id): 
    """ Imports a anatomy folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory (anatomical)
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file 
        participant_id (str): ID of participant
    """
    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    src=[]
    dest =[] 
    src.append(f'/{anat_name}.nii')
    dest.append(f'/{participant_id}_T1w.nii')
    src.append(f'/c1{anat_name}.nii')
    dest.append(f'/{participant_id}_label-GM_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/{participant_id}_label-WM_probseg.nii')
    src.append(f'/c3{anat_name}.nii')
    dest.append(f'/{participant_id}_label-CSF_probseg.nii')
    src.append(f'/y_{anat_name}.nii')
    dest.append(f'/{participant_id}_space-MNI_xfm.nii')
    for i in range(len(src)): 
        try: 
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except: 
            print('skipping ' + src[i])

def import_freesurfer(source_dir,dest_dir,old_id,new_id):
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src=[]
    dest =[] 
    src.append(f'/{old_id}.L.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_pial.surf.nii')
    src.append(f'/{old_id}.L.white.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_white.surf.nii')
    src.append(f'/{old_id}.R.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_pial.surf.nii')
    src.append(f'/{old_id}.R.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_white.surf.nii')
    src.append(f'/{old_id}.L.sulc.32k.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_sulc.shape.nii')
    src.append(f'/{old_id}.R.sulc.32k.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_sulc.shape.nii')
    for i in range(len(src)): 
        try: 
            shutil.copyfile(source_dir+src[i],dest_dir+dest[i])
        except: 
            print('skipping ' + src[i])

def import_spm_glm(source_dir,dest_dir,new_id):
    



if __name__ == "__main__":
    T= pd.read_csv(target_dir + '/participants.tsv',delimiter='\t')
    for s in T.participant_id:
        old_id = s.replace('sub-','s',1)
        # dir1 = orig_dir + f'/sc1/suit/anatomicals/{old_id}'
        # dir2 = target_dir + f'/derivatives/{s}/suit'
        # import_suit(dir1,dir2,'anatomical',s)
        # dir1 = orig_dir + f'/sc1/anatomicals/{old_id}'
        # dir2 = target_dir + f'/derivatives/{s}/anat'
        # import_anat(dir1,dir2,'anatomical',s)
        dir1 = orig_dir + f'/sc1/surfaceWB/{old_id}'
        dir2 = target_dir + f'/derivatives/{s}/anat'
        import_freesurfer(dir1,dir2,old_id,s)

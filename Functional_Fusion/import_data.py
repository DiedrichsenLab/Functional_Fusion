import shutil
import pandas as pd
from pathlib import Path
import numpy as np
import scipy.io as sio
import mat73
from copy import deepcopy

def import_suit(source_dir, dest_dir, anat_name, participant_id):
    """
    Imports a suit folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file (submitted to)
        participant_id (str): ID of participant
    """

    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    # Cerebellar Gray matter and white matter probability maps (full anatomical space) 
    src.append(f'/c1{anat_name}.nii')
    dest.append(f'/{participant_id}_label-GMc_probseg.nii')
    src.append(f'/c2{anat_name}.nii')
    dest.append(f'/{participant_id}_label-WMc_probseg.nii')
    
    # Gray-matter mask image in functional space
    src.append('/maskbrainSUITGrey.nii')
    dest.append(f'/{participant_id}_desc-cereb_mask.nii')
    
    # Deformation file: See Documentation for import 
    src.append(f'/y_{anat_name}_suitdef.nii')
    dest.append(f'/{participant_id}_space-SUIT_xfm.nii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir + src[i], dest_dir + dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])


def import_anat(source_dir, dest_dir, anat_name, participant_id):
    """
    Imports a anatomy folder into a BIDS/derivtive structure

    Args:
        source_dir (str): source directory (anatomical)
        dest_dir (str): destination directory
        anat_name (str): Name of the anatomical main file
        participant_id (str): ID of participant
    """

    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)

    src = []
    dest = []
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
            shutil.copyfile(source_dir + src[i], dest_dir + dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])

def import_freesurfer(source_dir, dest_dir, old_id, new_id):
    """
    Imports the output of a freesurfer reconstruction (and subsequent
    workbench import).
    Args:
        source_dir (str): Directory of the SPM GLM
        dest_dir (str): Destination directory for that subject / session
        old_id (str): Old subject name
        new_id (str): New name for the subject
    """

    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    src.append(f'/{old_id}.L.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_pial.surf.gii')
    src.append(f'/{old_id}.L.white.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_white.surf.gii')
    src.append(f'/{old_id}.R.pial.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_pial.surf.gii')
    src.append(f'/{old_id}.R.white.32k.surf.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_white.surf.gii')
    src.append(f'/{old_id}.L.sulc.32k.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-L_sulc.shape.gii')
    src.append(f'/{old_id}.R.sulc.32k.shape.gii')
    dest.append(f'/{new_id}_space-32k_hemi-R_sulc.shape.gii')
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i], dest_dir+dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])


def import_spm_glm(source_dir, dest_dir, sub_id, sess_id, info_dict):
    """
    Imports the output of the SPM GLM with an SPM_info.mat
    structure into BIDS deriviatie (Functional Fusion) framework.
    It assumes that a single GLM corresponds to single session.

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that
                           subject / session
        new_id (_type_): New name for the subject
        info_dict (_type_): Dictionary with the old field names and the
			                new field names for the information
    """

    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    src = []
    dest = []
    # Generate new dictionary from SPM info
    D = mat73.loadmat(source_dir + '/SPM_info.mat')
    T = {}
    for i in info_dict.items():
        series=D[i[0]]
        if isinstance(series[0],list):
            series=[series[i][0] for i in range(len(series))]
        T[i[1]]=series
    N = len(T[i[1]])
    if 'reg_id' not in T.keys():
        n = sum(T['run'] == 1)
        T['reg_num'] = np.arange(N)
        T['reg_id'] = T['reg_num'] % n

    # Ensure that run number is an integer value
    T['run'] = [int(j) for j in T['run']]
    D = pd.DataFrame(T)
    D.to_csv(dest_dir + f'/{sub_id}_{sess_id}_reginfo.tsv', sep='\t')

    # Prepare beta files for transfer
    src = []
    dest = []
    for i in range(N):
        src.append(f'/beta_{i+1:04d}.nii')
        dest.append(f'/{sub_id}_{sess_id}_run-{D.run[i]:02}_' +
                    'reg-{D.reg_id[i]:02d}_beta.nii')
    # Mask
    src.append('/mask.nii')
    dest.append(f'/{sub_id}_{sess_id}_mask.nii')

    # ResMS
    src.append('/resms.nii')
    dest.append(f'/{sub_id}_{sess_id}_resms.nii')

    # Copy those files over
    for i in range(len(src)):
        try:
            shutil.copyfile(source_dir+src[i], dest_dir+dest[i])
        except FileNotFoundError:
            print('skipping ' + src[i])


def import_spm_designmatrix(source_dir, dest_dir, sub_id, sess_id):
    """
    Imports the SPM design matrix for optimal contrast recombination
    at a later stage. Because python gives some errors when trying to
    read an SPM.mat structure, this requires the design matrix
    information to be extracted from the SPM.mat before, using the
    following matlab code (for every subject):
 
    load('SPM.mat');
    X = SPM.xX.xKXs.X
    save design_matrix.mat X

    See readme for output structure.
    Args:
        source_dir (_type_): Directory of the SPM GLM
        dest_dir (_type_): Destination directory for that
                           subject / session
        sub_id (_type_): New name for the subject
        sess_id (_type_): ID of the session to import
    """

    X = sio.loadmat(source_dir + '/design_matrix.mat')
    DM = X['X']
    filename = dest_dir + f'/{sub_id}_{sess_id}_designmatrix.npy'
    np.save(filename, DM)


def create_reginfo(dest_dir, participant_id, ses_id='ses-rest', reginfo_general='sub-02'):
    """ Creates a reginfo.tsv file for a given participant
    Args:
        reginfo_general (str): file name of general reginfo file
        dest_dir (str): destination directory
        participant_id (str): ID of participant
        ses_id (str): ID of session

        N.B. General reginfo file for rest is in the following format and needs to be created for ONE subject only in the subject's estimates folder:

            run	timepoint	task	time_id
            1	T0001	    rest	1
            1	T0002	    rest	2
            1	T0003	    rest	3
            1	T0004	    rest	4
                    .
                    .
                    .
            2	T0601	    rest	601
            2	T0602	    rest	602
            2	T0603	    rest	603
            2	T0604	    rest	604
                    .
                    .

    """

    # Import general info
    reginfo_general_file = f'{reginfo_general}/estimates/{ses_id}/{reginfo_general}_{ses_id}_reginfo.tsv'
    info = pd.read_csv(dest_dir.split(participant_id)[0] + reginfo_general_file, sep='\t')

    print(f'Creating reginfo for {participant_id}')

    # Ammend the reginfo.tsv file from the general file
    reginfo = deepcopy(info)


    # Make folder
    dest = dest_dir + \
        f'/{participant_id}_{ses_id}_reginfo.tsv'
    Path(dest).parent.mkdir(parents=True, exist_ok=True)

    # Save reginfo.tsv file
    reginfo.to_csv(dest, sep='\t', index=False)


def import_rest(src, dest, sub_id, ses_id, info_dict, mask_file=None):
    """Imports the resting state files
       into a BIDS/derivative structure
    Args:
        src (str): source name structure
        dest (str): destination name structure
        sub_id (str): ID of participant
        ses_id (str): ID of session
        info_dict (dict): Dictionary with run information and name of subject that stores the pre-created general timeseries reginfo.tsv file
        fix (bool): If True, then import the data that is FIX cleaned
    """
    run_names = info_dict['runs']

    # get source directory with Path
    source_dir = Path(src).parent
    dest_dir = Path(dest).parent

    # Make the destination directory
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    for run in run_names:

        # move data into the corresponding session folder
        src_file = src.format(run=run)
        dest_file = dest.format(run=run)

        # try:
        #     shutil.copyfile(src_file,
        #                     dest_file)
        # except:
        #     print('skipping ' + src)
        
        # # Make reginfo
        # create_reginfo(dest_dir, sub_id, ses_id=ses_id, reginfo_general=info_dict['reginfo_general'])

    
    # import mask
    if mask_file is None:
        mask_file = str(source_dir) + f'/{sub_id}_{ses_id}_mask.nii'
    dest_file = f'/{sub_id}_{ses_id}_mask.nii'
    try:
        shutil.copyfile(mask_file,
                        str(dest_dir) + dest_file)
    except:
        print('skipping ' + source_dir + src_file)


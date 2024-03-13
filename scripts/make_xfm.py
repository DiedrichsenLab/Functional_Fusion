# Script for distorting a new SUIT surface for the MNISymC template
import numpy as np
import nibabel as nb
import SUITPy as suit
import nitools as nt
import Functional_Fusion.util as ut

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
atlas_dir = base_dir + '/Atlases'


if __name__ == "__main__":
    # reslice_SUIT()
    ut.templateflow_xfm_h5_to_nii(atlas_dir+'/tpl-MNI152NLin6Asym'+'/tpl-MNI152NLin6Asym_from-MNI152NLin2009cAsym_mode-image_xfm.h5')
    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
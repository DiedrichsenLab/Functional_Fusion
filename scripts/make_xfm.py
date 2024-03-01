# Script for distorting a new SUIT surface for the MNISymC template
import numpy as np 
import nibabel as nb
import SUITPy as suit
import nitools as nt
import Functional_Fusion.util as ff

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
atlas_dir = base_dir + '/Atlases'


if __name__ == "__main__":
    # reslice_SUIT()
    deform_suit_surfaces()

    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
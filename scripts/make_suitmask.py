# Script for importing the MDTB data set from super_cerebellum to general format.
import numpy as np 
import nibabel as nb
from atlas_map import AtlasVolumetric, AtlasMapDeform, get_data
from dataset import DataSetMDTB
import SUITPy as suit

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'
hcp_dir = base_dir + '/HCP'

def reslice_SUIT():
    adir = atlas_dir +'/tpl-MNI152NLin6AsymC'
    img_name = adir + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
    def_name = adir + '/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
    a = suit.reslice.reslice_image(img_name,def_name,interp=1,voxelsize=(2,2,2))
    X = (a.get_fdata()>0.2).astype('i2')
    b = nb.Nifti1Image(X, affine=a.affine)                            
    nb.save(b,atlas_dir + '/tpl-SUIT/tpl-SUIT_res-2_gmcmas_uncorr.nii')
    pass 


def downsample_mask():
    '''
        
    '''

if __name__ == "__main__":
    reslice_SUIT()


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
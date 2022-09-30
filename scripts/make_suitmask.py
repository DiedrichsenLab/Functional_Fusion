# Script for importing the MDTB data set from super_cerebellum to general format.
import numpy as np 
import nibabel as nb
from nibabel.processing import resample_from_to, resample_to_output
from atlas_map import AtlasVolumetric, AtlasMapDeform
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


def downsample_mask(res=2):
    '''
        Downsamples a graymatter mask to a lower functional resolution 
    '''
    adir = atlas_dir +'/tpl-MNI152NLIn2000cSymC'
    img_name = adir + '/tpl-MNISymC_res-1.nii'
    out_name = adir + f'/tpl-MNISymC_res-{res:d}.nii'
    in_img = nb.load(img_name)
    #     out_img = nb.processing.resample_from_to(in_img,)
    temp_img = resample_to_output(in_img,voxel_sizes=[res,res,res])
    X=temp_img.get_fdata()
    X = (X+np.flip(X,axis=0))/2
    X=np.array(X>0.1,dtype=np.uint8)
    out_img = nb.Nifti1Image(X,temp_img.affine)
    nb.save(out_img,out_name);

if __name__ == "__main__":
    # reslice_SUIT()
    downsample_mask(2)
    downsample_mask(3)

    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
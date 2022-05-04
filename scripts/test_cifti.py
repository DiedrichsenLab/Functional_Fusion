# Script for importing the MDTB data set from super_cerebellum to general format.
import numpy as np 
import nibabel as nb
from atlas_map import AtlasVolumetric, AtlasMapDeform, get_data
from dataset import DataSetMDTB
import surfAnalysisPy as surf

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'
hcp_dir = base_dir + '/HCP'

def explore_cifti():
    dir = hcp_dir + '/derivatives/100307/estimates'
    A=nb.load(dir+'/sub-100307_ses-01_task-rest_space-fsLR32k_run-01_bold.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)

    for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
        print(idx,str(nam),slc)
    pass

def get_cereb():
    dir = hcp_dir + '/derivatives/100307/estimates'
    A=nb.load(dir+'/sub-100307_ses-01_task-rest_space-fsLR32k_run-01_bold.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)
    bmcl = bmf[bmf.name == 'CIFTI_STRUCTURE_CEREBELLUM_LEFT']
    bmcr = bmf[bmf.name == 'CIFTI_STRUCTURE_CEREBELLUM_RIGHT']
    bmc= bmcl + bmcr
    ijk = bmc.voxel
    X = np.zeros((90,50,50))
    X[ijk[:,0],ijk[:,1],ijk[:,2]] = 1
    N = nb.Nifti1Image(X,bmc.affine)
    nb.save(N,atlas_dir + '/tpl-MNI152AsymC_res-2' + '/tpl-MNI152AsymC_res-2_gmcmask.nii')
    pass


def get_cortex():
    dir = hcp_dir + '/derivatives/100307/estimates'
    A=nb.load(dir+'/sub-100307_ses-01_task-rest_space-fsLR32k_run-01_bold.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)
    bmcl = bmf[bmf.name == 'CIFTI_STRUCTURE_CORTEX_LEFT']
    bmcr = bmf[bmf.name == 'CIFTI_STRUCTURE_CORTEX_RIGHT']
    
    mask=np.zeros(32492,)
    mask[bmcr.vertex]=1
    colorM=np.array([[1,1,1,1],[0,0,1,1]])
    mask_gii = surf.map.make_label_gifti(data=mask,label_names=['medwall','cortex'],label_RGBA=colorM)
    nb.save(mask_gii,atlas_dir + '/tpl-fs32k' + '/tpl-fs32k_hemi-L_mask.label.gii')
    mask_gii = surf.map.make_label_gifti(data=mask,label_names=['medwall','cortex'],label_RGBA=colorM,anatomical_struct='CortexRight')
    nb.save(mask_gii,atlas_dir + '/tpl-fs32k' + '/tpl-fs32k_hemi-R_mask.label.gii')
    pass


if __name__ == "__main__":
    explore_cifti()


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
# Script for importing the MDTB data set from super_cerebellum to general format.
import numpy as np 
import nibabel as nb
from atlas_map import AtlasVolumetric, AtlasMapDeform, get_data3D
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
    dir2 = hcp_dir + '/derivatives/100307/estimates'
    A=nb.load(dir2+'/sub-100307_ses-01_task-rest_space-fsLR32k_run-01_bold.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)
    bmcl = bmf[bmf.name == 'CIFTI_STRUCTURE_CEREBELLUM_LEFT']
    bmcr = bmf[bmf.name == 'CIFTI_STRUCTURE_CEREBELLUM_RIGHT']
    bmc= bmcl + bmcr
    print(dir(bmc))
    ijk = bmc.voxel
    # X = np.zeros((90,50,50))
    X = np.zeros((bmc.volume_shape[0],bmc.volume_shape[1], bmc.volume_shape[2]))
    X[ijk[:,0],ijk[:,1],ijk[:,2]] = 1
    N = nb.Nifti1Image(X,bmc.affine)
    nb.save(N,atlas_dir + '/tpl-MNI152NLin6AsymC' + '/tpl-MNI152AsymC_res-2_gmcmask2.nii')
    return bmc


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

def get_ts_nii():
    """
    getting a cifti image, converting it to 4D vol and returning a list of 3D images (for each time point)
    """
    dir = hcp_dir + '/derivatives/100307/func'
    ts_cifti = nb.load(dir+'/ses-01'+'/sub-100307_ses-01_space-fsLR32k_run-01.dtseries.nii')
    # get brain axis models
    bmf = ts_cifti.header.get_axis(1)
    # get the data array with all the time points, all the structures
    ts_array = ts_cifti.get_fdata()

    # initialize a matrix representing 4D data (x, y, z, time_point)
    subcorticals_vol = np.zeros([bmf.volume_shape[0], bmf.volume_shape[1], bmf.volume_shape[2], ts_array.shape[0]])
    for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
        # get the values corresponding to the brain model
        bm_vals = ts_array[:, slc]

        # get the voxels/vertices corresponding to the current brain model
        ijk = bm.voxel
        # fill in data
        if (idx != 0) & (idx != 1): # indices 0 and 1 are cortical hemispheres
            # print(str(nam))
            subcorticals_vol[ijk[:, 0], ijk[:, 1], ijk[:, 2], :] = bm_vals.T 

    # save as nii 
    N = nb.Nifti1Image(subcorticals_vol,bmf.affine)

    Vs = nb.funcs.four_to_three(N)

    # ts_nifti = dir+'/sub-100307_ses-01_task-rest_space-subcortex_run-01_bold.nii'
    # nb.save(N,ts_nifti)

    return Vs


if __name__ == "__main__":
    explore_cifti()


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
# Script for importing the MDTB data set from super_cerebellum to general format.
import numpy as np
import nibabel as nb
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import surfAnalysisPy as surf
import nitools as nt

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'
hcp_dir = base_dir + '/HCP'

def explore_cifti():
    dir = hcp_dir + '/derivatives/100307/func'
    A=nb.load(dir+'/sub-100307_run-0_space-MSMSulc.dtseries.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)

    for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
        print(idx,str(nam),slc)
    pass

def explore_pconn():
    # B = nb.load('RSN-networks.32k_fs_LR.dlabel.nii')
    A = nb.load('./test_data/HCP_Yeo2011_17Networks.32k_fs_LR.pconn.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)

    for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
        print(idx,str(nam),slc)
    pass

def test_pdscalar():
    pass


def get_cereb_mask():
    dir2 = hcp_dir + '/derivatives/100307/func'
    A=nb.load(dir2+'/sub-100307_run-0_space-MSMSulc.dtseries.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)
    bmcl = bmf[bmf.name == 'CIFTI_STRUCTURE_CEREBELLUM_LEFT']
    bmcr = bmf[bmf.name == 'CIFTI_STRUCTURE_CEREBELLUM_RIGHT']
    bmc= bmcl + bmcr
    ijk = bmc.voxel
    X = np.zeros((bmc.volume_shape))
    X[ijk[:,0],ijk[:,1],ijk[:,2]] = 1
    N = nb.Nifti1Image(X,bmc.affine)
    nb.save(N,atlas_dir + '/tpl-MNI152NLin6AsymC' + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii')
    return bmc


def get_cortex():
    dir = hcp_dir + '/derivatives/100307/func'
    A=nb.load(dir+'/sub-100307_ses-01_task-rest_space-fsLR32k_run-01_bold.nii')
    ser = A.header.get_axis(0)
    bmf = A.header.get_axis(1)
    bmcl = bmf[bmf.name == 'CIFTI_STRUCTURE_CORTEX_LEFT']
    bmcr = bmf[bmf.name == 'CIFTI_STRUCTURE_CORTEX_RIGHT']

    maskl=np.zeros(32492,)
    maskl[bmcl.vertex]=1
    maskr=np.zeros(32492,)
    maskr[bmcr.vertex]=1
    colorM=np.array([[1,1,1,1],[0,0,1,1]])
    mask_gii = surf.map.make_label_gifti(data=maskl,label_names=['medwall','cortex'],label_RGBA=colorM)
    nb.save(mask_gii,atlas_dir + '/tpl-fs32k' + '/tpl-fs32k_hemi-L_mask.label.gii')
    mask_gii = surf.map.make_label_gifti(data=maskr,label_names=['medwall','cortex'],label_RGBA=colorM,anatomical_struct='CortexRight')
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


def reduce_cifti():
    wdir = '/Users/jdiedrichsen/Dropbox/projects/Pontine7T/'
    cimg = nb.load(wdir + 'beta_glm2_cereb_gray_S01.dscalar.nii')
    D = cimg.get_fdata()[0:2,:]
    bm = cimg.header.get_axis(1)
    row_axis = [f"row {r:03}" for r in range(D.shape[0])]
    row_axis = nb.cifti2.ScalarAxis(row_axis)
    header = nb.Cifti2Header.from_axes((row_axis, bm))
    cifti_img = nb.Cifti2Image(dataobj=D, header=header)
    nb.save(cifti_img,wdir + 'beta.dscalar.nii')
    pass

def cifti_to_nifti():
    wdir = '/Users/jdiedrichsen/Dropbox/projects/Pontine7T/'
    cimg = nb.load(wdir + 'beta.dscalar.nii')
    nimg = nt.volume_from_cifti(cimg)
    nb.save(nimg,wdir + 'beta.nii.gz')
    pass


def test_read_cifti():
    # Test a read-out of a cifti file in a different resolution than the atlas
    at,_ = am.get_atlas('MNISymC2')
    wdir = '/Users/jdiedrichsen/Dropbox/projects/Pontine7T/'
    data = at.read_data(wdir + 'beta.dscalar.nii',interpolation=1)
    X = at.data_to_nifti(data)
    nb.save(X,wdir + 'beta_res_2.nii.gz')
    pass



if __name__ == "__main__":
    # get_cereb_mask()
    # explore_cifti()
    test_read_cifti()
    # reduce_cifti()
    # sample_cifti()
    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
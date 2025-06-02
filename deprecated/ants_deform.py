def templateflow_xfm_h5_to_nii(filename):
    """This is a first attempt to understand the _xfm.h5 files from templateflow
    and convert them to nifti files. Not easy without documentation - so this is work
    in progress for now

    Args:
        filename (str): h5-filename
    """
    with h5py.File(filename, 'r') as f:
        P2 = np.array(f['TransformGroup']['2']['TranformParameters'])
        P2f = np.array(f['TransformGroup']['2']['TranformFixedParameters'])
        P1 = np.array(f['TransformGroup']['1']['TranformParameters'])
        P1f = np.array(f['TransformGroup']['1']['TranformFixedParameters'])

        # Wild guess on the paramaters for the second transform
        dim = P2f[:3].astype(int)
        trans = P2f[3:6]
        voxsize = P2f[6:9]
        rot = P2f[9:].reshape(3,3)
        affine = np.zeros((4,4))
        affine[:3,:3] = rot
        affine[:3,3] = trans
        affine[3,3] = 1
        P2 = P2.reshape(np.r_[dim,(3)])
        deform_img = nb.Nifti1Image(P2,affine)
        nb.save(deform_img, filename.replace('.h5','.nii'))
        pass

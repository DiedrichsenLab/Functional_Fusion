import numpy as np
import nibabel as nb
import pandas as pd


# download the atlas from FSL (MNI152NonLinear Asym 6 or sym) and make a mask for all the basal ganglia structures
# List of indices for Basal ganglia structures from HarvardOxford subcortical atlas"
# Left Caudate: 5
# Right Caudate: 16
# Left Pallidum: 7
# Right Pallidum: 18
# Left Putamen: 6
# Right Putamen: 17  
# Left Accumbens: 11
# Right Accumbens: 21

# list of indices for Thalamus structures from HarvardOxford subcortical atlas"
# Left Thalamus: 4
# Right Thalamus: 15

def make_mask_harvOx(atlas_file, mask_index = [5, 16, 7, 18, 6, 17, 11, 21]):
    """_uses thresholded max probability image created for HarvardOxford atlas from FSL to create a mask_

    Args:
        atlas_file (str): _full path (including name) to the atlas file_ Example name: "HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"
        mask_index (list, optional): _list of indices corresponding to bg structures from xml file_. 
                                    Defaults to [5, 16, 7, 18, 6, 17, 11, 21] for basal ganglia (gives one mask including left and right structures)
                                    Use [4, 15] for Thalamus (gives one mask including left and right structures).
    Returns:
        nb_mask (nb.nifti object): _nifti image of the mask_
    """
    
    # load the atlas
    nb_atl = nb.load(atlas_file)
    # get the data from the atlas and convert it to integer (it is saved as gloating point numbers)
    dat_atl = nb_atl.get_fdata().astype(int)
    # make a mask for thalamus only and test it
    dat_mask = np.isin(dat_atl, mask_index)*1 # multiply by 1 to convert boolean to integer

    # make a nifti imaging affine matrix and header from the original image of the atlas
    nb_mask = nb.Nifti1Image(dat_mask, affine=nb_atl.affine, header = nb_atl.header)

    return nb_mask

def make_thalamus_mask_harvOx(atlas_file, mask_index = [4, 15]):
    """_uses thresholded max probability image created for HarvardOxford atlas from FSL to create a mask_

    Args:
        atlas_file (str): _full path (including name) to the atlas file_ Example name: "HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"
        mask_index (list, optional): _list of indices corresponding to thalamus structures from xml file_. Defaults to [4, 15].
    Returns:
        th_mask (nifti image): _nifti image of the mask for the thalamus structures_
    """
    
    # load the atlas
    nb_atl = nb.load(atlas_file)
    # get the data from the atlas and convert it to integer (it is saved as gloating point numbers)
    dat_atl = nb_atl.get_fdata().astype(int)
    # make a mask for thalamus only and test it
    Xth = np.isin(dat_atl, mask_index)

    # make a nifti imaging affine matrix and header from the original image of the atlas
    th_mask = nb.Nifti1Image(Xth, affine=nb_atl.affine, header = nb_atl.header)
    # nb.save(new_image, '/Users/lshahsha/Desktop/th_left_test.nii')

    return th_mask


if __name__ == "__main__":
    path_to_atlas = "/Users/lshahsha/Desktop/HarvardOxford/HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"
    bg_mask = make_mask_harvOx(atlas_file = path_to_atlas, mask_index = [5, 16, 7, 18, 6, 17, 11, 21])
    th_mask = make_mask_harvOx(atlas_file = path_to_atlas, mask_index = [4, 15])

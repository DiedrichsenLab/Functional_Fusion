
import os
import Functional_Fusion.atlas_map as am

import numpy as np
import nibabel as nb

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

def save_bg_mask_maxprob(path_to_atlas, thr = 50, res = 1):
    """_code to make and save basal ganglia masks based on max probability images for HarvardOxford cortical/subcortical atlas_
    Args:
        path_to_scAtlas (str): _full path (including name) to the atlas file_ Example name: "HarvardOxford-sub-maxprob-thr50-1mm.nii.gz"
        thr (int, optional): _threshold for the max probability image_. Defaults to 50.
        res (int, optional): _resolution of the atlas (1 or 2 mm)_. Defaults to 1.
    """
    # directory to save the atlas mask
    atlas_dir = os.path.dirname(am.__file__) + '/Atlases'


    # this can be found in $FSLDIR/data/atlases
    path_to_nii = f"{path_to_atlas}/HarvardOxford-sub-maxprob-thr{thr}-{res}mm.nii.gz"
    
    # create mask for basal ganglia resolution 1 mm
    bg_mask = make_mask_harvOx(atlas_file = path_to_nii, mask_index = [5, 16, 7, 18, 6, 17, 11, 21])
    path_to_bg_atlas = f"{atlas_dir}/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-{res}_bgmask.nii"
    nb.save(bg_mask, path_to_bg_atlas)
    return


if __name__ == "__main__":

    # save masks for basal ganglia res 1 and 2 mm
    path_to_atlas = "/Users/lshahsha/Documents/myAtlases/HarvardOxford" # this can be found in $FSLDIR/data/atlases
    save_bg_mask_maxprob(path_to_atlas=path_to_atlas, res=1)
    save_bg_mask_maxprob(path_to_atlas=path_to_atlas, res=2)

    # testing functions from atlas map
    atlas, ainfo = am.get_atlas(atlas_str="MNIAsymBg2")

    print(atlas.structure)
    pass

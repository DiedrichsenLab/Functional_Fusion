# Test script for using atlases and atlasmaps for region definition and data extraction without the use of a Dataset object
import numpy as np
import nibabel as nb
import Functional_Fusion.atlas_map as am
import nitools as nt


# Data directory on the server for this example
data_dir = '/Volumes/diedrichsen_data$/data/Articulation/Articulotopy1'
glm_dir = data_dir + '/glm/glm_1'
surf_dir = data_dir + '/surfaceWB'


def test_define_region_surface():

    # Define the the region, get only left hemisphere
    atlas,_ = am.get_atlas('fs32k')
    atlas_left = atlas.get_hemisphere(0)
    # Equivalently you could have used
    # left,_ = am.get_atlas('fs32k_L')

    # Set the Gifti file for the region (func.gii or label.gii)
    roi_img = surf_dir + '/group/fs_LR.32k.L.ventralmotor.func.gii'
    subatlas = atlas_left.get_subatlas_image(roi_img)
    return subatlas

def test_atlas_map(subatlas):
    white = surf_dir + '/sub-01/sub-01.L.white.32k.surf.gii'
    pial = surf_dir + '/sub-01/sub-01.L.pial.32k.surf.gii'
    mask = glm_dir + '/sub-01/mask.nii'
    amap = am.AtlasMapSurf(subatlas.vertex[0],white,pial,mask)
    amap.build()
    return amap

def test_extract_data(amap):
    dnames = ['beta_0001.nii','beta_0002.nii','beta_0003.nii']
    fnames =[glm_dir + '/sub-01/' + d for d in dnames]
    data = amap.extract_data_native(fnames)
    g_data = amap.map_native_to_group(data)
    return data, g_data

if __name__ == "__main__":
    # get_cereb_mask()
    # explore_cifti()
    reg =  test_define_region_surface()
    amap = test_atlas_map(reg)
    data, g_data = test_extract_data(amap)
    pass
    # reduce_cifti()
    # sample_cifti()
    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
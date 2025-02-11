# Test script for different region definition and data extraction in Functional_Fusion.
import numpy as np
import nibabel as nb
import Functional_Fusion.atlas_map as am
import Functional_Fusion.region as reg
import nitools as nt


# Data directory on the server for this example 
data_dir = '/Volumes/diedrichsen_data$/data/Articulation/Articulotopy1'
glm_dir = data_dir + '/glm/glm_1'
surf_dir = data_dir + '/surfaceWB'


def test_define_region_surface(): 
    atlas,_ = am.get_atlas('fs32k')
    roi = surf_dir + '/fs_LR.32k.L.ventralmotor.func.gii'
    R1 = reg.RegionSurfaceImage(atlas,roi,hem=0,value=1)
    return R1 

if __name__ == "__main__":
    # get_cereb_mask()
    # explore_cifti()
    data =  test_define_region_surface()
    pass
    # reduce_cifti()
    # sample_cifti()
    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
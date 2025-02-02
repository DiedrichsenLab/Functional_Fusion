# Test for different forms of the atlas map
import pandas as pd
import shutil
from pathlib import Path
# import mat73
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import nibabel as nb
from Functional_Fusion.matrix import indicator
import sys
import os


def test_atlasmap_deform():
    at1,_ = am.get_atlas('MNISymC2')
    at2,_ = am.get_atlas('SUIT2')
    deform = am.default_atlas_dir + '/tpl-SUIT/tpl-SUIT_from-MNI152NLin2009cSymC_mode-image_xfm.nii'
    mask = am.default_atlas_dir + '/tpl-MNI152NLin2009cSymC/tpl-MNISymC_res-2_gmcmask.nii'
    amap = am.AtlasMapDeform(at2.world,deform,mask)
    amap.build(interpolation=1)

    []=am.get_data_nifti(fname,[amap])



if __name__ == "__main__":
    test_atlasmap_deform()
    # make_mdtb_suit()
    # test_atlas_sym()
    # test_atlas_deform()
    pass


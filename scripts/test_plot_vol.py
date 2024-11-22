import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as ut
import Functional_Fusion.plot as fplt
import nibabel as nb
import nilearn.plotting as nlp
import matplotlib.pyplot as plt
import nitools as nt

if __name__=="__main__":
    base_dir = ut.get_base_dir()
    # dn,ainf = am.get_atlas('MNISymDentate1')
    # dat,dinf,_ = ds.get_dataset(base_dir,'Language',atlas='MNISymDentate1', sess='ses-localizer',subj=[1],type='cond_fm_CondRun')
    adir = ut.default_atlas_dir
    bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_dentate.nii')
    # Project the functional data into the atlas space
    # fcn_img = dn.data_to_nifti(dat[0,:])

    c1 = [-25,-70,-43] # Lower left corner of image 
    c2 = [25,-40,-20] # Upper right corner of image

    pass

import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as ut
import Functional_Fusion.plot as fplt
import nibabel as nb
import nilearn.plotting as nlp
import matplotlib.pyplot as plt
import nitools as nt



def test_plot_roi(): 
    atlas_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Atlases'
    # Generate some random data in dentate atlas space
    dn,ainf = am.get_atlas('MNISymDentate1')
    data = np.random.normal(0,1,(3,dn.P))
    data = np.random.randint(1,33,(1,dn.P)).astype(np.int8) # Note that 0 is no assignment
    indx, color, labels = nt.read_lut(atlas_dir + '/tpl-MNI152NLin2009cSymC/atl-NettekovenSym32.lut')
    fplt.plot_dentate(data[0],cscale=[0,32],cmap=color)

if __name__=="__main__":
    test_plot_roi()
    pass
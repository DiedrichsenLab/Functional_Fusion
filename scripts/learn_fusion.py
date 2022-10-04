# Script for importing the MDTB data set from super_cerebellum to general format.
from time import gmtime
from pathlib import Path
import pandas as pd
import numpy as np
import atlas_map as am
from dataset import *
from scipy.linalg import block_diag
import nibabel as nb
import SUITPy as suit
import generativeMRF.full_model as fm
import generativeMRF.spatial as sp
import generativeMRF.arrangements as ar
import generativeMRF.emissions as em
import generativeMRF.evaluation as ev
import torch as pt
import Functional_Fusion.matrix as matrix
import matplotlib.pyplot as plt
import seaborn as sb
import sys

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

if sys.platform == "win32":
    base_dir = 'Y:\data\FunctionalFusion'
    sys.path.append('../')



if __name__ == "__main__":
    # mask = base_dir + '/Atlases/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    # suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)

    mdtb_dataset = DataSetMDTB(base_dir + '/MDTB')
    hcp_dataset = DataSetHcpResting(base_dir + '/HCP')
    pt7_dataset = DataSetPontine(base_dir + '/pontine7T')
    nishi_dataset = DataSetHcpResting(base_dir + '/Nishimoto')
    nn_dataset = DataSetNishi(base_dir + '/Nishimoto_103Task')
    fiel = ['run','task_name','reg_id','half']
    data_nn1,info_nn1 = nn_dataset.get_data('SUIT3','ses-01',
                                            'CondSes',fields=fiel)
    data_nn2,info_nn2 = nn_dataset.get_data('SUIT3','ses-02',
                                            'CondSes',fields=fiel)
    pass

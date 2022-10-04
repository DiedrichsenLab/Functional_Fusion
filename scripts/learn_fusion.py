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
    mask = base_dir + '/Atlases/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)

    mdtb_dataset = DataSetMDTB(base_dir + '/MDTB')
    hcp_dataset = DataSetHcpResting(base_dir + '/HCP')
    pt7_dataset = DataSetHcpResting(base_dir + '/pontine7T')
    nishi_dataset = DataSetHcpResting(base_dir + '/Nishimoto')

    X,D = mdtb_dataset.get_data('SUIT3','ses-s1','CondSes')
    r1 = reliability_within_subj(X,part_vec=D.half,cond_vec=D.cond_name)
    r2 = reliability_between_subj(X,cond_vec=D.cond_name)

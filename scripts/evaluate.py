# Script for evaluating DCBC on fused parcellation
from time import gmtime
from pathlib import Path
import pandas as pd
import numpy as np
import atlas_map as am
from dataset import *
from scipy.linalg import block_diag
import nibabel as nb
import SUITPy as suit
from learn_fusion import load_batch_fit
from learn_fusion import plot_parcel_flat
import torch as pt
import Functional_Fusion.matrix as matrix
import matplotlib.pyplot as plt
import seaborn as sb
import sys
import pickle

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

if __name__ == "__main__":
    fit_all([0,1,2])

    # load fit
    # fit_all([0, 1, 2])

    # 

    # mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    # sym_atlas = am.AtlasVolumeSymmetric('MNISymC3', mask_img=mask)
    # sess = [['ses-s1'],['ses-01'],['ses-01','ses-02']]
    # design_ind= ['cond_num_uni','task_id',',..']
    # info,models,Prop,V = load_batch_fit('SingleMDTB','MNISymC3',10)
    # parcel = pt.argmax(Prop,dim=1) # Get winner take all
    # parcel=parcel[:,sym_atlas.indx_reduced] # Put back into full space
    # plot_parcel_flat(parcel[0:3, :], sym_atlas, grid=[
    #                  1, 3], map_space='MNISymC')
    # pass
    # pass
    # Prop, V = fit_niter(data,design,K,n_iter)
    # r1 = ev.calc_consistency(Prop,dim_rem=0)
    # r2 = ev.calc_consistency(V[0],dim_rem=2)

    # parcel = pt.argmax(Prop,dim=1)
    # plot_parcel_flat(parcel,suit_atlas,(1,4))

    pass

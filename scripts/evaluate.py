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
import torch as pt
import Functional_Fusion.matrix as matrix
import matplotlib.pyplot as plt
import seaborn as sb
import sys
import pickle
import DCBC.DCBC_vol as dcbc
from learn_mdtb import get_all_mdtb

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'


def load_batch_fit(fname):
    wdir = base_dir + '/Models/'
    info = pd.read_csv(wdir + fname + '.tsv',sep='\t')
    with open(wdir + fname + '.pickle','rb') as file:
        models = pickle.load(file)
    n_iter = len(models)
        # Intialize data arrays
    Prop = pt.zeros((n_iter,models[0].arrange.K,models[0].arrange.P))
    V = []

    for i,M in enumerate(models):
        Prop[i,:,:] = M.arrange.logpi.softmax(axis=0)

        # Now switch the emission models accordingly:
        for j,em in enumerate(M.emissions):
            if i==0:
                V.append(pt.zeros((n_iter,em.M,info.K[i])))
            V[j][i,:,:]=em.V
    return info,models,Prop,V

def plot_parcel_flat(data,suit_atlas,grid,map_space='SUIT'):
    color_file = base_dir + '/Atlases/tpl-SUIT/atl-MDTB10.lut'
    color_info = pd.read_csv(color_file, sep=' ', header=None)
    MDTBcolors = np.zeros((11, 3))
    MDTBcolors[1:11, :] = color_info.iloc[:, 1:4].to_numpy()
    Nifti = suit_atlas.data_to_nifti(data)
    surf_data = suit.flatmap.vol_to_surf(Nifti, stats='mode',space=map_space)

    plt.figure
    for i in range(surf_data.shape[1]):
        plt.subplot(grid[0],grid[1],i+1)
        suit.flatmap.plot(surf_data[:,i], render='matplotlib',cmap=MDTBcolors, new_figure=False,overlay_type='label')

def plot_parcel_flat_best(model_names,grid):
    mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    atlas = am.AtlasVolumetric('MNISymC3',mask_img=mask)
    sym_atlas = am.AtlasVolumeSymmetric('MNISymC3',mask_img=mask)

    parcel=np.empty((len(model_names),atlas.P))

    for i,mn in enumerate(model_names):
        info,models,Prop,V = load_batch_fit(mn)
        j=np.argmax(info.loglik)
        par = pt.argmax(Prop[j,:,:],dim=0)+1 # Get winner take all 
        # If symmetric - project back to full map: 
        if mn[0:3]=='sym':
            par=par[sym_atlas.indx_reduced] # Put back into full space
        parcel[i,:]=par
    plot_parcel_flat(parcel,atlas,grid=grid,map_space='MNISymC') 
        

def eval_dcbc(parcels, suit_atlas, func_data=None, resolution=3):
    """DCBC: evaluate the resultant parcellation using DCBC
    Args:
        parcels (np.ndarray): the input parcellation, shape
        suit_atlas (<AtlasVolumetric>): the class object of atlas
        func_data (np.ndarray): the functional data,
                                shape (num_sub, N, P)
        resolution (np.float or int): the resolution of atlas in mm
    Returns:
        dcbc_values (np.ndarray): the DCBC values of subjects
    """
    dist = dcbc.compute_dist(suit_atlas.vox.T, resolution=resolution)
    if func_data is None:
        func_data, _, _ = get_sess_mdtb(atlas='SUIT3', ses_id='ses-s2')

    dcbc_values = []
    for sub in range(func_data.shape[0]):
        D = dcbc.compute_DCBC(parcellation=parcels,
                              dist=dist, func=func_data[sub].T)
        dcbc_values.append(D['DCBC'])

    return np.asarray(dcbc_values)

if __name__ == "__main__":
    model_name = ['asym_Md_space-MNISymC3_K-10',
                  'asym_Po_space-MNISymC3_K-10',
                  'asym_Ni_space-MNISymC3_K-10',
                  'asym_MdPoNi_space-MNISymC3_K-10']
    
    plot_parcel_flat_best(model_name,[2,2])

    # Evaluate DCBC on left out dataset
    mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    atlas = am.AtlasVolumetric('MNISymC3',mask_img=mask)

    sess = [['ses-s1'],['ses-01'],['ses-01','ses-02']]
    design_ind= ['cond_num_uni','task_id',',..']
    info,models,Prop,V = load_batch_fit('asym_Md','MNISymC3',10)
    parcel = pt.argmax(Prop,dim=1)+1 # Get winner take all 
    parcel=parcel[:,atlas.indx_reduced] # Put back into full space


    # Evaluate case
    # T, gbase, lb, parcellation = learn_half(
    #     K=10, e='VMF', runs=np.arange(1, 17))
    # T.to_csv('coserrs_wVMF.csv')
    # plot_parcel_flat(parcellation, suit_atlas, grid=[
                    #  1, 1], save_nii=False)  # Plot flat map
    data_eval, _, _ = get_all_mdtb(atlas='MNISymC3')
    dcbc_values = eval_dcbc(parcellation.numpy(), atlas,
                            func_data=data_eval, resolution=3)


    # parcel = pt.argmax(Prop,dim=1)
    # plot_parcel_flat(parcel,suit_atlas,(1,4))


    pass

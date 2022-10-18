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

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'


def load_batch_fit(name,atl_name,K):
    wdir = base_dir + '/Models'
    fname = f'/{name}_space-{atl_name}_K-{K}'
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
                V.append(pt.zeros((n_iter,em.M,K)))
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

if __name__ == "__main__":


    mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    atlas = am.AtlasVolumetric('MNISymC3',mask_img=mask)

    #sess = [['ses-s1'],['ses-01'],['ses-01','ses-02']]
    #design_ind= ['cond_num_uni','task_id',',..']
    info,models,Prop,V = load_batch_fit('asym_Md','MNISymC3',10)
    parcel = pt.argmax(Prop,dim=1)+1 # Get winner take all 
    # parcel=parcel[:,sym_atlas.indx_reduced] # Put back into full space
    plot_parcel_flat(parcel[0:3,:],atlas,grid=[2,3],map_space='MNISymC') 
    # pass
    # pass
    # Prop, V = fit_niter(data,design,K,n_iter)
    # r1 = ev.calc_consistency(Prop,dim_rem=0)
    # r2 = ev.calc_consistency(V[0],dim_rem=2)


    # parcel = pt.argmax(Prop,dim=1)
    # plot_parcel_flat(parcel,suit_atlas,(1,4))


    pass

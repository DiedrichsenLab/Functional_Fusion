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
from learn_mdtb import get_mdtb_parcel
import matplotlib.pyplot as plt
import seaborn as sb
import sys

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

if sys.platform == "win32":
    base_dir = 'Y:\data\FunctionalFusion'
    sys.path.append('../')

def get_all_mdtb(atlas='SUIT3'):
    mdtb_dataset = DataSetMDTB(base_dir + '/MDTB')
    fiel = ['study','half','common','cond_name','cond_num','cond_num_uni','common']
    data_mdtb1,info_mdtb1 = mdtb_dataset.get_data(atlas,'ses-s1',
                                                  'CondHalf',fields=fiel)
    data_mdtb2,info_mdtb2 = mdtb_dataset.get_data(atlas,'ses-s2',
                                                  'CondHalf',fields=fiel)
    info_mdtb1['sess']=np.ones((info_mdtb1.shape[0],))*1
    info_mdtb2['sess']=np.ones((info_mdtb2.shape[0],))*2
    info_mdtb = pd.concat([info_mdtb1,info_mdtb2],ignore_index=True,sort=False)
    data_mdtb=np.concatenate([data_mdtb1,data_mdtb2],axis=1)
    return data_mdtb, info_mdtb, mdtb_dataset

def get_all_pontine(atlas='SUIT3'):
    pt7_dataset = DataSetPontine(base_dir + '/pontine7T')
    fiel = ['task_name','task_num','half']
    data_pt,info_pt = pt7_dataset.get_data(atlas,'ses-01',
                                           'TaskHalf',fields=fiel)
    return data_pt, info_pt, pt7_dataset

def get_all_nishi(atlas='SUIT3'):
    nn_dataset = DataSetNishi(base_dir + '/Nishimoto_103Task')
    fiel = ['task_name','reg_id','half']
    data_nn1,info_nn1 = nn_dataset.get_data(atlas,'ses-01',
                                            'CondHalf',fields=fiel)
    data_nn2,info_nn2 = nn_dataset.get_data(atlas,'ses-02',
                                            'CondHalf',fields=fiel)
    info_nn1['sess']=np.ones((info_nn1.shape[0],))*1
    info_nn2['sess']=np.ones((info_nn2.shape[0],))*2
    info_nn = pd.concat([info_nn1,info_nn2],ignore_index=True,sort=False)
    data_nn=np.concatenate([data_nn1,data_nn2],axis=1)
    return data_nn, info_nn, nn_dataset

def fit_fusion(data,design,K=10,intialization='random'):
    
    P = data[0].shape[2]
    ar_model = ar.ArrangeIndependent(K=K, P=P, spatial_specific=True,
                                         remove_redundancy=False)
    pt.normal(0.0,1.0,ar_model.logpi.shape,out=ar_model.logpi)

    # Initialize emission models 
    em_models=[]
    for i,ds in enumerate(data):
        em_model = em.MixVMF(K=K, N=40, P=P, X=design[i], uniform_kappa=True)
        em_model.initialize(ds)
        em_models.append(em_model)
 
    # Use random prior
    M = fm.FullMultiModel(ar_model, em_models)

    # Step 5: Estimate the parameter thetas to fit the new model using EM
    M, ll, theta, U_hat = M.fit_em(Y=data, iter=20, tol=0.00001, fit_arrangement=True)
    return M,ll,theta,U_hat

def fit_niter(data,design,K,n_iter):
    # Initialize container
    Prop = pt.zeros((n_iter,K,n_vox))
    V = []
    for s in range(n_sets):
        V.append(pt.zeros((n_iter,design[s].shape[1],K)))

    for i in range(n_iter):
        print(f'iter: {i}')
        M,ll,theta,U_hat = fit_fusion(data,design,K=K)
        pp = M.arrange.logpi.softmax(axis=0)
        if i == 0: 
            indx = np.arange(K)
        else:
            indx = ev.matching_greedy(Prop[0,:,:],pp)
        Prop[i,:,:]=pp[indx,:]
        for s in range(n_sets):
            V[s][i,:,:]=M.emissions[s].V[:,indx]
    return Prop, V

def plot_parcel_flat(data,suit_atlas,grid):
    color_file = base_dir + '/Atlases/tpl-SUIT/atl-MDTB10.lut'
    color_info = pd.read_csv(color_file, sep=' ', header=None)
    MDTBcolors = np.zeros((11, 3))
    MDTBcolors[1:11, :] = color_info.iloc[:, 1:4].to_numpy()
    Nifti = suit_atlas.data_to_nifti(data)
    surf_data = suit.flatmap.vol_to_surf(Nifti, stats='mode')

    plt.figure
    for i in range(surf_data.shape[1]):
        plt.subplot(grid[0],grid[1],i+1)
        suit.flatmap.plot(surf_data[:,i], render='matplotlib',cmap=MDTBcolors, new_figure=False)


if __name__ == "__main__":
    mask = base_dir + '/Atlases/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    data_mdtb, info_mdtb, mdtb_dataset = get_all_mdtb(atlas='SUIT3')
    # data_pt, info_pt, pt7_dataset = get_all_pontine(atlas='SUIT3')
    # data_nn, info_nn, nn_dataset = get_all_nishi(atlas='SUIT3')
    # design matrix for emission model 
    X_mdtb = matrix.indicator(info_mdtb.cond_num_uni)
    n_vox = data_mdtb.shape[2]
    n_iter = 4
    data = [data_mdtb]
    design = [X_mdtb]
    n_sets = len(data)
    K=10  # Number of parcels 

    Prop, V = fit_niter(data,design,K,n_iter)
    parcel = pt.argmax(Prop,dim=1)
    plot_parcel_flat(parcel,suit_atlas,(1,4))


    pass

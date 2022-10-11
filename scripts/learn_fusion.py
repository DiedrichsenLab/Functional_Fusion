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
import pickle

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

if sys.platform == "win32":
    base_dir = 'Y:\data\FunctionalFusion'
    sys.path.append('../')

def get_all_any(dataset,atlas='SUIT3',sess='all',type='CondHalf'):
    if dataset == 'MDTB':
        data,info,dataset = get_all_mdtb(atlas,sess,type)
    if dataset == 'pontine7T':
        data,info,dataset = get_all_mdtb(atlas,sess,type)
    if dataset == 'nishimotor':
        data,info,dataset = get_all_mdtb(atlas,sess,type)
    return data,info,dataset

def get_all_mdtb(atlas='SUIT3',sess='all',type='CondHalf'):
    mdtb_dataset = DataSetMDTB(base_dir + '/MDTB')
    fiel = ['study','half','common','cond_name','cond_num','cond_num_uni','common']
    info_mdtb = []
    data_mdtb = []
    if sess=='all':
        sess=['ses-s1','ses-s2']

    for s in sess:
        dat,info = mdtb_dataset.get_data(atlas,s,type,fields=fiel)
        data_mdtb.append(dat)
        info['sess']=[s]*info.shape[0]
        info_mdtb.append(info)
    info_mdtb = pd.concat(info_mdtb,ignore_index=True,sort=False)
    data_mdtb=np.concatenate(data_mdtb,axis=1)
    return data_mdtb, info_mdtb, mdtb_dataset

def get_all_pontine(atlas='SUIT3',type='TaskHalf'):
    pt7_dataset = DataSetPontine(base_dir + '/pontine7T')
    fiel = ['task_name','task_num','half']
    data_pt,info_pt = pt7_dataset.get_data(atlas,'ses-01',
                                           'TaskHalf',fields=fiel)
    return data_pt, info_pt, pt7_dataset

def get_all_nishi(atlas='SUIT3',sess='all',type='CondHalf'):
    nn_dataset = DataSetNishi(base_dir + '/Nishimoto_103Task')
    fiel = ['task_name','reg_id','half']
    info_nn = []
    data_nn = []
    if sess=='all':
        sess=['ses-01','ses-02']

    for s in sess:
        dat,info = nn_dataset.get_data(atlas,s,type,fields=fiel)
        data_nn.append(dat)
        info['sess']=[s]*info.shape[0]
        info_nn.append(info)

    info_nn = pd.concat(info_nn,ignore_index=True,sort=False)
    data_nn=np.concatenate(data_nn,axis=1)
    return data_nn, info_nn, nn_dataset

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

def align_fits(models,inplace=True):
    """Aligns the prior probabilities and emission models
    across different model fits, returns parameters in aligned form
    if Inplace==True, it also aignes the model parameters in the models themselves.

    Args:
        models (list): List of full models
        inplace (bool): If true (default), it al
    Returns:
        Prop: Prior Probabilites as 3d-arrays (aligned)
        V: List of Mean vectors as 3d-arrays (aligned)
    """
    n_iter = len(models)
    K = models[0].arrange.K
    K_em = models[0].emissions[0].K
    n_vox = models[0].arrange.P

    # Intialize data arrays
    Prop = pt.zeros((n_iter,K,n_vox))
    V = []

    for i,M in enumerate(models):

        pp = M.arrange.logpi.softmax(axis=0)
        if i == 0:
            indx = np.arange(K)
        else:
            indx = ev.matching_greedy(Prop[0,:,:],pp)
        Prop[i,:,:]=pp[indx,:]
        if inplace:
            models[i].arrange.logpi=models[i].arrange.logpi[indx,:]

        # Now switch the emission models accordingly:
        for j,em in enumerate(M.emissions):
            if i==0:
                V.append(pt.zeros((n_iter,em.M,K_em)))
            if K == K_em: # non-symmetric model
                V[j][i,:,:]=em.V[:,indx]
            else:
                V[j][i,:,:]=em.V[:,np.concatenate([indx,indx+K])]

            if inplace:
                em.V=V[j][j,:,:]
    return Prop, V


def batch_fit(datasets,sess,design_ind,subj=None,
                atlas=None,K=10,arrange='independent',emission='VMF',
                n_iter=10,save=True,name=None):
    """ Executes a set of fits starting from random starting values
    saves the

    Args:
        datasets (list): _description_
        sess (list): _description_
        design_ind (list): _description_
        subj (list, optional): _description_. Defaults to None.
        atlas (Atlas): Atlas to be used. Defaults to None.
        K (int): Number of parcels. Defaults to 10.
        arrange (str): Type of arangement model. Defaults to 'independent'.
        emission (list / strs): Type of emission models. Defaults to 'VMF'.
        n_iter (int, optional): Number of fits random starting values. Defaults to 10.
        save (bool): Save the resulting fits? Defaults to True.
        name (str): Name of model (for filename). Defaults to None.

    Returns:
        info (pd.DataFrame):
    """
    # Load all necessary data and designs
    n_sets = len(datasets)
    data = []
    design = []
    if sess is None:
        sess = ['all']*n_sets

    for i in range(n_sets):
        dat,info,ds = get_all_any(datasets[i],atlas=atlas.name,sess=sess[i])
        if subj is None:
            data.append(dat)
        else:
            data.append(dat[subj[i],:,:])
        X = matrix.indicator(info[design_ind].values.reshape(-1,))
        design.append(X)

    # Collect info and fits and iterate
    models=[]
    info = pd.DataFrame({'name':[name]*n_iter,
                         'atlas':[atlas.name]*n_iter,
                         'K':[K]*n_iter,
                         'datasets':[datasets]*n_iter,
                         'sess':[sess]*n_iter,
                         'subj':[subj]*n_iter,
                         'arrange':[arrange]*n_iter,
                         'emission':[emission]*n_iter});

    # Check for size of Atlas + whether symmetric
    if isinstance(atlas,am.AtlasVolumeSymmetric):
        P_arrange = atlas.Psym
        K_arrange = np.ceil(K/2).astype(int)
    else:
        P_arrange = atlas.P
        K_arrange = K

    for i in range(n_iter):
        print(f'iter: {i}')

        # Initialize arrangement model
        if arrange=='independent':
            ar_model = ar.ArrangeIndependent(K=K_arrange, P=P_arrange,
                                             spatial_specific=True,
                                             remove_redundancy=False)
        else:
            raise(NameError(f'unknown arrangement model:{arrange}'))
        # Intialize randomly
        pt.normal(0.0,1.0,ar_model.logpi.shape,out=ar_model.logpi)

        # Initialize emission models
        em_models=[]
        for j,ds in enumerate(data):
            if emission=='VMF':
                em_model = em.MixVMF(K=K, N=40, P=atlas.P,
                                     X=design[j], uniform_kappa=True)
            else:
                raise((NameError(f'unknown emission model:{emission}')))
            em_model.initialize(ds)
            em_models.append(em_model)

        # Make a full fusion model
        if isinstance(atlas,am.AtlasVolumeSymmetric):
                M = fm.FullMultiModelSymmetric(ar_model, em_models,
                                               atlas.indx_full,atlas.indx_reduced,
                                               same_parcels=False)
        else:
                M = fm.FullMultiModel(ar_model, em_models)

        # Step 5: Estimate the parameter thetas to fit the new model using EM
        M, ll, theta, U_hat = M.fit_em(Y=data, iter=20,
                        tol=0.00001, fit_arrangement=True)
        models.append(M)

    # Align the different models
    align_fits(models)

    # Save the fits and information
    if save is True:
        wdir = base_dir + '/Models'
        fname = f'/{name}_space-{atlas.name}_K-{K}'
        info.to_csv(wdir + fname + '.tsv',sep='\t')
        with open(wdir + fname + '.pickle','wb') as file:
            pickle.dump(models,file)
    return info,models

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

if __name__ == "__main__":
    # mask = base_dir + '/Atlases/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    # suit_atlas = am.AtlasVolumetric('SUIT3',mask_img=mask)

    mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    sym_atlas = am.AtlasVolumeSymmetric('MNISymC3',mask_img=mask)

    #datasets = ['MDTB']
    #sess = [['ses-s1']]
    #design_ind= ['cond_num_uni']
    #batch_fit(datasets,sess,design_ind,atlas=sym_atlas,
    #           K=10,name='SingleMDTB',n_iter=10, save=True)
    info,models,Prop,V = load_batch_fit('SingleMDTB','MNISymC3',10)
    parcel = pt.argmax(Prop,dim=1) # Get winner take all 
    parcel=parcel[:,sym_atlas.indx_reduced] # Put back into full space
    plot_parcel_flat(parcel[0:3,:],sym_atlas,grid=[1,4]) 
    pass
    pass
    # Prop, V = fit_niter(data,design,K,n_iter)
    # r1 = ev.calc_consistency(Prop,dim_rem=0)
    # r2 = ev.calc_consistency(V[0],dim_rem=2)


    # parcel = pt.argmax(Prop,dim=1)
    # plot_parcel_flat(parcel,suit_atlas,(1,4))


    pass

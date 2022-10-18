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


def batch_fit(datasets,sess,type=None,design_ind=None,part_ind=None,subj=None,
                atlas=None,K=10,arrange='independent',emission='VMF',
                n_starts=10,n_iter=20,save=True,name=None):
    """ Executes a set of fits starting from random starting values
    saves the

    Args:
        datasets (list): List of dataset names to be used as training
        sess (list): List of list of sessions to be used for each
        type (list): List the type  
        design_ind (list): Name of the info-field that indicates the condition
        part_ind (list): Name of the field indicating independent partitions of the data
        subj (list, optional): _description_. Defaults to None
        atlas (Atlas): Atlas to be used. Defaults to None.
        K (int): Number of parcels. Defaults to 10.
        arrange (str): Type of arangement model. Defaults to 'independent'.
        emission (list / strs): Type of emission models. Defaults to 'VMF'.
        n_starts (int): Number of random starting values. default: 10
        n_iter (int): Maximal number of iterations per fit: default: 20
        save (bool): Save the resulting fits? Defaults to True.
        name (str): Name of model (for filename). Defaults to None.

    Returns:
        info (pd.DataFrame):
    """
    # Load all necessary data and designs
    n_sets = len(datasets)
    data = []
    design = []
    part_vec = [] 

    # Set defaults for data sets:
    if sess is None:
        sess = ['all'] * n_sets
    if part_ind is None: 
        part_ind = [None] * n_sets
    if type is None: 
        type = [None] * n_sets

    
    # Run over datasets get data + design 
    for i in range(n_sets):
        dat,info,ds = get_dataset(base_dir,datasets[i],atlas=atlas.name,sess=sess[i],type=type[i])
        if subj is None:
            data.append(dat)
        else:
            data.append(dat[subj[i],:,:])
        X = matrix.indicator(info[design_ind[i]].values.reshape(-1,))
        design.append(X)
        if part_ind[i] is None:
            part_vec.append(None)
        else:
            part_vec.append(info[part_ind[i]].values)

    # Initialize data frame for results
    models=[]
    info = pd.DataFrame({'name':[name]*n_starts,
                         'atlas':[atlas.name]*n_starts,
                         'K':[K]*n_starts,
                         'datasets':[datasets]*n_starts,
                         'sess':[sess]*n_starts,
                         'type':[type]*n_starts,
                         'subj':[subj]*n_starts,
                         'arrange':[arrange]*n_starts,
                         'emission':[emission]*n_starts,
                         'loglik':[np.nan]*n_starts});

    # Check for size of Atlas + whether symmetric
    if isinstance(atlas,am.AtlasVolumeSymmetric):
        P_arrange = atlas.Psym
        K_arrange = np.ceil(K/2).astype(int)
    else:
        P_arrange = atlas.P
        K_arrange = K

    # Iterate over the number of fits
    ll = np.empty((n_starts,n_iter))
    for i in range(n_starts):
        print(f'start: {i}')

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
                em_model = em.MixVMF(K=K, N=40, 
                                     P=atlas.P,
                                     X=design[j], 
                                     part_vec=part_vec[j],
                                     uniform_kappa=True)
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
        M, ll[i,:], theta, U_hat = M.fit_em(Y=data, iter=n_iter,
                        tol=0.00001, fit_arrangement=True)
        info.loglik.iloc[i] = ll[i,~np.isnan(ll[i,:])][-1]
        M.clear()
        models.append(M)

    # Align the different models
    Prop,V = align_fits(models)

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

def fit_all(set_ind=[0,1,2]):
    # Data sets need to numpy arrays to allow indixing by list
    datasets = np.array(['Mdtb','Pontine','Nishimoto'],
                    dtype = object)
    sess = np.array([['ses-s1','ses-s2'],
            ['ses-01'],
            ['ses-01','ses-02']],
            dtype = object)
    type = np.array(['CondHalf','TaskHalf','CondHalf'],
            dtype = object)
    design_ind= np.array(['cond_num_uni','task_num','reg_id'],
            dtype = object)
    part_ind = np.array(['half','half','half'],
            dtype = object)

    # Use specific mask / atlas. 
    mask = base_dir + '/Atlases/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    atlas = [am.AtlasVolumetric('MNISymC3',mask_img=mask),
             am.AtlasVolumeSymmetric('MNISymC3',mask_img=mask)]

    # Give a overall name for the type of model
    mname =['asym','sym']
    
    #Generate a dataname from first two letters of each training data set 
    dataname = [datasets[i][0:2] for i in set_ind]
    
    for i in [0]:
        name = mname[i] + '_' + ''.join(dataname) 
        batch_fit(datasets[set_ind],
              sess = sess[set_ind],
              type = type[set_ind],
              design_ind = design_ind[set_ind],
              part_ind = part_ind[set_ind],
              atlas=atlas[i],
              K=10,name=name,n_starts=10, save=True)

if __name__ == "__main__":
    # fit_all([0])
    # fit_all([1])
    fit_all([0])
    fit_all([0,1,2])

    #sess = [['ses-s1'],['ses-01'],['ses-01','ses-02']]
    #design_ind= ['cond_num_uni','task_id',',..']
    # info,models,Prop,V = load_batch_fit('SingleMDTB','MNISymC3',10)
    # parcel = pt.argmax(Prop,dim=1) # Get winner take all 
    # parcel=parcel[:,sym_atlas.indx_reduced] # Put back into full space
    # plot_parcel_flat(parcel[0:3,:],sym_atlas,grid=[1,3],map_space='MNISymC') 
    # pass
    # pass
    # Prop, V = fit_niter(data,design,K,n_iter)
    # r1 = ev.calc_consistency(Prop,dim_rem=0)
    # r2 = ev.calc_consistency(V[0],dim_rem=2)


    # parcel = pt.argmax(Prop,dim=1)
    # plot_parcel_flat(parcel,suit_atlas,(1,4))


    pass

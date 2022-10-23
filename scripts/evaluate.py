# Script for evaluating DCBC on fused parcellation
from time import gmtime
from pathlib import Path
import pandas as pd
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.matrix as matrix
from Functional_Fusion.dataset import *
import generativeMRF.emissions as em 
import generativeMRF.arrangements as ar 
import generativeMRF.full_model as fm
import generativeMRF.evaluation as ev

from scipy.linalg import block_diag
import nibabel as nb
import SUITPy as suit
import torch as pt
import matplotlib.pyplot as plt
import seaborn as sb
import sys
import pickle
import DCBC.DCBC_vol as dcbc
from learn_mdtb import get_sess_mdtb

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:\data\FunctionalFusion'
if not Path(base_dir).exists():
    raise(NameError('Could not find base_dir'))

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
    """Load a bunch of model fits, selects the best from 
    each of them and plots the flatmap of the parcellation
    ToDo: Align colors across different parcellations- 
    Pick good color schemes for different K + symmetric parcels  
    """
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

def cross_prediction_error(M,tdata,U_hats):
    """Evaluates the predictions from a trained full model on some testdata. 
    The full model consists of a trained arrangement model and is 
    combined with untrained emission model for the test data. 
    The function then trains the emission model based on N-1 subjects
    (arrangement fixed) and evaluates the left-out subjects using different
    parcellations (Uhats)

    Args:
        M (full model): Full model including emission model for test data
        tdata (ndarray): (numsubj x N x P) array of test data
        U_hats (list): List of strings or tensors, the indi. parcellation
            'group':   Group-parcellation 
            'floor':   Noise-floor (E-step on left-out subject)
            pt.tensor: Any arbitrary Individual parcellation based on outside data
    Returns: 
        A num_eval x num_subj matrix of cosine errors
    """
    num_subj = tdata.shape[0]
    subj = np.arange(num_subj)
    group_parc = M.arrange.marginal_prob()
    pred_err = np.empty((len(U_hats),num_subj))
    for s in range(num_subj):
        print(f'Subject:{s}')
        # initialize the emssion model using all but one subject 
        M.emissions[0].initialize(tdata[subj != s,:,:])
        # For fitting an emission model witout the arrangement model,
        # We can not without multiple starting values 
        M,ll,theta,Uhat = M.fit_em(
                iter=200, tol=0.1,
                fit_emission=True, 
                fit_arrangement=False,
                first_evidence=False)
        X = M.emissions[0].X
        dat = pt.linalg.pinv(X) @ tdata[subj==s,:,:]
        for i,crit in enumerate(U_hats):
            if crit=='group':
                U = group_parc
            elif crit=='floor':
                U,ll = M.Estep(Y=pt.tensor(tdata[subj==s,:,:]).unsqueeze(0))
            else:
                U = crit[subj==s,:,:]
            a=ev.coserr(dat, M.emissions[0].V, U,
                        adjusted=True, soft_assign=True)
            pred_err[i,s] = a
    return pred_err

def calc_prediction_error(model_names,test_data,test_sess,
                    design_ind,part_ind=None,
                    eval_types=['group','floor'],
                    indivtrain_ind=None,indivtrain_values=[0]):
    """ Calculates a prediction error using a test_data set 
    and test_sess. 
    if indivtrain_ind is given, it splits the test_data set
    again and uses one half to derive an individual parcellation 
    (using the model) and the other half to evaluate it. 
    The Means of the parcels are always estimated on N-1 subjects 
    and evaluated on the Nth left-out subject   

    Args:
        model_names (list or str): Name of model fit (tsv/pickle file)
        test_data (str): Name of test data set 
        test_sess (list): List or sessions to include into test_data 
        design_ind (str): Fieldname of the condition vector in test-data info
        part_ind (str): Fieldname of partition vector in test-data info
        eval_types (list): Defaults to ['group','floor'].
        indivtrain_ind (str): If given, data will be split for individual
             training along this field in test-data info. Defaults to None.
        indivtrain_values (list): Values of field above to be taken as
             individual training sets.  

    Returns:
        data-frame with model evalution
    """
    wdir = base_dir + '/Models/'
    tdata,tinfo,tds = get_dataset(base_dir,test_data,
                              atlas='MNISymC3',sess=test_sess)
    # For testing: tdata=tdata[0:5,:,:]
    num_subj = tdata.shape[0]
    results = pd.DataFrame()
    if not isinstance(model_names,list):
        model_names = [model_names]
    
    # Get condition and partition vector of test data 
    cond_vec = tinfo[design_ind].values.reshape(-1,)
    if part_ind is None:
        part_vec = np.zeros((tinfo.shape[0],),dtype=int)
    else:
        part_vec = tinfo[part_ind].values

    # Decide how many splits we need 
    if indivtrain_ind is None:
        n_splits = 1
    else:
        n_splits = len(indivtrain_values)

    # Now loop over possible models we want to evaluate 
    for model_name in model_names:
        minfo = pd.read_csv(wdir + model_name + '.tsv',sep='\t')
        with open(wdir + model_name + '.pickle','rb') as file:
            models = pickle.load(file)
        n_iter = len(models)

        for i,m in enumerate(models[0:2]):
            # Loop over the splits - if split then train a individual model
            for n in range(n_splits):
                # ------------------------------------------
                # Train an emission model on the individual training data and get a Uhat (individual parcellation) from it.
                if indivtrain_ind is not None:
                    train_indx = tinfo[indivtrain_ind]==indivtrain_values[n]
                    test_indx = tinfo[indivtrain_ind]!=indivtrain_values[n]
                    indivtrain_em = em.MixVMF(K=minfo.K[i], N=40, 
                             P = m.emissions[0].P,
                             X = matrix.indicator(cond_vec[train_indx]), 
                             part_vec=part_vec[train_indx],
                             uniform_kappa=True)
                    indivtrain_em.initialize(tdata[:,train_indx,:])
                    m.emissions = [indivtrain_em]
                    m.nparams = m.arrange.nparams + indivtrain_em.nparams 
                    m,ll,theta,U_indiv = m.fit_em(
                        iter=200, tol=0.1,
                        fit_emission=True, 
                        fit_arrangement=False,
                        first_evidence=False)
                    all_eval = eval_types + [U_indiv]
                else:
                    test_indx =  np.ones((tinfo.shape[0],),dtype=bool)
                    all_eval = eval_types
                # ------------------------------------------
                # Now build the model for the test data and crossvalidate
                # across subjects
                em_model = em.MixVMF(K=minfo.K[i], N=40, 
                             P=m.emissions[0].P,
                             X=matrix.indicator(cond_vec[test_indx]), 
                             part_vec=part_vec[test_indx],
                             uniform_kappa=True)
                # Add this single emission model
                m.emissions = [em_model] 
                # recalculate total number parameters
                m.nparams = m.arrange.nparams + em_model.nparams 
                # To CARO: You could copy the function and then replace this prediction_error function with the DCBC claculation for group and individual parcellation:  
                res = cross_prediction_error(m,tdata[:,test_indx,:],all_eval)
                # ------------------------------------------
                # Collect the information from the evaluation 
                # in a data frame
                ev_df = pd.DataFrame({'model_name':[minfo.name[i]]*num_subj,
                                'atlas':[minfo.atlas[i]]*num_subj,
                                'K':[minfo.K[i]]*num_subj,
                                'model_num':[i]*num_subj,
                                'train_data':[minfo.datasets[i]]*num_subj,
                                'train_loglik':[minfo.loglik[i]]*num_subj,
                                'test_data':[test_data]*num_subj,
                                'indivtrain_ind':[indivtrain_ind]*num_subj,
                                'indivtrain_val':[indivtrain_values[n]]*num_subj,
                                'subj_num':np.arange(num_subj)})
                # Add all the evaluations to the data frame
                for e,ev in enumerate(all_eval):
                    if isinstance(ev,str):
                        ev_df['coserr_' + ev]=res[e,:]
                    else:
                        ev_df[f'coserr_ind{e}']=res[e,:]
                results = pd.concat([results, ev_df], ignore_index=True)
    return results



def calc_dcbc():
    """ Calculates dcbc using a test_data set 
    and test_sess. 
    """
    pass


def eval_dcbc_group(model_names, space, testdata):
    s = space.split('C')[0]
    r = int(space.split('C')[1])
    mask = base_dir + \
        f'/Atlases/tpl-MNI152NLIn2000cSymC/tpl-{s}C_res-{r}_gmcmask.nii'
    atlas = am.AtlasVolumetric(space, mask_img=mask)

    if testdata is None:
        tdata, _, _ = get_sess_mdtb(atlas=space, ses_id='ses-s2')
    elif testdata == 'Md':
        tdata, tinfo, _ = get_dataset(base_dir, 'MDTB',
                                      atlas=space, type='CondHalf')
    elif testdata == 'Po':
        tdata, tinfo, _ = get_dataset(base_dir, 'Pontine',
                                      atlas=space, type='TaskHalf')
    elif testdata == 'Ni':
        tdata, tinfo, _ = get_dataset(base_dir, 'Nishimoto',
                                      atlas=space, type='CondHalf')
    
    wdir = base_dir + '/Models/'
    

    # parcel = np.empty((len(model_names), atlas.P))
    results = pd.DataFrame()
    for i, mn in enumerate(model_names):
        
        info, models, Prop, V = load_batch_fit(mn)
        j = np.argmax(info.loglik)
        par = pt.argmax(Prop[j, :, :], dim=0) + 1  # Get winner take all
        # parcel[i, :] = par
        if i == 0:
            dcbc = np.zeros((len(model_names), tdata.shape[0]))
        dcbc[i, :] = eval_dcbc(par, tdata, atlas, r)
        minfo = pd.read_csv(wdir + mn + '.tsv', sep='\t')
        num_subj = tdata.shape[0]

        ev_df = pd.DataFrame({'model_name': [minfo.name[i]] * num_subj,
                            'atlas': [minfo.atlas[i]] * num_subj,
                            'K': [minfo.K[i]] * num_subj,
                            'model_num': [i] * num_subj,
                            'train_data': [minfo.datasets[i]] * num_subj,
                            'train_loglik': [minfo.loglik[i]] * num_subj,
                            'test_data': [testdata] * num_subj,
                            'subj_num': np.arange(num_subj),
                            'dcbc': dcbc[i,:]
                            })
        
        
        results = pd.concat([results, ev_df], ignore_index=True)
    
    return results


def eval_dcbc(parcels, testdata, atlas, resolution=3, trim_nan=False):
    """DCBC: evaluate the resultant parcellation using DCBC
    Args:
        parcels (np.ndarray): the input parcellation, shape
        atlas (<AtlasVolumetric>): the class object of atlas
        testdata (np.ndarray): the functional test dataset,
                                shape (num_sub, N, P)
        resolution (np.float or int): the resolution of atlas in mm
        trim_nan (boolean): if true, make the nan voxel label will be
                            removed from DCBC calculation. Otherwise,
                            we treat nan voxels are in the same parcel
                            which is label 0 by default.
    Returns:
        dcbc_values (np.ndarray): the DCBC values of subjects
    """

    dist = dcbc.compute_dist(atlas.vox.T, resolution=resolution)

    if trim_nan:  # mask the nan voxel pairs distance to nan
        dist[np.where(np.isnan(parcels))[0], :] = np.nan
        dist[:, np.where(np.isnan(parcels))[0]] = np.nan

    dcbc_values = []
    for sub in range(testdata.shape[0]):
        D = dcbc.compute_DCBC(parcellation=parcels,
                              dist=dist, func=testdata[sub].T)
        dcbc_values.append(D['DCBC'])

    return np.asarray(dcbc_values)

def eval1():
    model_name = ['asym_Md_space-MNISymC3_K-10',
                   'asym_Po_space-MNISymC3_K-10',
                   'asym_Ni_space-MNISymC3_K-10',
                   'asym_MdPoNi_space-MNISymC3_K-10']
    
    # plot_parcel_flat_best(model_name,[2,2])
    R = calc_prediction_error(model_name,
                         test_data='Mdtb',
                         test_sess=['ses-s1','ses-s2'],
                         design_ind='cond_num_uni',
                         part_ind='half',
                         indivtrain_ind='sess',
                         indivtrain_values=['ses-s1','ses-s2'])
    R.to_csv(base_dir + '/Models/eval_Mdtb.tsv',sep='\t')

def eval2():
    space = 'MNISymC3'
    model_name = [f'asym_Md_space-{space}_K-10',
                  f'asym_Po_space-{space}_K-10',
                   f'asym_Ni_space-{space}_K-10',
                   f'asym_MdPoNi_space-{space}_K-10']

    allR = pd.DataFrame()
    for testdata in ['Md', 'Po', 'Ni']:
        R = eval_dcbc_group(model_name, space,
                                    testdata)
        # R.to_csv(base_dir + f'/Models/eval2_{testdata}.tsv', sep='\t')
        allR = pd.concat([allR, R], ignore_index=True)

    allR.to_csv(base_dir + f'/Models/eval_dcbc_group.tsv', sep='\t')

    # test_data = 'Mdtb'
    # R = calc_dcbc(model_name,
    #               test_data=test_data,
    #               test_sess=['ses-s1', 'ses-s2'],
    #               design_ind='cond_num_uni',
    #               part_ind='half',
    #               indivtrain_ind='sess',
    #               indivtrain_values=['ses-s1', 'ses-s2'])
    # R.to_csv(base_dir + f'/Models/eval2.tsv', sep='\t')
    
    # print(R)

    pass



def eval_generative_SNMF(model_names = ['asym_Md_space-SUIT3_K-10']):
    """This is the evaluation case of the parcellation comparison
    between the new fusion model vs. convex semi non-negative matrix
    factorization (King et. al, 2019).

    Args:
        model_names (list): the list of model names to be evaluated
    Returns:
        the plot
    Note:
        I'm just simply curious whether the fusion model fit on mdtb
        standalone (indepent ar. + VMF em.) can beat NMF algorithm
        or not. So nothing hurt the script.    -- dzhi
    """
    # Use specific mask / atlas.
    mask = base_dir + '/Atlases/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    atlas = am.AtlasVolumetric('SUIT3', mask_img=mask)

    # get original mdtb parcels (nmf)
    from learn_mdtb import get_mdtb_parcel
    mdtb_par, _ = get_mdtb_parcel(do_plot=False)
    mdtb_par = np.where(mdtb_par == 0, np.nan, mdtb_par)

    parcel = np.empty((len(model_names), atlas.P))

    for i, mn in enumerate(model_names):
        info, models, Prop, V = load_batch_fit(mn)
        j = np.argmax(info.loglik)
        # Get winner take all
        par = pt.argmax(Prop[j, :, :], dim=0) + 1
        parcel[i, :] = np.where(np.isnan(mdtb_par), np.nan, par.numpy())

    # Evaluate case: use all MDTB data
    # It kinda of overfitting but still fair comparison
    data_eval, _, _ = get_all_mdtb(atlas='SUIT3')
    dcbc_base = eval_dcbc(mdtb_par, atlas, func_data=data_eval,
                          resolution=3, trim_nan=True)

    dcbc_compare = []
    for p in range(parcel.shape[0]):
        this_dcbc = eval_dcbc(parcel[p], atlas, func_data=data_eval,
                              resolution=3, trim_nan=True)
        dcbc_compare.append(this_dcbc)

    plt.figure()
    plt.bar(['NMF', 'ind+vmf'], [dcbc_base.mean(), dcbc_compare[0].mean()],
            yerr=[dcbc_base.std() / np.sqrt(24),
                  dcbc_compare[0].std() / np.sqrt(24)])
    plt.show()




if __name__ == "__main__":
    eval2()
    pass

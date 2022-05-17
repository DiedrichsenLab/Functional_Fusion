# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
from pathlib import Path
import numpy as np
import atlas_map as am
from dataset import DataSetMDTB, DataSetHcpResting
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
import sys

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

if sys.platform == "win32":
    base_dir = 'Y:\data\FunctionalFusion'
    sys.path.append('../')

hcp_dir = base_dir + '/HCP'
data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

# Globals are ugly, but somewhat usegul here 
mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)

mdtb_dataset = DataSetMDTB(data_dir)
hcp_dataset = DataSetHcpResting(base_dir + '/HCP')

def show_mdtb_suit(subj,sess,cond): 
    T = mdtb_dataset.get_participants()
    s = T.participant_id[subj]
    ses = f'ses-s{sess}'
    C = nb.load(mdtb_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_CondSes.dscalar.nii')
    D = pd.read_csv(mdtb_dataset.data_dir.format(s) + f'/{s}_{ses}_info-CondSes.tsv',sep='\t')
    X = C.get_fdata()
    Nifti = suit_atlas.data_to_nifti(X)
    surf_data = suit.flatmap.vol_to_surf(Nifti)
    fig = suit.flatmap.plot(surf_data[:,cond],render='plotly')
    fig.show()
    print(f'Showing {D.cond_name[cond]}')
    pass 

def get_mdtb_data(ses_id='ses-s1', type='CondSes'):
    T = mdtb_dataset.get_participants()
    
    # Assemble the data 
    Data = None 
    for i,s in enumerate (T.participant_id): 

        C = nb.load(mdtb_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses_id}_{type}.dscalar.nii')
        if Data is None:
            Data = np.zeros((len(T.participant_id),C.shape[0],C.shape[1]))
        Data[i,:,:] = C.get_fdata()

    # Get design matrix 
    D = pd.read_csv(mdtb_dataset.data_dir.format(s) + f'/{s}_{ses_id}_info-{type}.tsv',sep='\t')
    Xcond = matrix.indicator(D.cond_num)

    # center the data for each voxel for each half of the experiment 
    if type=='CondSes':
        Xmean = matrix.indicator(D.half)
        Data -= (Xmean @ np.linalg.pinv(Xmean) @ Data)
    elif type=='CondRun':
        Xmean = matrix.indicator(D.run)
        Data -= (Xmean @ np.linalg.pinv(Xmean) @ Data)

    return Data, Xcond, D

def get_hcp_data(tessel=162, ses_id=['ses-01'], range=None, save=False):
    """Get the HCP resting-state connnectivity profile
       The running time of this function is very slow (~1m per subject/sess)
    Args:
        tessel: the cortical map used to generate connectivity
        ses_id: session id
        range: the index of subject among unrelated 100 dataset
               default - None, which means get all participants
    Returns:
        the HCP rs-FC, shape (n_subj, N, P) - N is based on tessel
    """
    suit_atlas = am.AtlasVolumetric('SUIT', mask_img=mask)
    deform = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
    MNI_mask = base_dir + '/Atlases/tpl-MNI152NLin6AsymC/tpl-MNI152AsymC_res-2_gmcmask2.nii'

    # get the tessellation file
    tessel_dir = atlas_dir + '/tpl-fs32k'
    # get the gifti file for the label
    gii_labels = [nb.load(tessel_dir + f'/Icosahedron-{tessel}.32k.L.label.gii'),
                  nb.load(tessel_dir + f'/Icosahedron-{tessel}.32k.R.label.gii')]

    # get the gifti of the mask
    gii_mask = [nb.load(atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii'),
                nb.load(atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii')]

    T = hcp_dataset.get_participants()

    if range is not None:
        id = T.participant_id.values[range]
    else:
        id = T.participant_id.values

    output = []
    for s in id:
        con_fp = []  # connectivity finger print list
        for sess in ses_id:
            # create a mapping based on atlas deformation
            atlas_map = am.AtlasMapDeform(hcp_dataset, suit_atlas, s, deform, MNI_mask)
            # print(f"building atlas map")
            atlas_map.build(smooth=2.0)  # smoothing level?

            # get the connectivity matrix and put it in a list
            data = hcp_dataset.get_data(s, [atlas_map], gii_labels, gii_mask, sess)

            if save:
                # save as np.ndarry as .csv under the same folder in fusion
                print(f'Saving sub-{s}, {sess} rs-FC')
                target_dir = hcp_dir + f'/derivatives/{s}/func/{sess}'
                pd.DataFrame(data.T).to_csv(target_dir+f'/sub-{s}_{sess}_tessel-{tessel}_conn.csv')

            con_fp.append(pt.tensor(data.T))

        sub_data = pt.vstack(con_fp)
        output.append(sub_data)

    return pt.stack(output)  # output is torch tensor

def get_hcp_data_from_csv(tessel=162, ses_id=['ses-01'], range=None):

    T = hcp_dataset.get_participants()

    if range is not None:
        id = T.participant_id.values[range]
    else:
        id = T.participant_id.values

    output = []
    for s in id:
        con_fp = []  # connectivity finger print list
        for sess in ses_id:
            print(f'Loading sub-{s}, {sess} rs-FC')
            target_dir = hcp_dir + f'/derivatives/{s}/func/{sess}'
            data = pd.read_csv(target_dir+f'/sub-{s}_{sess}_tessel-{tessel}_conn.csv', index_col=0)

            con_fp.append(pt.tensor(np.asarray(data)))

        sub_data = pt.vstack(con_fp)
        output.append(sub_data)

    return pt.stack(output)  # output is torch tensor

def get_mdtb_parcel(do_plot=True):
    """Samples the existing MDTB10 parcellation
    Then displays it as check 
    """
    parcel = nb.load(atlas_dir + '/tpl-SUIT/atl-MDTB10_space-SUIT_dseg.nii')
    T = mdtb_dataset.get_participants()
    data = suit.reslice.sample_image(parcel,
            suit_atlas.world[0],
            suit_atlas.world[1],
            suit_atlas.world[2],0)
    
    # Read the MDTB colors: Add additional row for parcel 0 
    color_file = atlas_dir + '/tpl-SUIT/atl-MDTB10.lut'
    color_info = pd.read_csv(color_file, sep = ' ',header=None)
    MDTBcolors = np.zeros((11,3))
    MDTBcolors[1:11,:]  = color_info.iloc[:,1:4].to_numpy()
    
    # Map Plot if requested (for a check) 
    if do_plot: 
        Nifti = suit_atlas.data_to_nifti(data)
        surf_data = suit.flatmap.vol_to_surf(Nifti,stats='mode')
        fig = suit.flatmap.plot(surf_data,render='plotly',overlay_type='label',cmap= MDTBcolors)
        fig.show()
    return data,MDTBcolors

def _plot_maps(data, sub=None, stats='mode', render_type='plotly',
               overlay='label', color=None, save=None):
    # Read the MDTB colors: Add additional row for parcel 0
    color_file = atlas_dir + '/tpl-SUIT/atl-MDTB10.lut'
    color_info = pd.read_csv(color_file, sep=' ', header=None)
    MDTBcolors = np.zeros((11, 3))
    MDTBcolors[1:11, :] = color_info.iloc[:, 1:4].to_numpy()

    if sub is not None:
        Nifti = suit_atlas.data_to_nifti(data[sub])
    else:
        Nifti = suit_atlas.data_to_nifti(data)

    surf_data = suit.flatmap.vol_to_surf(Nifti, stats=stats)
    if color is not None:
        fig = suit.flatmap.plot(surf_data, render=render_type,
                                overlay_type=overlay, cmap=MDTBcolors)
    else:
        fig = suit.flatmap.plot(surf_data, render=render_type,
                                overlay_type=overlay)

    if save is not None:
        fig.write_image(save)

    fig.show()

def learn_single(ses_id='ses-s1', max_iter=100, fit_arr=False):
    """Learn a single data set 
    """
    # Data_HCP = get_hcp_data(ses_id=['ses-01','ses-02'], range=np.arange(2), save=True)
    Data_HCP = get_hcp_data_from_csv(ses_id=['ses-01','ses-02'], range=np.arange(2))
    Data,Xdesign = get_mdtb_data(ses_id)
    P = Data.shape[2]
    K = 10 

    # Make arrangement model and initialize the prior from the MDTB map
    prior_w = 7.0 # Weight of prior
    mdtb_parcel,mdtb_colors = get_mdtb_parcel()
    logpi = ar.expand_mn(mdtb_parcel.reshape(1,P)-1,K)
    logpi = logpi.squeeze()*prior_w
    # Set parcel 0 to unassigned 
    logpi[:,mdtb_parcel==0]=0

    # Making emission model for fitting, the input N here
    # doesn't matter and will be overwrite internally by X
    ar_model = ar.ArrangeIndependent(K=K, P=P, spatial_specific=True,
                                         remove_redundancy=False)
    em_model = em.MixVMF(K=K, N=40, P=P, X=Xdesign, uniform_kappa=True)
    em_model.initialize(Data)

    # Initilize parameters from prior
    group_prior = logpi.softmax(dim=0).unsqueeze(0).repeat(em_model.num_subj,1,1)
    em_model.Mstep(group_prior)
    M = fm.FullModel(ar_model, em_model)

    # Step 5: Estimate the parameter thetas to fit the new model using EM
    M, ll, theta, U_hat = M.fit_em(Y=Data, iter=max_iter, tol=0.00001, fit_arrangement=True)

    # plot emission log-likelihood
    plt.plot(ll, color='b')
    plt.show()

    # PLOT 1 - group logpi
    _plot_maps(pt.argmax(M.arrange.logpi, dim=0)+1, save="group_logpi.pdf")

    # Make inference on first half of the data - now Data shape (n_sub, 29, P)
    M.emission.X = pt.tensor(Xdesign[0:29,:])
    U_hat_em = M.emission.Estep(Y=Data[:,0:29,:])
    U_hat_complete, _ = M.Estep(Y=Data[:,0:29,:])

    for s in [8,9,10,11]:
        # PLOT 2 - top row: the indiviudal Uhat without group logpi
        _plot_maps(pt.argmax(U_hat_em, dim=1)+1, sub=s)

        # PLOT 3 - bottom row: the indiviudal Uhat + group logpi
        _plot_maps(pt.argmax(U_hat_complete, dim=1)+1, sub=s)

    pass

def learn_runs(K=10, max_iter=100, run_test=np.arange(58, 122),
               runs=np.arange(1, 17), sub=None, do_plot=True):
    Data_s1, Xdesign_s1, _ = get_mdtb_data('ses-s1')
    Data_s2, Xdesign_s2, _ = get_mdtb_data('ses-s2')
    Data = np.concatenate((Data_s1, Data_s2), axis=1)

    # Create combined X design matrix
    # Z1 = np.zeros((Xdesign_s1.shape[0], Xdesign_s2.shape[1]), dtype=Xdesign_s1.dtype)
    # Z2 = np.zeros((Xdesign_s2.shape[0], Xdesign_s1.shape[1]), dtype=Xdesign_s2.dtype)
    # Xdesign = np.asarray(np.bmat([[Xdesign_s1, Z1], [Z2, Xdesign_s2]]))
    Xdesign = pt.tensor(block_diag(Xdesign_s1, Xdesign_s2))
    del Data_s1, Data_s2, Xdesign_s1, Xdesign_s2

    # Make arrangement model and initialize the prior from the MDTB map
    P = Data.shape[2]
    prior_w = 7.0  # Weight of prior
    mdtb_parcel, mdtb_colors = get_mdtb_parcel(do_plot=False)
    logpi = ar.expand_mn(mdtb_parcel.reshape(1, P) - 1, K)
    logpi = logpi.squeeze() * prior_w
    # Set parcel 0 to unassigned
    logpi[:, mdtb_parcel == 0] = 0

    # Making emission model for fitting, the input N here
    # doesn't matter and will be overwrite internally by X
    ar_model = ar.ArrangeIndependent(K=K, P=P, spatial_specific=True,
                                         remove_redundancy=False)
    em_model = em.MixVMF(K=K, N=40, P=P, X=Xdesign, uniform_kappa=True)
    em_model.initialize(Data)

    # Initilize parameters from group prior and train the model
    group_prior = logpi.softmax(dim=0).unsqueeze(0).repeat(em_model.num_subj,1,1)
    em_model.Mstep(group_prior)
    M = fm.FullModel(ar_model, em_model)
    M, ll, theta, U_hat = M.fit_em(Y=Data, iter=max_iter, tol=0.00001, fit_arrangement=True)

    # plot emission log-likelihood and group prior
    plt.plot(ll, color='b')
    plt.show()
    # _plot_maps(pt.argmax(M.arrange.logpi, dim=0) + 1, save="group_logpi.pdf")
    prior = pt.softmax(M.arrange.logpi, dim=0).unsqueeze(0).repeat(Data.shape[0], 1, 1)
    M.emission.X = Xdesign[run_test, :]
    group_baseline = ev.coserr(pt.tensor(Data[:, run_test, :]),
                               pt.matmul(Xdesign[run_test], M.emission.V), prior,
                               adjusted=True, soft_assign=True)
    lower_bound = ev.coserr(pt.tensor(Data[:, run_test, :]),
                            pt.matmul(Xdesign[run_test], M.emission.V),
                            pt.softmax(M.emission.Estep(Y=Data[:,run_test,:]), dim=1),
                            adjusted=True, soft_assign=True)

    ############################# Cross-validation starts here #############################
    # Make inference on the selected run's data - run_infer
    Data_cv, Xdesign_cv, D = get_mdtb_data(ses_id='ses-s1', type='CondRun')
    Xdesign_cv = pt.tensor(np.concatenate((Xdesign_cv,np.zeros((Xdesign_cv.shape[0], 32))), axis=1))
    indices, cos_em, cos_complete, uhat_em_all, uhat_complete_all = [],[],[],[],[]
    for i in runs:
        indices.append(np.asarray(np.where(D.run==i)).reshape(-1))
        acc_run_idx = np.concatenate(indices).ravel()
        M.emission.X = Xdesign_cv[acc_run_idx,:]

        # Infer on current accumulated runs data
        U_hat_em = M.emission.Estep(Y=Data_cv[:,acc_run_idx,:])
        U_hat_complete, _ = M.Estep(Y=Data_cv[:,acc_run_idx,:])
        uhat_em_all.append(U_hat_em)
        uhat_complete_all.append(U_hat_complete)

        # Calculate cosine error/u abs error between another test data
        # and U_hat inferred from testing
        coserr_Uem = ev.coserr(pt.tensor(Data[:,run_test,:]),
                               pt.matmul(Xdesign[run_test], M.emission.V),
                               pt.softmax(U_hat_em, dim=1),
                               adjusted=True, soft_assign=True)
        coserr_Uall = ev.coserr(pt.tensor(Data[:,run_test,:]),
                                pt.matmul(Xdesign[run_test], M.emission.V), U_hat_complete,
                                adjusted=True, soft_assign=True)

        cos_em.append(coserr_Uem)
        cos_complete.append(coserr_Uall)
        # uerr_Uem = ev.rmse_YUhat(U_hat_em, pt.tensor(Data[:, 58:90, :]),
        #                          M.emission.V[29:61, :])
        # uerr_Uall = ev.rmse_YUhat(U_hat_complete, pt.tensor(Data[:, 58:90, :]),
        #                          M.emission.V[29:61, :])

    return group_baseline, lower_bound, cos_em, cos_complete, uhat_em_all, uhat_complete_all

if __name__ == "__main__":
    # Data_HCP = get_hcp_data(ses_id=['ses-01', 'ses-02'], range=np.arange(63, 64), save=True)
    # data = get_hcp_data_from_csv(tessel=162, ses_id=['ses-01'], range=[0])
    # data = data[0, :, :].cpu().detach().numpy()
    # _plot_maps(data, sub=1, stats='nanmean', overlay='func', color=None, save=None)
    gbase, lb, cos_em, cos_complete, uhat_em_all, uhat_complete_all = learn_runs(K=10,
                                                                          runs=np.arange(1, 17))

    plt.figure()
    x = np.arange(len(np.arange(1, 17)))

    plt.errorbar(x, pt.stack(cos_em).mean(dim=1),
             yerr=pt.stack(cos_em).std(dim=1)/np.sqrt(24), capsize=10, label='emission only')
    plt.errorbar(x, pt.stack(cos_complete).mean(dim=1),
             yerr=pt.stack(cos_complete).std(dim=1)/np.sqrt(24), capsize=10, label='emi + prior')
    plt.axhline(y=lb.mean(), color='r', linestyle=':', label='lower bound')
    plt.axhline(y=gbase.mean(), color='k', linestyle=':', label='goup prior')
    # plt.xticks(x, ['s1-r1', 's1-r2', 's1-rall', 's2-rall'])
    plt.xlabel('Inferred on')
    plt.ylabel('Adjusted cosine error')
    # plt.ylim(0, 0.35)
    plt.legend(loc='upper right')
    plt.title('test on session 2 - all runs')

    plt.show()
    pass

# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
from pathlib import Path
import numpy as np
import atlas_map as am
from dataset import DataSetMDTB, DataSetHcpResting
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

def get_mdtb_data(ses_id='ses-s1'):
    T = mdtb_dataset.get_participants()
    
    # Assemble the data 
    Data = None 
    for i,s in enumerate (T.participant_id): 

        C = nb.load(mdtb_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses_id}_CondSes.dscalar.nii')
        if Data is None:
            Data = np.zeros((len(T.participant_id),C.shape[0],C.shape[1]))
        Data[i,:,:] = C.get_fdata()

    # Get design matrix 
    D = pd.read_csv(mdtb_dataset.data_dir.format(s) + f'/{s}_{ses_id}_info-CondSes.tsv',sep='\t')
    Xcond = matrix.indicator(D.cond_num)

    # center the data for each voxel for each half of the experiment 
    Xhalf = matrix.indicator(D.half)
    Data -=(Xhalf@np.linalg.pinv(Xhalf)@Data)

    return Data, Xcond

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

def _plot_maps(data, sub=None, render_type='plotly', save=None):
    # Read the MDTB colors: Add additional row for parcel 0
    color_file = atlas_dir + '/tpl-SUIT/atl-MDTB10.lut'
    color_info = pd.read_csv(color_file, sep=' ', header=None)
    MDTBcolors = np.zeros((11, 3))
    MDTBcolors[1:11, :] = color_info.iloc[:, 1:4].to_numpy()

    if sub is not None:
        Nifti = suit_atlas.data_to_nifti(data[sub])
    else:
        Nifti = suit_atlas.data_to_nifti(data)

    surf_data = suit.flatmap.vol_to_surf(Nifti, stats='mode')
    fig = suit.flatmap.plot(surf_data, render=render_type,
                            overlay_type='label', cmap=MDTBcolors)

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

if __name__ == "__main__":
    learn_single()
    pass

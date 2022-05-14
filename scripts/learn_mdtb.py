# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
from pathlib import Path
import numpy as np
import atlas_map as am
from dataset import DataSetMDTB
import nibabel as nb
import SUITPy as suit
import generativeMRF.full_model as fm
import generativeMRF.spatial as sp
import generativeMRF.arrangements as ar
import generativeMRF.emissions as em
import torch as pt
import Functional_Fusion.matrix as matrix 

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

# Globals are ugly, but somewhat usegul here 
mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
mdtb_dataset = DataSetMDTB(data_dir)

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

def learn_single(ses_id = 'ses-s1'):
    """Learn a single data set 
    """

    Data,Xdesign = get_mdtb_data(ses_id)
    P = Data.shape[2]
    K = 10 

    # Make arrangement model and initialize the prior from the MDTB map 
    armodel = ar.ArrangementModel(K,P)
    prior_w = 3.0 # Weight of prior 
    mdtb_parcel,mdtb_colors = get_mdtb_parcel()
    logpi = ar.expand_mn(mdtb_parcel.reshape(1,P)-1,K)
    logpi = logpi.squeeze()*prior_w
    # Set parcel 0 to unassigned 
    logpi[:,mdtb_parcel==0]=0
    armodel.logpi = logpi

    pass

if __name__ == "__main__":
    learn_single()
    pass

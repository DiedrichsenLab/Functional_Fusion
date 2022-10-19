# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetMDTB
import nibabel as nb
import SUITPy as suit


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

data_dir = base_dir + '/IBC'
atlas_dir = base_dir + '/Atlases'

sess =['archi','hcp1','hcp2']

def extract_ibc_suit(ses_id='ses-archi',type='condHalf',atlas='SUIT3'):
    ibc_dataset = DataSetIBC(data_dir)
    ibc_dataset.extract_all_suit(ses_id,type,atlas)

def extract_mdtb_fs32k(ses_id='ses-s1',type='condHalf'):
    ibc_dataset = DataSetIBC(data_dir)
    ibc_dataset.extract_all_fs32k(ses_id,type)

def show_mdtb_suit(subj,sess,cond):
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    mdtb_dataset = DataSetMDTB(data_dir)
    T = mdtb_dataset.get_participants()
    ses = f'ses-s{sess}'

    if subj=='all':
        for i,s in enumerate(T.participant_id):
            C = nb.load(mdtb_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_CondAll.dscalar.nii')
            if i==0:
                X = np.zeros((C.shape[0],C.shape[1],T.shape[0]))
            X[:,:,i] = C.get_fdata()
        X = X.mean(axis=2)
        D = pd.read_csv(mdtb_dataset.data_dir.format(s) + f'/{s}_{ses}_info-CondAll.tsv',sep='\t')

    else:
        s = T.participant_id[subj]
        C = nb.load(mdtb_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_CondAll.dscalar.nii')
        D = pd.read_csv(mdtb_dataset.data_dir.format(s) + f'/{s}_{ses}_info-CondAll.tsv',sep='\t')
        X = C.get_fdata()
    Nifti = suit_atlas.data_to_nifti(X[cond,:])
    surf_data = suit.flatmap.vol_to_surf(Nifti)
    fig = suit.flatmap.plot(surf_data,render='plotly',cscale=[-0.3,0.3])
    fig.show()
    print(f'Showing {D.cond_name[cond]}')
    pass





if __name__ == "__main__":
    # parcel_mdtb_fs32k()
    for i in sess_ids:
        extract_ibc_suit(ses_id=i,type='CondHalf',atlas='MNISymC3')
        extract_ibc_fs32k(ses_id='ses-s1',type='CondAll')
    # extract_mdtb_suit(ses_id='ses-s1',type='CondAll')
    # extract_mdtb_suit(ses_id='ses-s2',type='CondAll')
    # 
    # extract_mdtb_fs32k(ses_id='ses-s2',type='condHalf')
    # show_mdtb_suit('all',1,0)
    pass


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
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
import matplotlib.pyplot as plt

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:\data\FunctionalFusion'
if not Path(base_dir).exists():
    raise(NameError('Could not find base_dir'))

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'


def group_mtdb(ses_id='ses-s1', type='CondHalf', atlas='SUIT3'):
    mdtb_dataset = DataSetMDTB(data_dir)
    mdtb_dataset.group_average_data(ses_id,type, atlas)

def extract_mdtb(ses_id='ses-s1', type='CondHalf', atlas='SUIT3'):
    mdtb_dataset = DataSetMDTB(data_dir)
    mdtb_dataset.extract_all(ses_id,type,atlas)

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


def show_mdtb_group(type='CondHalf', atlas='SUIT3', cond=0, info_column='cond_name', savefig=False):
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum', mask_img=mask)
    mdtb_dataset = DataSetMDTB(data_dir)
    C = nb.load(mdtb_dataset.data_dir.split('/{0}')[0] +
                f'/group/group_space-{atlas}_{type}.dscalar.nii')
    D = pd.read_csv(mdtb_dataset.data_dir.split('/{0}')[0] +
                    f'/group/group_info-{type}.tsv', sep='\t')
    X = C.get_fdata()

    if cond=='all':
        conditions = D[info_column]
        # -- as subplot --
        # dim = int(np.ceil(np.sqrt(len(conditions))))
        # fig, axs = plt.subplots(dim,dim)
        # for i, c in enumerate(conditions):
        #     Nifti = suit_atlas.data_to_nifti(X[i, :])
        #     surf_data = suit.flatmap.vol_to_surf(Nifti)
        #     axs[int(i % dim), int(i / dim)]=suit.flatmap.plot(
        #         surf_data, render='matplotlib', new_figure=False)
        #     axs[int(i % dim), int(i / dim)].set_title(c)
        # fig.show()
        # -- each in seperate figures --
        dest_dir = mdtb_dataset.data_dir.split('/{0}')[0] + f'/group/figures/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        for i, c in enumerate(conditions):
            Nifti = suit_atlas.data_to_nifti(X[i, :])
            surf_data = suit.flatmap.vol_to_surf(Nifti)
            fig = suit.flatmap.plot(
                surf_data, render='matplotlib', new_figure=True)
            fig.set_title(c)
            # save figure
            if savefig:
                plt.savefig(dest_dir + f'group_{c}.png')
            plt.clf()
            pass

    else:
        Nifti = suit_atlas.data_to_nifti(X[cond, :])
        surf_data = suit.flatmap.vol_to_surf(Nifti)
        fig = suit.flatmap.plot(surf_data, render='plotly')
        fig.show()
        print(f'Showing {D.cond_name[cond]}')
        pass


def parcel_mdtb_fs32k(res=162,ses_id='ses-s1',type='condHalf'):
    # Make the atlas object
    surf_parcel =[]
    hem_name = ['cortex_left','cortex_right']
    # Get the parcelation
    for i,h in enumerate(['L','R']):
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i],gifti))

    # initialize the data set object
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    for s in T.participant_id[0:2]:
        print(f'Average {s}')
        s_dir = mdtb_dataset.data_dir.format(s)
        C = nb.load(s_dir + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
        bmf = C.header.get_axis(1)
        bmp = []
        R = []
        for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
            D = np.asanyarray(C.dataobj[:,slc])
            X = np.zeros((D.shape[0],surf_parcel[0].label_vec.shape[0]))
            X[:,bm.vertex]=D
            R.append(surf_parcel[idx].agg_data(X))
            bmp.append(surf_parcel[idx].get_parcel_axis())
        header = nb.Cifti2Header.from_axes((C.header.get_axis(0),bmp[0]+bmp[1]))
        cifti_img = nb.Cifti2Image(dataobj=np.c_[R[0],R[1]],header=header)
        nb.save(cifti_img,s_dir + f'/{s}_space-fs32k_{ses_id}_{type}_Iso-{res}.pscalar.nii')
        pass



if __name__ == "__main__":
    extract_mdtb(ses_id='ses-s1', type='CondHalf',atlas='MNISymC2')
    # show_mdtb_group(type='CondHalf', atlas='SUIT3', cond='all', savefig=True)
    pass



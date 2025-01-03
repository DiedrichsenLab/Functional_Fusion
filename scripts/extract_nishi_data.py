# Script for importing the nishi data set from super_cerebellum to general format.
import os, time
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.util as fut
from Functional_Fusion.dataset import DataSetNishi
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import subprocess

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    raise(NameError('Could not find base_dir'))

data_dir = base_dir + '/Nishimoto'
atlas_dir = base_dir + '/Atlases'


def extract_nishi_group(type='CondHalf', atlas='SUIT3', info_column='task_name'):
    nishi_dataset = DataSetNishi(data_dir)
    nishi_dataset.group_average_data(type, atlas, info_column)

def extract_nishi_suit(ses_id='ses-01',type='CondHalf', atlas= 'SUIT3'):
    nishi_dataset = DataSetNishi(data_dir)
    nishi_dataset.extract_all_suit(ses_id,type,atlas)
    
def extract_nishi_fs32k(ses_id='ses-01',type='CondHalf'):
    nishi_dataset = DataSetNishi(data_dir)
    nishi_dataset.extract_all_fs32k(ses_id,type)

def show_nishi_suit(subj,ses,cond):
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    nishi_dataset = DataSetNishi(data_dir)
    T = nishi_dataset.get_participants()
    s = T.participant_id[subj]
    ses = f'ses-{ses:02d}'
    C = nb.load(nishi_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_CondHalf.dscalar.nii')
    D = pd.read_csv(nishi_dataset.data_dir.format(s) + f'/{s}_{ses}_info-CondHalf.tsv',sep='\t')
    X = C.get_fdata()
    Nifti = suit_atlas.data_to_nifti(X)
    surf_data = suit.flatmap.vol_to_surf(Nifti)
    fig = suit.flatmap.plot(surf_data[:,cond],render='plotly')
    fig.show()
    print(f'Showing {D.task_name[cond]}')
    pass

def parcel_nishi_fs32k(res=162,ses_id='ses-01',type='condHalf'):
    # Make the atlas object
    surf_parcel =[]
    hem_name = ['cortex_left','cortex_right']
    # Get the parcelation
    for i,h in enumerate(['L','R']):
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i],gifti))

    # initialize the data set object
    nishi_dataset = DataSetNishi(data_dir)

    # create and calculate the atlas map for each participant
    T = nishi_dataset.get_participants()
    for s in T.participant_id:
        print(f'Average {s}')
        s_dir = nishi_dataset.data_dir.format(s)
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


def smooth_nishi_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1, kernel='gaussian'):
    dataset = DataSetNishi(data_dir)
    # get the surfaces for smoothing
    surf_L = dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.L.midthickness.surf.gii'
    surf_R = dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.R.midthickness.surf.gii'

    T = dataset.get_participants()
    for s in T.participant_id:
        print(f'Smoothing data for {s} fs32k {ses_id} in {smooth}mm {kernel} ...')

        start = time.perf_counter()
        file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        fut.smooth_fs32k_data(file, surf_L, surf_R, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")

if __name__ == "__main__":
    smooth_nishi_fs32k(ses_id='ses-01', type='CondHalf', smooth=4, kernel='fwhm')
    smooth_nishi_fs32k(ses_id='ses-02', type='CondHalf', smooth=4, kernel='fwhm')
    # extract_nishi_group(type='CondHalf', atlas='SUIT3')
    # show_nishimoto_group(type='CondHalf', atlas='SUIT3', cond='all', savefig=True)
    # parcel_nishi_fs32k(res=162,ses_id='ses-01',type='condHalf')
    # parcel_nishi_fs32k(res=362,ses_id='ses-01',type='condHalf')
    # parcel_nishi_fs32k(res=642,ses_id='ses-01',type='condHalf')
    # parcel_nishi_fs32k(res=162,ses_id='ses-02',type='condHalf')
    # parcel_nishi_fs32k(res=362,ses_id='ses-02',type='condHalf')
    # parcel_nishi_fs32k(res=642,ses_id='ses-02',type='condHalf')
    # extract_nishi_suit(ses_id='ses-01',type='CondHalf',atlas='MNISymC3')
    # extract_nishi_suit(ses_id='ses-02',type='CondHalf',atlas='MNISymC3')
    # extract_nishi_fs32k(ses_id='ses-01',type='CondHalf')
    # extract_nishi_fs32k(ses_id='ses-02',type='CondHalf')

    dataset = DataSetNishi(data_dir)
    dataset.extract_all(ses_id='ses-01', type='CondAll', atlas='MNISymC3')
    dataset.extract_all(ses_id='ses-02', type='CondAll', atlas='MNISymC3')
    # dataset.group_average_data(type='CondHalf', atlas='MNISymC3')
    # dataset.plot_cerebellum(savefig=True, atlas='MNISymC3', colorbar=True)

    ## fs32k smoothing ##
    # for s in [1,2,3,4,5,6,7]:
    #     smooth_nishi_fs32k(ses_id='ses-01', type='CondHalf', smooth=s)
    #     smooth_nishi_fs32k(ses_id='ses-02', type='CondHalf', smooth=s)
    pass

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
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    raise(NameError('Could not find base_dir'))

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'


def group_mtdb(ses_id='ses-s1', type='CondHalf', atlas='SUIT3'):
    mdtb_dataset = DataSetMDTB(data_dir)
    mdtb_dataset.group_average_data(ses_id,type, atlas)

def extract_mdtb(ses_id='ses-s1', type='CondHalf', atlas='SUIT3', smooth=2.0):
    mdtb_dataset = DataSetMDTB(data_dir)
    mdtb_dataset.extract_all(ses_id,type,atlas, smooth=smooth)

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
    extract_mdtb(ses_id='ses-s1', type='CondHalf', atlas='MNISymC3', smooth=None)
    extract_mdtb(ses_id='ses-s2', type='CondHalf', atlas='MNISymC3', smooth=None)
    # show_mdtb_group(type='CondHalf', atlas='SUIT3', cond='all', savefig=True)
    

    # dataset = DataSetMDTB(data_dir)
    # # dataset.group_average_data(atlas='MNISymC3')
    # dataset.group_average_data(atlas='MNISymC3', ses_id='ses-s2')
    #
    # dataset.plot_cerebellum(savefig=True, atlas='MNISymC3',
    #                         colorbar=True, sessions=['ses-s2'])
    #
    # pass



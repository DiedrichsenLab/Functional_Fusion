# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetMDTB
import nibabel as nb
import subprocess
import scripts.paths as pt
import Functional_Fusion.connectivity as conn


base_dir = pt.set_base_dir()
atlas_dir = pt.set_atlas_dir(base_dir)
dname = 'MDTB'
data_dir = base_dir + '/' + dname
atlas_dir = base_dir + '/Atlases'


def group_mtdb(ses_id='ses-s1', type='CondHalf', atlas='SUIT3'):
    mdtb_dataset = DataSetMDTB(data_dir)
    mdtb_dataset.group_average_data(ses_id,type, atlas)


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

def smooth_mdtb_fs32k(ses_id='ses-s1', type='CondHalf', smooth=1):
    myatlas, _ = am.get_atlas('fs32k', atlas_dir)
    mdtb_dataset = DataSetMDTB(data_dir)
    T = mdtb_dataset.get_participants()

    # get the surfaces for smoothing
    surf_L = mdtb_dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.L.midthickness.surf.gii'
    surf_R = mdtb_dataset.atlas_dir + f'/tpl-fs32k/fs_LR.32k.R.midthickness.surf.gii'

    for s in T.participant_id:
        print(f'- Smoothing data for {s} fs32k {ses_id} in {smooth}mm ...')
        # Load the unsmoothed data and fill nan with zeros
        C = nb.load(mdtb_dataset.data_dir.format(s)
                    + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
        mask = np.isnan(C.get_fdata())
        C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
        nb.save(C, 'tmp.dscalar.nii')

        dest_dir = mdtb_dataset.data_smooth_dir.format(s, smooth)
        cifti_out = dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        # Write in smoothed surface data (filled with 0)
        smooth_cmd = f"wb_command -cifti-smoothing tmp.dscalar.nii " \
                     f"{smooth} {smooth} COLUMN {cifti_out} " \
                     f"-left-surface {surf_L} -right-surface {surf_R} " \
                     f"-fix-zeros-surface"
        subprocess.run(smooth_cmd, shell=True)
        os.remove("tmp.dscalar.nii")

        # Replace 0s back to NaN (we don't want the 0s impact model learning)
        C = nb.load(cifti_out)
        data = C.get_fdata()
        data[mask] = np.nan
        C = nb.Cifti2Image(dataobj=data, header=C.header)
        nb.save(C, cifti_out)

        # Copy info to the corresponding /smoothed folder
        if not Path(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv').exists():
            info = pd.read_csv(mdtb_dataset.data_dir.format(s)
                               + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)


if __name__ == "__main__":
    mdtb_dataset = DataSetMDTB(data_dir)
    # mdtb_dataset.extract_all(ses_id='ses-s1', type='CondRun', atlas='MNISymC2',smooth=2.0)
    # mdtb_dataset.extract_all(ses_id='ses-s1', type='CondHalf', atlas='fs32k', smooth=None)
    # mdtb_dataset.extract_all(ses_id='ses-s2', type='CondHalf', atlas='fs32k', smooth=None)
    # show_mdtb_group(type='CondHalf', atlas='SUIT3', cond='all', savefig=True)

    # for s in [1,2,3]:
    #     smooth_mdtb_fs32k(ses_id='ses-s1', type='CondHalf', smooth=s)
    #     smooth_mdtb_fs32k(ses_id='ses-s2', type='CondHalf', smooth=s)

    # # mdtb_dataset.group_average_data(atlas='MNISymC3')
    # mdtb_dataset.group_average_data(atlas='MNISymC3', ses_id='ses-s2')
    #
    # mdtb_dataset.plot_cerebellum(savefig=True, atlas='MNISymC3',
    #                         colorbar=True, sessions=['ses-s2'])
    #

    # -- Extract resting-state timeseries --
    # mdtb_dataset.extract_all(ses_id='ses-rest', type='Tseries', atlas='MNISymC2', smooth=2.0)
    # mdtb_dataset.extract_all(ses_id='ses-rest', type='Tseries', atlas='fs32k', smooth=2.0)

    # -- Get connectivity fingerprint --
    conn.get_connectivity_fingerprint(dname,
        type='Net69Run', space='MNISymC2', ses_id='ses-rest')
    conn.get_connectivity_fingerprint(dname,
        type='Net69Run', space='SUIT3', ses_id='ses-rest')
    conn.get_connectivity_fingerprint(dname,
        type='Net300Run', space='MNISymC2', ses_id='ses-rest')
    conn.get_connectivity_fingerprint(dname,
        type='Net300Run', space='SUIT3', ses_id='ses-rest')
    # mdtb_dataset.extract_all(ses_id='ses-rest', type='Net69Run', atlas='MNISymC2', smooth=2.0)


# Script for getting all the HCP data for cerebellar-cortical connectivity
import pandas as pd
import shutil, time, subprocess
from pathlib import Path
import mat73
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetHcpResting
import Functional_Fusion.dataset as ds
import nibabel as nb
import SUITPy as suit
import os
import sys
import matplotlib.pyplot as plt
import re
import Functional_Fusion.connectivity as conn

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/data/tge/Tian/UKBB_full/imaging'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

hcp_dir = '/data/tge/Tian/HCP_img'
atlas_dir = base_dir + '/Atlases'
hem_name = ['cortex_left', 'cortex_right']


def extract_hcp_timeseries(ses_id='ses-rest1', type='Tseries', atlas='MNISymC3'):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    hcp_dataset.extract_all(ses_id, type, atlas)


def make_info(type='Tseries', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file

    Args:
        type: file type
        ses_id: session id

    Returns:
        Write in the modified tsv file
    """
    hcp_dataset = DataSetHcpResting(hcp_dir)

    T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    for p, participant_id in enumerate(T.participant_id):
        # Make info
        dest_dir = hcp_dataset.base_dir + \
            f'/derivatives/{participant_id}/data/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        info = pd.read_csv(dest_dir + f'{participant_id}_{ses_id}_info-{type}.tsv',
                           sep='\t')

        info['run_id'] = info['run'].copy()
        if ses_id == 'ses-rest2':
            info['run_id'] = info['run_id'] + 2

        info.to_csv(
            dest_dir + f'{participant_id}_{ses_id}_info-{type}.tsv', sep='\t', index=False)

def rename_infofile():
    hcp_dataset = DataSetHcpResting(hcp_dir)

    T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    for p, participant_id in enumerate(T.participant_id):

        dest_dir = hcp_dataset.base_dir + \
                   f'/derivatives/{participant_id}/data/'
        for file_name in os.listdir(dest_dir):
            if "_info" in file_name and file_name.endswith(".tsv"):
                new_name = file_name.replace("info-", "")
                os.rename(
                    os.path.join(dest_dir, file_name),
                    os.path.join(dest_dir, new_name)
                )

def group_average_hcp(type='Net69Run', atlas='MNISymC3'):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    hcp_dataset.group_average_data(
        ses_id='ses-rest1', type=type, atlas=atlas)
    hcp_dataset.group_average_data(
        ses_id='ses-rest2', type=type, atlas=atlas)
    hcp_dataset.plot_cerebellum(subject='group', sessions=[
                                'ses-rest1', 'ses-rest2'], type=type, atlas=atlas, savefig=True, colorbar=True)
    # # get figures for each subject
    # T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    # for p, participant_id in enumerate(T.participant_id):
    #     hcp_dataset.plot_cerebellum(subject=participant_id, sessions=[
    #         'ses-rest1', 'ses-rest2'], type=type, atlas=atlas, savefig=True, colorbar=True)



def get_hcp_fs32k_rsfc(type='Net69Run', space='MNISymC3', ses_id='ses-rest1',
                       subj_list=None, smooth=None, kernel=None, thres=None, keeptop=False):
    # Load dataset
    if subj_list is None:
        dset = DataSetHcpResting(hcp_dir)
    else:
        dset = DataSetHcpResting(hcp_dir, subj_id_file=subj_list)
    T = dset.get_participants()
    
    # Load the cortical networks
    target, type = re.findall(r'[A-Z]+[a-z0-9]*', type)      
    res = target[3:]
    
    if target[:3] == 'Net':
        net = nb.load(hcp_dir + '/derivatives/group' +
                      f'/{target}_space-fs32k.dscalar.nii')
    elif target[:3] == 'Ico':
        net = [atlas_dir + f'/tpl-fs32k/Icosahedron{res}.L.label.gii',
            atlas_dir + f'/tpl-fs32k/Icosahedron{res}.R.label.gii']
    elif target[:3] == 'Fus':
        net = nb.load(hcp_dir + '/derivatives/group' +
                      f'/{target}_space-fs32k.pscalar.nii')
    atlas, _ = am.get_atlas(space)
        
    for i, s in enumerate(T.participant_id):
        dest_dir = dset.data_dir.format(s)
        if smooth is not None:
            file_name = f'{dest_dir}/{s}_space-{space}_{ses_id}'\
                        f'_{target+type}_desc-sm{smooth}'\
                        f'{kernel if kernel is not None else ""}'
        else:
            file_name = f'{dest_dir}/{s}_space-{space}_{ses_id}'\
                        f'_{target+type}'
        
        file_name = file_name + '_binarized' if thres is not None else file_name
        file_name = file_name + '_kt' if keeptop else file_name
        file_name += '.dscalar.nii'
        # if os.path.exists(file_name):
        #     print(f'Already extracted {file_name}, skipping...')
        # else:
        print(f"-- Start extracting rsFC {s} {ses_id} {target+type} "
                f"smooth={smooth} kernel={kernel} binarize={thres} keeptop={keeptop} --")
        try:
            start = time.perf_counter()
            # Get the subject's data
            if smooth is None:
                data_cortex_subj, info = dset.get_data(
                    space='fs32k', ses_id=ses_id, type='Tseries', subj=[i])
                data_cortex_subj = data_cortex_subj.squeeze()
            else:
                _, info = dset.get_data(space='fs32k', ses_id=ses_id,
                                        type='Tseries', subj=[i])
                tmp_subjlist = pd.DataFrame({'participant_id': [s]})
                tmp_subjlist.to_csv(f'{dest_dir}/{s}.tsv', sep='\t', index=False)
                data_cortex_subj = smooth_hcp_fs32k(f'{dest_dir}/{s}.tsv', ses_id=ses_id, 
                                                    type='Tseries', smooth=smooth, 
                                                    kernel=kernel, return_data_only=True)
                os.remove(f'{dest_dir}/{s}.tsv')

            if target[:3] == 'Net' or target[:3] == 'Fus':
                names = [f'Network_{i}' for i in range(1, int(res)+1)]
                if target[:3] == 'Fus':
                    icos = [atlas_dir + f'/tpl-fs32k/Icosahedron1002.L.label.gii',
                        atlas_dir + f'/tpl-fs32k/Icosahedron1002.R.label.gii']
                    # Average the subject's cortical data within each icosahedron 
                    # if using the Fusion connectivity maps (which are given at 
                    # icosahedron1002 resolution)
                    data_cortex_subj, _ = conn.average_within_Icos(
                        icos, data_cortex_subj)
                    names = net.header.get_axis(0).name.tolist()
                
                # Regress each network into the fs32k cortical 
                # data to get a run-specific network timecourse
                network_timecourse = conn.regress_networks(
                    net.get_fdata(), data_cortex_subj)
                        
            elif target[:3] == 'Ico':
                # Average within each parcel
                network_timecourse, names = conn.average_within_Icos(
                    net, data_cortex_subj)
                network_timecourse = network_timecourse.T
                sides = np.repeat(['L', 'R'], len(names) / 2)
                names = [f'Ico_{sides[i]}{name}' for i,name in enumerate(names)]  
                        
            elif target[:3] == 'ICA':
                assert smooth is None, "ICA components were calculated on unsmoothed data, \
                                cannot apply functional connectvity profile on smoothed data!"
                ses_dic = {'ses-rest1': 0, 'ses-rest2': 1}
                ica_dir = dset.base_dir + f'/derivatives/group/node_timeseries/3T_HCP1200_MSMAll_d{res}_ts2'
                network_timecourse = np.loadtxt(ica_dir + f'/{s}.txt').T
                network_timecourse = np.hsplit(network_timecourse,2)[ses_dic[ses_id]]
                names = [f'Network_{i}' for i in range(1, int(res)+1)]

            # Calculate the connectivity fingerprint
            coef = conn.connectivity_fingerprint(data_cortex_subj, network_timecourse, info,
                                                type, threshold=thres, keeptop=keeptop)
            # Make info
            runs = np.repeat([info.run.unique()], len(names))
            net_id = np.tile(np.arange(len(names)),
                            int(coef.shape[0] / len(names))) + 1
            info = pd.DataFrame({'sn': [s] * coef.shape[0],
                                'sess': [ses_id] * coef.shape[0],
                                'run': runs,
                                'half': 2 - (runs < runs[-1]),
                                'net_id': net_id,
                                'names': names * int(coef.shape[0] / len(names))})

            # Save the data
            C = atlas.data_to_cifti(coef, info.names)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)

            nb.save(C,  file_name)
            info.to_csv(f'{dest_dir}/{s}_{ses_id}_info-{target+type}.tsv', 
                        sep='\t', index=False)
            
            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - time {elapse}")
        except:
            print(f"   Incomplete - {s} {ses_id} {target+type}!")

def smooth_hcp_fs32k(bulk, ses_id='ses-s1', type='Tseries', smooth=1, kernel=None,
                     return_data_only=True):
    hcp_dataset = ds.DataSetUkbResting(hcp_dir)
    T = pd.read_csv(bulk, sep='\t')

    # get the surfaces for smoothing
    surf_L = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_midthickness.surf.gii'
    surf_R = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_midthickness.surf.gii'

    for s in T.participant_id.drop_duplicates():
        # Make smoothed file name
        dest_dir = hcp_dataset.data_dir.format(s)
        cifti_out = dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_desc' \
                    f'-sm{smooth}{kernel if kernel is not None else ""}.dscalar.nii'
        
        if os.path.exists(cifti_out):
            print(f"Already smoothed for {s} fs32k {ses_id} {smooth}")
            if return_data_only:
                data = nb.load(cifti_out).get_fdata()
                return data
        else:
            print(f'- Smoothing data for {s} fs32k {ses_id} in {smooth}mm ...')
            contain_nan = False
            # Load the unsmoothed data
            input_file = hcp_dataset.data_dir.format(s) \
                        + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
            C = nb.load(input_file)

            # fill nan with zeros if unsmoothed data contains any
            if np.isnan(C.get_fdata()).any():
                contain_nan = True
                mask = np.isnan(C.get_fdata())
                C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
                nb.save(C, f'{s}_tmp.dscalar.nii')
                input_file = f'{s}_tmp.dscalar.nii'

            # Write in smoothed surface data (filled with 0)
            smooth_cmd = f"wb_command -cifti-smoothing {input_file} " \
                        f"{smooth} {smooth} COLUMN {cifti_out} " \
                        f"{f'-{kernel} ' if kernel is not None else ''}" \
                        f"-left-surface {surf_L} -right-surface {surf_R} " \
                        f"-fix-zeros-surface"
            subprocess.run(smooth_cmd, shell=True)
            
            C = nb.load(cifti_out)
            data = C.get_fdata()
            
            if contain_nan:
                os.remove(f'{s}_tmp.dscalar.nii')
                # Replace 0s back to NaN (we don't want the 0s impact model learning)
                data[mask] = np.nan
                C = nb.Cifti2Image(dataobj=data, header=C.header)
                nb.save(C, cifti_out)
            
            if return_data_only:
                os.remove(cifti_out)
                return data

def binarize_rsfc_fs32k(bulk, ses_id='ses-s1', type='Tseries', smooth=None, kernel=None,
                        thres=None):
    hcp_dataset = ds.DataSetUkbResting(hcp_dir)
    T = pd.read_csv(bulk, sep='\t')

    for s in T.participant_id:
        # Make smoothed file name
        dest_dir = hcp_dataset.data_dir.format(s)
        
        if smooth is not None:
            file_name = f'{dest_dir}/{s}_space-fs32k_{ses_id}'\
                        f'_{type}_desc-sm{smooth}'\
                        f'{kernel if kernel is not None else ""}'
        else:
            file_name = f'{dest_dir}/{s}_space-fs32k_{ses_id}'\
                        f'_{type}'

        print(f'- Banarizing data for {s} fs32k {ses_id} {smooth} ...')
        # Load the data to be binarized
        C = nb.load(file_name + '.dscalar.nii')
        coef = conn.binarize_top_percent(C.get_fdata(), percent=thres, keep_top=False)

        C = nb.Cifti2Image(dataobj=coef, header=C.header)
        nb.save(C, file_name + f'_binarized-{thres}.dscalar.nii')
        
        print(f'- Done.')

if __name__ == "__main__":
    # make_info(type='Tseries', ses_id='ses-rest1')
    # make_info(type='Tseries', ses_id='ses-rest2')
    #  -- Extract timeseries --
    # extract_hcp_timeseries(
    #     ses_id='ses-rest1', type='Tseries', atlas='MNISymC2')
    # extract_hcp_timeseries(
    #     ses_id='ses-rest2', type='Tseries', atlas='MNISymC2')

    # extract_hcp_timeseries(ses_id='ses-rest1', type='Tseries', atlas='fs32k')
    # extract_hcp_timeseries(ses_id='ses-rest2', type='Tseries', atlas='fs32k')

    # -- Get connectivity fingerprint --
    dname = 'HCP'
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest1')
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest2')

    for t in [162]:
        get_hcp_fs32k_rsfc(type=f'Ico{t}Run', space='fs32k', ses_id='ses-rest1', 
                        subj_list='/subj_list/HCP40_training_set.tsv',
                        smooth=None, kernel='fwhm', thres=None, keeptop=False)
        get_hcp_fs32k_rsfc(type=f'Ico{t}Run', space='fs32k', ses_id='ses-rest2',
                        subj_list='/subj_list/test_split/HCP923_test_set_split_1.tsv',
                        smooth=4, kernel='fwhm', thres=0.1, keeptop=False)
        # get_hcp_fs32k_rsfc(type='Ico42Run', space='fs32k', ses_id='ses-rest1',
        #                    subj_list='/subj_list/HCP40_validation_set.tsv',
        #                    smooth=4, kernel='fwhm')
        # get_hcp_fs32k_rsfc(type='Ico42Run', space='fs32k', ses_id='ses-rest2',
        #                    subj_list='/subj_list/HCP40_validation_set.tsv',
        #                    smooth=4, kernel='fwhm')
        # get_hcp_fs32k_rsfc(type='Ico42Run', space='fs32k', ses_id='ses-rest1', 
        #                    subj_list='/subj_list/HCP923_test_set.tsv',
        #                    smooth=4, kernel='fwhm')
        # get_hcp_fs32k_rsfc(type='Ico42Run', space='fs32k', ses_id='ses-rest2', 
        #                    subj_list='/subj_list/HCP923_test_set.tsv',
        #                    smooth=4, kernel='fwhm')

    #  -- fs32k smoothing (cortex)
    # smooth_hcp_fs32k(hcp_dir + '/subj_list/HCP40_validation_set.tsv', ses_id='ses-rest1',
    #                 type=f'Tseries', smooth=4, kernel='fwhm', return_data_only=False)
    # smooth_hcp_fs32k(hcp_dir + '/subj_list/HCP40_validation_set.tsv', ses_id='ses-rest2',
    #                 type=f'Tseries', smooth=4, kernel='fwhm', return_data_only=False)
    
    #  -- fs32k binarizing (cortex)
    for t in [0.25, 0.5, 0.75]:
        binarize_rsfc_fs32k(hcp_dir + '/subj_list/HCP40_training_set.tsv', ses_id='ses-rest1',
                            type='ICA50Run', smooth=4, kernel='fwhm', thres=t)
        binarize_rsfc_fs32k(hcp_dir + '/subj_list/HCP40_training_set.tsv', ses_id='ses-rest2',
                            type='ICA50Run', smooth=4, kernel='fwhm', thres=t)
    name = 'HCP'
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='fs32k', ses_id='ses-rest1')
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='fs32k', ses_id='ses-rest2')
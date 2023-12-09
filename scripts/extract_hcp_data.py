# Script for getting all the HCP data for cerebellar-cortical connectivity
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetHcpResting
import nibabel as nb
import SUITPy as suit
import os
import sys
import matplotlib.pyplot as plt
from ProbabilisticParcellation.util import plot_multi_flat, plot_data_flat
import re

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

hcp_dir = base_dir + '/HCP'
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


def extract_connectivity_fingerprint(type='Net69Run', space='MNISymC3', ses_id='ses-rest1'):
    """Extracts the connectivity fingerprint for each network in the HCP data
    Steps:  Step 1: Regress each network into the fs32k cortical data to get a run-specific network timecourse
            Step 2: Get the correlation of each voxel with each network timecourse (connectivity fingerprint)
            Step 3: Save the data.
    """

    hcp_dataset = DataSetHcpResting(hcp_dir)

    # Load the networks
    target, type = re.findall('[A-Z][^A-Z]*', type)
    net = nb.load(hcp_dataset.base_dir +
                  f'/targets/{target}_space-fs32k.dscalar.nii')

    atlas, _ = am.get_atlas(space, hcp_dataset.atlas_dir)

    T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    for p, participant_id in enumerate(T.participant_id):
        # Get cortical data
        data_cortex, _ = hcp_dataset.get_data(
            space='fs32k', ses_id=ses_id, type='Tseries', subj=[p])

        # Regress each network into the fs32k cortical data to get a run-specific network timecourse
        network_timecourse = hcp_dataset.regress_networks(
            net.get_fdata(), data_cortex)

        # Calculate the connectivity fingerprint
        data_cereb, info = hcp_dataset.get_data(
            space=space, ses_id=ses_id, type='Tseries', subj=[p])
        data_cereb = data_cereb.squeeze()

        coef = hcp_dataset.connectivity_fingerprint(
            data_cereb, network_timecourse, info, type)
        # Make info
        names = [f'Network_{i}' for i in range(1, 70)]
        runs = np.repeat([info.run.unique()], len(names))
        net_id = np.tile(np.arange(len(names)),
                         int(coef.shape[0] / len(names))) + 1
        info = pd.DataFrame({'sn': [participant_id] * coef.shape[0],
                             'sess': [ses_id] * coef.shape[0],
                             'run': runs,
                             'half': 2 - (runs < runs[-1]),
                             'net_id': net_id,
                             'names': names * int(coef.shape[0] / len(names))})

        # Save the data

        C = atlas.data_to_cifti(coef, info.names)
        dest_dir = hcp_dataset.base_dir + \
            f'/derivatives/{participant_id}/data/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        nb.save(C, dest_dir +
                f'{participant_id}_space-{space}_{ses_id}_{target+type}.dscalar.nii')
        info.to_csv(
            dest_dir + f'{participant_id}_{ses_id}_info-{target+type}.tsv', sep='\t', index=False)


def extract_connectivity_fingerprint_da(type='Ico162Run', space='MNISymC3', ses_id='ses-rest1'):
    """Extracts the connectivity fingerprint for each network in the HCP data

    Args:
        type: data extraction type, 'IcoXXXRun', 'IcoXXXAll', 'NetXXXRun', or
              'NetXXXRun', where XXX indicates the number of networks
        space: the space of cerebellar time series
        ses_id: session ID

    Returns:
        Write in the extracted data to CIFTI format along with its .tsv info file

    Steps:  Step 1: Regress each network into the fs32k cortical data to get
                    a run-specific network timecourse
            Step 2: Get the correlation of each voxel with each network
                    timecourse (connectivity fingerprint)
            Step 3: Save the data.
    """

    hcp_dataset = DataSetHcpResting(hcp_dir)

    # Load the networks
    target, type = re.findall('[A-Z][^A-Z]*', type)

    # 1. Extract connectivity from ICA Network
    if target.startswith('Net'):
        net = nb.load(hcp_dataset.base_dir +
                      f'/targets/tpl-fs32k_{target}.dscalar.nii')
        names = [f'Network_{i}' for i in range(1, net.shape[0] + 1)]

    # 2. Extract connectivity from Icosahedrons
    elif target.startswith('Ico'):
        res = ''.join(re.findall('\d+', target))
        # Get cortical parcelation
        labels, masks = [], []
        for i, h in enumerate(['L', 'R']):
            dir = atlas_dir + '/tpl-fs32k'
            labels += [dir + f'/Icosahedron-{res}_Sym.32k.{h}.label.gii']
            masks += [dir + f'/tpl-fs32k_hemi-{h}_mask.label.gii']

        surf_parcel = am.AtlasSurface(
            'Coretex', masks, ['cortex_left', 'cortex_right'])

        net = surf_parcel.get_parcel(labels, None)[0]
        bpa = surf_parcel.get_parcel_axis()
        names = list(bpa.name)

    atlas, _ = am.get_atlas(space, hcp_dataset.atlas_dir)

    T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    for p, participant_id in enumerate(T.participant_id):
        print(
            f'-Extracting sub {participant_id} using Network: {target}, Type: {type} ...')
        # Get cortical data
        data_cortex, _ = hcp_dataset.get_data(
            space='fs32k', ses_id=ses_id, type='Tseries', subj=[p])

        if target.startswith('Net'):
            # Regress each network into the fs32k cortical data to get a run-specific network timecourse
            network_timecourse = hcp_dataset.regress_networks(
                net.get_fdata(), data_cortex)
        elif target.startswith('Ico'):
            # Average
            network_timecourse = hcp_dataset.average_within_Icos(
                net - 1, data_cortex.squeeze())

        # Calculate the connectivity fingerprint
        data_cereb, info = hcp_dataset.get_data(
            space=space, ses_id=ses_id, type='Tseries', subj=[p])
        data_cereb = data_cereb.squeeze()
        coef = hcp_dataset.connectivity_fingerprint(
            data_cereb, network_timecourse, info, type)

        # Make info
        runs = np.repeat([info.run.unique()], len(names))
        net_id = np.tile(np.arange(len(names)),
                         int(coef.shape[0] / len(names))) + 1
        info = pd.DataFrame({'sn': [participant_id] * coef.shape[0],
                             'sess': [ses_id] * coef.shape[0],
                             'run': runs,
                             'half': 2 - (runs < runs[-1]),
                             'net_id': net_id,
                             'names': names * int(coef.shape[0] / len(names))})

        # Save the data

        C = atlas.data_to_cifti(coef, info.names)
        dest_dir = hcp_dataset.base_dir + \
            f'/derivatives/{participant_id}/data/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        nb.save(C, dest_dir +
                f'{participant_id}_space-{space}_{ses_id}_{target+type}.dscalar.nii')
        info.to_csv(
            dest_dir + f'{participant_id}_{ses_id}_info-{target+type}.tsv', sep='\t', index=False)


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
    extract_connectivity_fingerprint(
        type='Net69Run', space='MNISymC2', ses_id='ses-rest1')
    extract_connectivity_fingerprint(
        type='Net69Run', space='MNISymC2', ses_id='ses-rest2')

    # extract_hcp_timeseries(
    #     ses_id='ses-rest1', type='Tseries', atlas='fs32k')

    # -- Group average --
    # group_average_hcp(type='Ico162Run', atlas='MNISymC3')
    # pass
    # extract_hcp_timeseries(
    #     ses_id='ses-rest1', type='Tseries', atlas='fs32k')
    # pass

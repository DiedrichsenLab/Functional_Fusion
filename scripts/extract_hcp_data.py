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


def get_connectivity(type='Net69Run', space='MNISymC3', ses_id='ses-rest1'):
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
        names = [axis[0] for axis in net.header.get_axis(0)]
        runs = np.repeat([info.run.unique()], len(names))
        net_id = np.tile(np.arange(len(names)),
                         int(coef.shape[0] / len(names)))
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
        fname = f'{participant_id}_{ses_id}_{target+type}'
        nb.save(C, dest_dir + fname + '.dscalar.nii')
        info.to_csv(
            dest_dir + fname + '.tsv', sep='\t', index=False)


if __name__ == "__main__":
    #  -- Extract timeseries --
    # extract_hcp_timeseries(
    #     ses_id='ses-rest1', type='Tseries', atlas='MNISymC3')
    # extract_hcp_timeseries(
    #     ses_id='ses-rest2', type='Tseries', atlas='MNISymC3')
    # extract_hcp_timeseries(ses_id='ses-rest1', type='Tseries', atlas='fs32k')
    # extract_hcp_timeseries(ses_id='ses-rest2', type='Tseries', atlas='fs32k')

    # -- Get connectivity fingerprint --
    # Step 1: Regress each network into the fs32k cortical data to get a run-specific network timecours
    get_connectivity()
    # Step 2: Get the correlation of each voxel with each network timecourse
    pass

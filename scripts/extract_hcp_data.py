# Script for getting all the HCP data for cerebellar-cortical connectivity
import pandas as pd
import shutil
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
from ProbabilisticParcellation.util import plot_multi_flat, plot_data_flat
import re
import Functional_Fusion.connectivity as conn

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
    dname = 'HCP'
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest1')
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest2')
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='SUIT3', ses_id='ses-rest', subj=subject_subset)


    conn.get_connectivity_fingerprint(dname,
                                      type='Fus06All', space='MNISymC3', ses_id='ses-rest1')
    conn.get_connectivity_fingerprint(dname,
                                      type='Fus06All', space='MNISymC3', ses_id='ses-rest2')

    # extract_hcp_timeseries(
    #     ses_id='ses-rest1', type='Tseries', atlas='fs32k')

    # -- Group average --
    # group_average_hcp(type='Ico162Run', atlas='MNISymC3')
    # pass
    # extract_hcp_timeseries(
    #     ses_id='ses-rest1', type='Tseries', atlas='fs32k')
    # pass

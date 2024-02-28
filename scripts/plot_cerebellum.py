
from pathlib import Path
import numpy as np
import pandas as pd
import os
import os.path as op
import sys

import Functional_Fusion.util as util
import Functional_Fusion.matrix as matrix
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetHcpResting
import scipy.linalg as sl
import nibabel as nb
import nitools as nt
import SUITPy as suit
import matplotlib.pyplot as plt

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')



def plot_cerebellum(dset, subject='group', sessions=None, atlas='SUIT3', type=None, cond='all', savefig=False, cmap='hot', colorbar=False):
    """Loads group data in SUIT3 space from a standard experiment structure
    averaged across all subjects and projects to SUIT flatmap. Saves the results as .png figures in the data/group/figures directory.

    Args:
        dset (dataset object): The dataset object containing the data.
        sub (str, optional): Subject string. Defaults to group to plot data averaged across all subjects.
        session (str, optional): Session string. Defaults to first session of in session list of dataset.
        atlas (str, optional): Atlas string. Defaults to 'SUIT3'.
        type (str, optional): Type - defined in ger_data. Defaults to 'CondHalf'.
        cond (str or list): List of condition indices (e.g. [0,1,2] for the first three conditions) or 'all'. Defaults to 'all'.
        savefig (str, optional): Boolean indicating whether figure should be saved. Defaults to 'False'.
        cmap (str, optional): Matplotlib colour map. Defaults to 'hot'.
        colorbar (str, optional): Boolean indicating whether colourbar should be plotted in figure. Defaults to 'False'.
    """
    if sessions is None:
        sessions = dset.sessions
    if type is None:
        type = dset.default_type
    if subject == 'all':
        subjects = dset.get_participants().participant_id.tolist()
    else:
        subjects = [subject]

    atlasmap, atlasinfo = am.get_atlas(atlas, dset.atlas_dir)

    for sub in subjects:
        print(f'Plotting {sub}')
        for session in sessions:
            info = dset.data_dir.split(
                '/{0}')[0] + f'/{sub}/data/{sub}_{session}_info-{type}.tsv'
            data = dset.data_dir.split(
                '/{0}')[0] + f'/{sub}/data/{sub}_space-{atlas}_{session}_{type}.dscalar.nii'

            # Load average
            C = nb.load(data)
            D = pd.read_csv(info, sep='\t')
            X = C.get_fdata()
            # limes = [X[np.where(~np.isnan(X))].min(), X[np.where(~np.isnan(X))].max()] # cannot use nanmax or nanmin because memmap does not have this attribute
            limes = [np.percentile(X[np.where(~np.isnan(X))], 5), np.percentile(
                X[np.where(~np.isnan(X))], 95)]

            if session == 'ses-rest':
                cond_col = 'names'
            else:
                cond_col = 'cond_names'

            if cond == 'all':
                conditions = D[cond_col]
            else:
                conditions = D[cond_col]

            # -- Plot each condition in seperate figures --
            dest_dir = dset.data_dir.split('/{0}')[0] + f'/{sub}/figures/'
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            for c in conditions:
                condition_name = c.strip()
                D[D[cond_col] == c].index
                Nifti = atlasmap.data_to_nifti(
                    X[D[D[cond_col] == c].index, :].mean(axis=0))
                surf_data = suit.flatmap.vol_to_surf(
                    Nifti, space=atlasinfo['normspace'])
                fig = suit.flatmap.plot(
                    surf_data, render='matplotlib', new_figure=True, cscale=limes, cmap=cmap, colorbar=colorbar)
                fig.set_title(condition_name)

                # save figure
                if savefig:
                    plt.savefig(
                        dest_dir + f'{sub}_{session}_{condition_name}.png')
                plt.clf()


if __name__ == "__main__":
    type='Net69Run'
    atlas='MNISymC3'

    hcp_dir = base_dir + '/HCP'
    atlas_dir = base_dir + '/Atlases'
    hem_name = ['cortex_left', 'cortex_right']
    
    hcp_dataset = DataSetHcpResting(hcp_dir)
    hcp_dataset.group_average_data(
        ses_id='ses-rest1', type=type, atlas=atlas)
    hcp_dataset.plot_cerebellum(subject='group', sessions=[
                                    'ses-rest1', 'ses-rest2'], type=type, atlas=atlas, savefig=True, colorbar=True)
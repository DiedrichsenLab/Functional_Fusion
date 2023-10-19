#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Individual parcellation using localizing data

Created on 10/16/2023 at 2:02 PM
Author: dzhi
"""
import numpy as np
import torch as pt
import nibabel as nb
import nitools as nt
import matplotlib.pyplot as plt
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import HierarchBayesParcel.arrangements as ar
import HierarchBayesParcel.full_model as fm
from FusionModel.util import plot_data_flat

# pytorch cuda global flag
# pt.cuda.is_available = lambda : False
pt.set_default_tensor_type(pt.cuda.FloatTensor
                           if pt.cuda.is_available() else
                           pt.FloatTensor)

def plot_multi_flat(data, atlas, grid, cmap='tab20b', dtype='label',
                    cscale=None, titles=None, colorbar=False,
                    save_fig=False):
    """ Plot multiple flatmaps in a grid

    Args:
        data: the input parcellations, shape(N, K, P) where N indicates
              the number of parcellations, K indicates the number of
              parcels, and P is the number of vertices.
        atlas: the atlas name used to plot the flatmap
        grid: the grid shape of the subplots
        cmap: the colormap used to plot the flatmap
        dtype: the data type of the input data, 'label' or 'prob'
        cscale: the color scale used to plot the flatmap
        titles: the titles of the subplots
        colorbar: whether to plot the colorbar
        save_fig: whether to save the figure, default format is png

    Returns:
        The plt figure plot
    """

    if isinstance(data, np.ndarray):
        n_subplots = data.shape[0]
    elif isinstance(data, list):
        n_subplots = len(data)

    if not isinstance(cmap, list):
        cmap = [cmap] * n_subplots

    for i in np.arange(n_subplots):
        plt.subplot(grid[0], grid[1], i + 1)
        plot_data_flat(data[i], atlas,
                       cmap=cmap[i],
                       dtype=dtype,
                       cscale=None,
                       render='matplotlib',
                       colorbar=(i == 0) & colorbar)

        plt.title(titles[i])
        plt.tight_layout()

    if save_fig:
        plt.savefig('/indiv_parcellations.png')

if __name__ == "__main__":
    K=17
    # Step 1.1: Load the atlas
    atlas, _ = am.get_atlas('MNISymC3', 'Y:\data\FunctionalFusion\Atlases')

    # Step 1.2: Load the group prior from a pre-trained model
    model_dir = 'Y:/data/Cerebellum/ProbabilisticParcellationModel/Models'
    model_name = f'/Models_03/asym_Md_space-MNISymC3_K-{K}'
    fname = model_dir + model_name
    U, _ = ar.load_group_parcellation(fname, device='cuda')

    # Step 1.3: Build the arrangement model
    ar_model = ar.build_arrangement_model(U, K, atlas)

    # Step 2.1: Load the individual localizing data / info
    tdata, cond_v, part_v, sub_ind = ds.build_dataset_from_fusionProject('MDTB', atlas,
                                            'Y:/data/FunctionalFusion',
                                            sess='all', cond_ind='cond_num_uni',
                                            type='CondHalf', part_ind='half',
                                            subj=None, join_sess=False,
                                            join_sess_part=False, smooth=None)

    # Step 2.2: Compute the individual parcellations
    indiv_par, _ = fm.get_indiv_parcellation(ar_model, tdata, atlas, cond_v, part_v,
                                             sub_ind, return_soft_parcel=True)

    # Step 3: Save the individual parcellations as a nifti/gifti file
    # Step 3.1: Convert the individual parcellations to gifti file
    # gii_file = nt.make_label_gifti(indiv_par.cpu().numpy().transpose(),
    #                                labels=["label_{}".format(i) for i in range(K)],
    #                                column_names=["subj_{}".format(i+1)
    #                                              for i in range(indiv_par.shape[0])])
    # nb.save(gii_file, '/Md_Asym_17.dlabel.gii')
    # TODO: here we need to write a function to convert the
    #       individual parcellations to nifti file

    # Step 4: Visualization
    # Step 4.1: plot the group parcellations for comparison
    plt.figure(figsize=(20,20))
    plot_multi_flat(indiv_par.cpu().numpy(), 'MNISymC3', grid=(6, 4),
                    cmap='tab20', dtype='prob',
                    titles=["subj_{}".format(i+1) for i in range(indiv_par.shape[0])])
    plt.show()
    # Step 4.2: plot the individual parcellations







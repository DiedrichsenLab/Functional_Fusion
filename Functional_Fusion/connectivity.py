import numpy as np
from numpy.linalg import inv, pinv
import nibabel as nb
import Functional_Fusion.util as ut
import Functional_Fusion.atlas_map as am
import pandas as pd
from pathlib import Path
import re
import Functional_Fusion.dataset as ds
import scripts.paths as pt


base_dir = pt.set_base_dir()
atlas_dir = pt.set_atlas_dir(base_dir)


def regress_networks(X, Y):
    """Regresses a spatial map (X) into data (Y).
    Returns the network timecourses.

    Args:
        X (np.arry): 4D Network data of the signal components
            (default input networks are in fs32k Space: 59518 vertices x nComponents )
        Y (<nibabel CIFTI image object>): fMRI timeseries in volume
            Has to be in the same space as networks (59518 vertices x nTimepoints )
    Returns:
        network_timecourse (np.ndarray):
            A numpy array (nTimepoints x nNetworks) with the fMRI timecourse for
            each resting-state network
    """
    X = X.T
    Y = Y.T.squeeze()
    d_excluded = np.where(np.isnan(Y))[0].shape[0]
    v_excluded = np.unique(np.where(np.isnan(Y))[0]).shape[0]
    if not v_excluded == 0:
        print(
            f'Setting nan datapoints ({v_excluded} unique vertices) to zero. Entire timeseries: {d_excluded/v_excluded == Y.shape[1]}')
    Y[np.isnan(Y)] = 0
    network_timecourse = np.matmul(np.linalg.pinv(X), Y)

    return network_timecourse


def average_within_Icos(label_file, data, atlas="fs32k"):
    """Average the raw time course for voxels within a parcel

    Args:
        label_file (str): cortical parcellation label file
        Y (np.ndarray): fMRI timeseries in volume
            Has to be in the same space as networks (59518 vertices x nTimepoints)
    Returns:
        A numpy array (nNetworks x nTimepoints) with the fMRI timecourse for
        each resting-state network
    """

    # create an instance of atlas to get the label vector
    atlas, ainfo = am.get_atlas(atlas)

    # create label_vector by passing on the label file
    # Set unite_struct to true if you want to integrate over left and right hemi
    atlas.get_parcel(label_file, unite_struct=False)

    # use agg_parcel to aggregate data over parcels and get the list of unique parcels
    parcel_data, parcels = agg_parcels(
        data, atlas.label_vector, fcn=np.nanmean)

    # fill nan value in Y to zero
    print("Setting nan datapoints (%d unique vertices) to zero"
          % np.unique(np.where(np.isnan(parcel_data))[1]).shape[0])
    # Y = np.nan_to_num(np.transpose(Y))
    parcel_data = np.nan_to_num(parcel_data)

    # return np.matmul(np.linalg.pinv(indicator_mat.T), Y)
    return parcel_data, parcels


def connectivity_fingerprint(source, target, info, type):
    """ Calculate the connectivity fingerprint of a target region

    Args:
        source (np.ndarray): Source data
        target (np.nzdarray): Target timecourse
        info (pandas.DataFrame): Information about the source data
        type (str): Type of fingerprint to calculate ('Run' or 'All').
                    Estimates fingerprint from each run seperately or from all concatenated runs.

    Returns:
        coef (np.ndarray): Connectivity fingerprint
    """
    coefs = []
    if type == 'Run':
        for run in info.run.unique():
            data_run = source[info.run == run]
            net_run = target.T[info.run == run]
            coef = ut.correlate(data_run, net_run)
            coefs.append(coef)

    elif type == 'All':
        coef = ut.correlate(source, target)
        coefs.append(coef)

    return np.vstack(coefs)


def get_connectivity_fingerprint(dname, type='Net69Run', space='MNISymC3', ses_id='ses-rest1', subj=None):
    """Extracts the connectivity fingerprint for each network in the HCP data
    Steps:  Step 1: Regress each network into the fs32k cortical data to get a run-specific network timecourse
            Step 2: Get the correlation of each voxel with each network timecourse (connectivity fingerprint)
            Step 3: Save the data.
    """

    dset = ds.get_dataset_class(base_dir, dname)

    # Load the networks
    target, type = re.findall('[A-Z][^A-Z]*', type)
    net = nb.load(dset.base_dir +
                  f'/../targets/{target}_space-fs32k.dscalar.nii')

    atlas, _ = am.get_atlas(space, dset.atlas_dir)

    T = pd.read_csv(dset.base_dir + '/participants.tsv', sep='\t')

    # Deal with subset of subjects
    if subj is not None:
        subj = [T.participant_id.tolist().index(s) for s in subj]
        T = T.iloc[subj]

        
    for p, row in T.iterrows():
        participant_id = row.participant_id

        # Get cortical data
        data_cortex, _ = dset.get_data(
            space='fs32k', ses_id=ses_id, type='Tseries', subj=[p])

        # Regress each network into the fs32k cortical data to get a run-specific network timecourse
        network_timecourse = regress_networks(
            net.get_fdata(), data_cortex)

        # Calculate the connectivity fingerprint
        data_cereb, info = dset.get_data(
            space=space, ses_id=ses_id, type='Tseries', subj=[p])
        data_cereb = data_cereb.squeeze()

        coef = connectivity_fingerprint(
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
        dest_dir = dset.base_dir + \
            f'/derivatives/{participant_id}/data/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        nb.save(C, dest_dir +
                f'{participant_id}_space-{space}_{ses_id}_{target+type}.dscalar.nii')
        info.to_csv(
            dest_dir + f'{participant_id}_{ses_id}_info-{target+type}.tsv', sep='\t', index=False)


    # orig = nb.load(base_dir +
    #               f'/targets/{target}')
    # # extract left hemisphere and right hemisphere from the cifti file
    # lh, rh = orig.dataobj[..., :59412], orig.dataobj[..., 59412:]
    # orig.header.get_axis(1).name

    # # Get all data from cifti structure cortex left
    # axis = orig.header.get_axis(1)
    # surf_name = 'CIFTI_STRUCTURE_CORTEX_LEFT'
    # # Get data indices for the structure
    # orig.header.get_axis(1).name == surf_name
    # data_indices = orig.header.get_axis(1).get_index(surf_name)
    # orig[orig.name == 'CIFTI_STRUCTURE_CORTEX_LEFT']
    # # get the data array with all the time points, for surf_name

    # for idx, (name, slc, bm) in enumerate(orig.header.get_axis(1).iter_structures()):
    #     print((str(name), slc))
    
    # data = []
    # assert isinstance(axis, nb.cifti2.BrainModelAxis)


    # for name, data_indices, model in axis.iter_structures():  # Iterates over volumetric and surface structures
    #     if name == surf_name:
    #         data.append(orig.dataobj[..., data_indices])

    # assert isinstance(axis, nb.cifti2.BrainModelAxis)
    # for name, data_indices, model in axis.iter_structures():  # Iterates over volumetric and surface structures
    #     if name == surf_name:
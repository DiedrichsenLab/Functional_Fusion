import numpy as np
from numpy.linalg import inv, pinv
import nibabel as nb
import util as ut
import atlas_map as am

def regress_networks(self, X, Y):
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
    print(
        f'Setting nan datapoints ({v_excluded} unique vertices) to zero. Entire timeseries: {d_excluded/v_excluded == Y.shape[1]}')
    Y[np.isnan(Y)] = 0
    network_timecourse = np.matmul(np.linalg.pinv(X), Y)

    return network_timecourse

def average_within_Icos(self, label_file, data, atlas="fs32k"):
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


def connectivity_fingerprint(self, source, target, info, type):
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
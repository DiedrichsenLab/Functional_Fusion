import numpy as np
from numpy.linalg import inv, pinv
import nibabel as nb
import Functional_Fusion.util as ut
import Functional_Fusion.atlas_map as am
import pandas as pd
from pathlib import Path
import re
import Functional_Fusion.dataset as ds
import paths as paths


base_dir = paths.set_base_dir()
atlas_dir = paths.set_atlas_dir(base_dir)


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
    parcel_data, parcels = ds.agg_parcels(
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
                    Alternatively, average cortical timecourse within each Icosahedron parcel to get network timecourse
            Step 2: Get the correlation of each voxel with each network timecourse (connectivity fingerprint)
            Step 3: Save the data.

    Args:
        dname (str): Name of the dataset
        type (str): Type of fingerprint to calculate (e.g. 'Net69Run' or 'Ico42All').
                    Estimates fingerprint from each run seperately or from all concatenated runs.
                    Net69Run: 69 networks, estimated from each run seperately
                    Net69All: 69 networks, estimated from all concatenated runs
                    Ico42All: 42 Icosahedron parcels, estimated from all concatenated runs
                    Ico162All: 162 Icosahedron parcels, estimated from all concatenated runs
        space (str): Space of the cortical data ('MNISymC2', 'MNISymC3', 'fs32k')
        ses_id (str): Session ID
        subj (list): List of subjects to extract the fingerprint for
    """
    # Load dataset
    dset = ds.get_dataset_class(base_dir, dname)

    T = pd.read_csv(dset.base_dir + '/participants.tsv', sep='\t')

    # Deal with subset of subjects
    if subj is not None:
        subj = [T.participant_id.tolist().index(s) for s in subj]
        T = T.iloc[subj]
    
    # Get cortical data
    data_cortex, _ = dset.get_data(
        space='fs32k', ses_id=ses_id, type='Tseries', subj=subj)

    # Get cerebellar data
    data_cereb, info_cereb = dset.get_data(
            space=space, ses_id=ses_id, type='Tseries', subj=subj)

    # Load the cortical networks
    target, type = re.findall('[A-Z][^A-Z]*', type)        
    res = target[3:]
    if target[:3] == 'Net':
        net = nb.load(dset.base_dir +
                    f'/../targets/{target}_space-fs32k.dscalar.nii')
    elif target[:3] == 'Ico':
        net = [atlas_dir + f'/tpl-fs32k/Icosahedron{res}.L.label.gii',
            atlas_dir + f'/tpl-fs32k/Icosahedron{res}.R.label.gii']

    atlas, _ = am.get_atlas(space, dset.atlas_dir)
        
    for p, row in enumerate(T.itertuples()):
        participant_id = row.participant_id

        # Get the subject's data
        data_cortex_subj = data_cortex[p,:,:]
        data_cereb_subj = data_cereb[p,:,:]

        if target[:3] == 'Net':
            # Regress each network into the fs32k cortical data to get a run-specific network timecourse
            network_timecourse = regress_networks(
                net.get_fdata(), data_cortex_subj)
            names = [f'Network_{i}' for i in range(1, int(res)+1)]
        elif target[:3] == 'Ico':
            # Average within each parcel
            network_timecourse, names = average_within_Icos(
                net, data_cortex_subj)
            network_timecourse = network_timecourse.T
            sides = np.repeat(['L', 'R'], len(names) / 2)
            names = [f'Ico_{sides[i]}{name}' for i,name in enumerate(names)]            

        # Calculate the connectivity fingerprint
        coef = connectivity_fingerprint(
            data_cereb_subj, network_timecourse, info_cereb, type)
        # Make info
        runs = np.repeat([info_cereb.run.unique()], len(names))
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


def get_cortical_target(target):

    orig = nb.load(base_dir +
                  f'/targets/{target}')
    
    lh, rh = ut.surf_from_cifti(orig)
    seed_names = ['Network_{}'.format(i)
                  for i in range(1, len(orig.header.get_axis(0).name) + 1)]
    
    bpa = nb.cifti2.ScalarAxis(seed_names)
    # lh = cifti2.Cifti2Image(lh, transforms.get_cifti2_axes('32k'))
    target_name = target.split('/')[-1].split('_space')[0]
    print(f'Writing {target_name} ...')

    # Remove medial wall
    atlas, _ = am.get_atlas('fs32k', atlas_dir)
    bmc = atlas.get_brain_model_axis()
    lh_masked = [data[atlas.mask[0]] for data in lh]
    rh_masked = [data[atlas.mask[1]] for data in rh]

    # --- Build a connectivity CIFTI-file and save ---
    # Make the object
    header = nb.Cifti2Header.from_axes((bpa, bmc))
    cifti_img = nb.Cifti2Image(
        dataobj=np.c_[lh_masked, rh_masked], header=header)
    dest_dir = base_dir + '/targets/'
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    nb.save(cifti_img, dest_dir + f'/{target_name}_space-fs32k.dscalar.nii')

        

if __name__ == "__main__":
    # get_cortical_target('orig/hcp_1200/Net300_space-fs32k.dscalar.nii')
    # get_cortical_target('orig/hcp_1200/Net100_space-fs32k.dscalar.nii')
    # get_cortical_target('orig/hcp_1200/Net50_space-fs32k.dscalar.nii')
    get_cortical_target('orig/hcp_1200/Net15_space-fs32k.dscalar.nii')
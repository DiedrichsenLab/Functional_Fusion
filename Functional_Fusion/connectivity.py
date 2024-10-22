import numpy as np
from numpy.linalg import inv, pinv
import nibabel as nb
import Functional_Fusion.util as ut
import Functional_Fusion.atlas_map as am
import pandas as pd
from pathlib import Path
import re
import Functional_Fusion.dataset as ds
# import scripts.fusion_paths as paths


# base_dir = paths.set_base_dir()
base_dir = '/data/tge/Tian/UKBB_full/imaging'
data_dir = base_dir
atlas_dir = base_dir + '/Atlases'


def regress_networks(X, Y):
    """Regresses a spatial map (X) into data (Y).
    Returns the network timecourses.

    Args:
        X (np.array): 4D Network data of the signal components
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


def binarize_top_percent(arr, percent=0.1, keep_top=False):
    """Binarization of the rsFC by giving a top percent
    For example, the top 10% values are 1, the rest are 0s

    Args:
        arr (np.array): the given rsFC matrix to biniarize
        percent (float): the number of top percent to keep

    Returns:
        result (np.int8): the binarized rsFC
    """
    # Ensure the input is a NumPy array / set nan to -1
    arr = np.nan_to_num(np.asarray(arr), nan=-1)
    threshold = np.percentile(arr.flatten(), (1-percent)*100)

    # Apply the threshold to keep the top `percent` values
    if keep_top:
        # 1. keep the original values
        result = np.where(arr >= threshold, arr, 0).astype(np.float32)
    else:
        # 2. set to 1
        result = np.where(arr >= threshold, 1, 0).astype(np.int8)
    
    return result

def keep_top_percent(arr, percent=0.1):
    """Keep the rsFC by giving a top percent
    For example, the top 10% values are kept (the original values),
    the rest are set tp 0s

    Args:
        arr (np.array): the given rsFC matrix
        percent (float): the number of top percent to keep

    Returns:
        result (np.float32): the kept rsFC for top <percent>
    """
    # Ensure the input is a NumPy array
    arr = np.asarray(arr)
    threshold = np.percentile(arr.flatten(), (1-percent)*100)

    # Apply the threshold to keep the top `percent` values
    result = np.where(arr >= threshold, arr, 0).astype(np.float32)
    
    return result

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


def connectivity_fingerprint(source, target, info, type, threshold=None, keeptop=False):
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
            if threshold is not None:
                coef = binarize_top_percent(coef, percent=threshold,
                                            keep_top=keeptop)
<<<<<<< Updated upstream
            coefs.append(coef)
    
    elif type == 'Half':
        for half in info.half.unique():
            data_half = source[info.half == half]
            net_half = target.T[info.half == half]
            coef = ut.correlate(data_half, net_half)
            coefs.append(coef)

    elif type == 'All':
        coef = ut.correlate(source, target.T)
=======
            coefs.append(coef)

    elif type == 'All':
        coef = ut.correlate(source, target)
>>>>>>> Stashed changes
        if threshold is not None:
            coef = binarize_top_percent(coef, percent=threshold,
                                        keep_top=keeptop)
        coefs.append(coef)

    return np.vstack(coefs) 


def get_connectivity_fingerprint(dname, type='Net69Run', space='MNISymC3', ses_id='ses-rest1', smooth=None, subj=None):
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
    dset = ds.get_dataset_class(data_dir, dname)

    T = pd.read_csv(f'{data_dir}/{dname}/participants.tsv', sep='\t')

    # Deal with subset of subjects
    if subj is not None:
        subj = [T.participant_id.tolist().index(s) for s in subj]
        T = T.iloc[subj]

    # Load the cortical networks
    type_parts = re.findall('[A-Z][^A-Z]*', type)        
    if len(type_parts) == 2:
        target, type = type_parts
        tseries_type = ''
    elif len(type_parts) == 3:
        target, tseries_type, type = type_parts
    
    
    if tseries_type == 'Fix':
        load_tseries_type = 'FixTseries'
    elif tseries_type == 'Res':
        load_tseries_type = 'Residuals'
    elif tseries_type == '':
        load_tseries_type = 'Tseries'
    

    res = target[3:]
    
    if target[:3] == 'Net':
        net = nb.load(data_dir +
                    f'/targets/{target}_space-fs32k.dscalar.nii')
    elif target[:3] == 'Ico':
        net = [atlas_dir + f'/tpl-fs32k/Icosahedron{res}.L.label.gii',
            atlas_dir + f'/tpl-fs32k/Icosahedron{res}.R.label.gii']
    elif target[:3] == 'Fus':
        net = nb.load(data_dir +
                    f'/targets/{target}_space-fs32k.pscalar.nii')
    atlas, _ = am.get_atlas(space, dset.atlas_dir)
        
    for p, row in enumerate(T.itertuples()):
        participant_id = row.participant_id

        # Get the subject's data
        # Get cortical data
        data_cortex_subj, _ = dset.get_data(
            space='fs32k', ses_id=ses_id, type=load_tseries_type, subj=[row.Index])


        # Get cerebellar data
        data_cereb_subj, info_cereb = dset.get_data(
                space=space, ses_id=ses_id, type=load_tseries_type, subj=[row.Index])
        data_cereb_subj = data_cereb_subj.squeeze()

        if target[:3] == 'Net' or target[:3] == 'Fus':
            names = [f'Network_{i}' for i in range(1, int(res)+1)]
            if target[:3] == 'Fus':
                icos = [atlas_dir + f'/tpl-fs32k/Icosahedron1002.L.label.gii',
                    atlas_dir + f'/tpl-fs32k/Icosahedron1002.R.label.gii']
                # Average the subject's cortical data within each icosahedron if using the Fusion connectivity maps (which are given at icosahedron1002 resolution)
                data_cortex_subj, _ = average_within_Icos(
                    icos, data_cortex_subj)
                names = net.header.get_axis(0).name.tolist()
            # Regress each network into the fs32k cortical data to get a run-specific network timecourse
            network_timecourse = regress_networks(
                net.get_fdata(), data_cortex_subj)
                    
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
        if type == 'Run':
            info = pd.DataFrame({'sn': [participant_id] * coef.shape[0],
                             'sess': [ses_id] * coef.shape[0],
                             'run': runs,
                             'half': 2 - (runs < runs[-1]),
                             'net_id': net_id,
                             'names': names * int(coef.shape[0] / len(names))})
            info['names'] = [f'{d.names.strip()}-run{d.run}' for i, d in info.iterrows()]
        elif type == 'Half':
            
            info = pd.DataFrame({'sn': [participant_id] * coef.shape[0],
                             'sess': [ses_id] * coef.shape[0],
                             'half': np.repeat([1,2],coef.shape[0]/2),
                             'net_id': net_id,
                             'names': names * int(coef.shape[0] / len(names))})
            info['names'] = [f'{d.names.strip()}-half{d.half}' for i, d in info.iterrows()]
        elif type == 'All':
            info = pd.DataFrame({'sn': [participant_id] * coef.shape[0],
                             'sess': [ses_id] * coef.shape[0],
                             'net_id': net_id,
                             'names': names * int(coef.shape[0] / len(names))})

        # Save the data
        C = atlas.data_to_cifti(coef, info.names)
        dest_dir = dset.data_dir.format(participant_id)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        nb.save(C,  f'{dest_dir}/{participant_id}_space-{space}_{ses_id}_{target+tseries_type+type}.dscalar.nii')
        info.to_csv(
            f'{dest_dir}/{participant_id}_{ses_id}_info-{target+tseries_type+type}.tsv', sep='\t', index=False)


def get_cortical_target(target):

    orig = nb.load(data_dir +
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
    dest_dir = data_dir + '/targets/'
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    nb.save(cifti_img, dest_dir + f'/{target_name}_space-fs32k.dscalar.nii')

def make_net67():
    """Script to remove two networks (first and 11th network) from the 69 network set
    First and 11th network are cerebellar-only networks and do not have substantial cortical weights. Therefore they are removed."""
    
    net = nb.load(data_dir +
                    f'/targets/Net69_space-fs32k.dscalar.nii')
    # Remove the first and 10th network
    net_data = np.delete(net.get_fdata(), [0, 10], axis=0)

    # Add the new networks to header
    bpa = nb.cifti2.ScalarAxis(['Network_{}'.format(i)
                                for i in range(1, net_data.shape[0] + 1)])
    bmc = net.header.get_axis(1)
    header = nb.cifti2.Cifti2Header.from_axes((bpa, bmc))
    cifti_img = nb.cifti2.Cifti2Image(dataobj=net_data, header=header)

    # Save the new network
    dest_dir = data_dir + '/targets/'
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    nb.save(cifti_img, dest_dir + f'/Net67_space-fs32k.dscalar.nii')

    pass
        

if __name__ == "__main__":
    # get_cortical_target('orig/hcp_1200/Net300_space-fs32k.dscalar.nii')
    # get_cortical_target('orig/hcp_1200/Net100_space-fs32k.dscalar.nii')
    # get_cortical_target('orig/hcp_1200/Net50_space-fs32k.dscalar.nii')
    # get_cortical_target('orig/hcp_1200/Net15_space-fs32k.dscalar.nii')
    make_net67()
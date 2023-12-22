


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




def show_hcp_group(ses_id='ses-s1', type='Run', atlas='MNISymC3', cond=0, info_column='names', savefig=True):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    C = nb.load(hcp_dataset.data_dir.split('/{0}')[0] +
                f'/group/group_{ses_id}_space-{atlas}_{type}.dscalar.nii')
    D = pd.read_csv(hcp_dataset.data_dir.split('/{0}')[0] +
                    f'/group/group_{ses_id}_info-{type}.tsv', sep='\t')
    X = C.get_fdata()
    limits = [X.min(), X.max()]
    conditions = D[info_column]
    dest_dir = hcp_dataset.data_dir.split('/{0}')[0] + f'/group/figures/'
    if cond == 'all':
        # -- plot all in one figure --
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        plot_multi_flat(X, atlas,
                        grid=(4, 3),
                        dtype='func',
                        cscale=limits,
                        colorbar=True,
                        titles=conditions)
        if savefig:
            plt.savefig(dest_dir + f'group_{ses_id}_{type}.png')
        plt.clf()

    elif cond == 'separate':
        # -- each in seperate figures --
        for i, c in enumerate(conditions):
            plot_data_flat(X[i, :], atlas,
                           dtype='func',
                           cscale=limits,
                           colorbar=False)
            plt.title(c)
            # save figure
            if savefig:
                plt.savefig(dest_dir + f'group_{ses_id}_{type}_{c}.png')
            plt.clf()


def extract_hcp_data(res=162):
    # Make the atlas object
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum', mask_img=mask)

    # initialize the data set object
    hcp_dataset = DataSetHcpResting(hcp_dir)

    # Get the deformation map from MNI to SUIT
    mni_atlas = atlas_dir + '/tpl-MNI152NLin6AsymC'
    deform = mni_atlas + '/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
    mask = mni_atlas + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
    atlas_map = am.AtlasMapDeform(
        hcp_dataset, suit_atlas, 'group', deform, mask)
    atlas_map.build(smooth=2.0)

    # Get the parcelation
    surf_parcel = []
    for i, h in enumerate(['L', 'R']):
        dir = atlas_dir + '/tpl-fs32k'
        gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
        surf_parcel.append(am.AtlasSurfaceParcel(hem_name[i], gifti))

    T = hcp_dataset.get_participants()
    for s in T.participant_id:
        print(f'Extract {s}')
        coef = hcp_dataset.get_cereb_connectivity(s, atlas_map, surf_parcel)
        # Average across runs
        coef = np.nanmean(coef, axis=0)

        # Build a connectivity CIFTI-file and save
        bmc = suit_atlas.get_brain_model_axis()
        bpa = surf_parcel[0].get_parcel_axis(
        ) + surf_parcel[1].get_parcel_axis()
        header = nb.Cifti2Header.from_axes((bpa, bmc))
        cifti_img = nb.Cifti2Image(dataobj=coef, header=header)
        dest_dir = hcp_dataset.data_dir.format(s)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        nb.save(cifti_img, dest_dir + f'/sub-{s}_tessel-{res}.dpconn.nii')



def parcel_hcp_dpconn_cortex(dpconn_file):
    """
    Args:
        dpconn_file (_type_): _description_
    """
    label_file = [atlas_dir + '/tpl-fs32k/ROI.32k.L.label.gii',
                  atlas_dir + '/tpl-fs32k/ROI.32k.R.label.gii', ]
    cifti_img = []
    mask = [atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii',
            atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii']

    for h in range(2):
        C = nb.load(dpconn_file[h])
        roi_atlas = am.AtlasSurfaceParcel(
            'ROI', label_gii=label_file[1 - h], mask_gii=mask[h])
        R = roi_atlas.agg_data(np.asanyarray(C.dataobj))
        bm_cortex = C.header.get_axis(0)
        names = [name.label for name in roi_atlas.label_gii.labeltable.labels[1:]]
        row_axis = nb.cifti2.ScalarAxis(names)
        this_cifti_img = nb.Cifti2Image(R.T, [row_axis, bm_cortex])
        cifti_img.append(this_cifti_img)

    return cifti_img


def indv_hcp_pscalar(res=162, index=range(0, 100), refix=False):
    hcp_dataset = DataSetHcpResting(hcp_dir)
    T = hcp_dataset.get_participants()
    subjects_id = T.participant_id[index]
    for s in T.participant_id[index]:
        data_dir = hcp_dataset.data_dir.format(s)
        if refix:
            C = parcel_hcp_dpconn(
                data_dir + f'/sub-{s}_tessel-{res}-ReFIX.dpconn.nii')
            nb.save(C, data_dir + f'/sub-{s}_tessel-{res}-ReFIX.pscalar.nii')
        else:
            C = parcel_hcp_dpconn(
                data_dir + f'/sub-{s}_tessel-{res}.dpconn.nii')
            nb.save(C, data_dir + f'/sub-{s}_tessel-{res}.pscalar.nii')

        print(f"-Saved scalar file for subject {s}, ReFIX={refix}")

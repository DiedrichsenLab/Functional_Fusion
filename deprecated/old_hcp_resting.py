

class DataSetHcpResting(DataSetCifti):
    def __init__(self, dir):
        super(DataSetHcpResting, self).__init__(base_dir=dir)
        # self.func_dir = self.base_dir + '/{0}/estimates'
        self.derivative_dir = self.base_dir + '/derivatives'
        self.sessions = ['ses-s1', 'ses-s2']
        self.hem_name = ['cortex_left', 'cortex_right']
        self.default_type = 'NetAutoRun'
        self.cond_ind = 'reg_id'
        self.cond_name = 'region_name'
        self.part_ind = 'half'

    def get_data_fnames(self, participant_id):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
        Returns:
            fnames (list): List of fnames
        """
        dirw = self.derivative_dir + f"/{participant_id}" + "/func"
        fnames = []
        for r in range(4):
            fnames.append(
                f'{dirw}/sub-{participant_id}_run-{r}_space-MSMSulc.dtseries.nii')
        return fnames

    def extract_all_suit(self, ses_id='ses-s1', type='IcoAll', atlas='SUIT3', res=162):
        """ MDTB extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'IcoAll': Single estimate per session, Correlation with cortical icosahedron parcels
                'IcoRun': Seperate estimates per run, Correlation with cortical icosahedron parcels
                'NetAll': Single estimate per session, Correlation with cortico-cerebellar resting-state networks estimated with dimensionality = 25.
                'NetRun': Seperate estimates per run, Correlation with cortico-cerebellar resting-state networks estimated with dimensionality = 25.
                'NetAutoAll': Single estimate per session, Correlation with cortico-cerebellar resting-state networks estimated with automatic dimensionality estimation.
                'NetAutoRun': Seperate estimates per run, Correlation with cortico-cerebellar resting-state networks estimated with automatic dimensionality estimation.

                    Defaults to 'IcoAll'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        suit_atlas, _ = am.get_atlas(atlas, self.atlas_dir)

        # Get the deformation map from MNI to SUIT
        mni_atlas = self.atlas_dir + '/tpl-MNI152NLin6AsymC'
        deform = mni_atlas + '/tpl-MNI152NLin6AsymC_space-SUIT_xfm.nii'
        mask = mni_atlas + '/tpl-MNI152NLin6AsymC_res-2_gmcmask.nii'
        atlas_map = am.AtlasMapDeform(suit_atlas.world, deform, mask)
        atlas_map.build(smooth=2.0)

        # Split type information on capital letters
        type_info = re.findall('[A-Z][^A-Z]*', type)
        surf_parcel = None
        if type_info[0] == 'Ico':
            networks = None
            # Get the parcelation
            for i, h in enumerate(['L', 'R']):
                dir = self.atlas_dir + '/tpl-fs32k'
                gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
                surf_parcel.append(
                    am.AtlasSurfaceParcel(self.hem_name[i], gifti))
            bpa = surf_parcel[0].get_parcel_axis(
            ) + surf_parcel[1].get_parcel_axis()
            seed_names = list(bpa.name)

        elif type_info[0] == 'Net':
            dimensionality = type_info[1].casefold()
            # Get the networks
            networkdir = self.base_dir + f'/group_ica/dim_{dimensionality}/'
            networkimg = nb.load(networkdir +
                                 'melodic_IC.nii.gz')
            networks = networkimg.get_fdata()
            net_selected = pd.read_csv(
                networkdir + 'classified_components.txt', sep=', ', skiprows=[0], skipfooter=1, engine='python', header=None, names=['Network', 'Classification', 'IsNoise'], dtype="category")
            networks = networks[:, :, :,
                                net_selected.Classification == 'Signal']
            seed_names = [
                f'Network-{n+1:02}' for n in np.arange(networks.shape[-1])]

        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}')
            if ses_id == 'ses-s1':
                runs = [0, 1]
            elif ses_id == 'ses-s2':
                runs = [2, 3]
            else:
                raise ValueError('Unknown session id.')

            coef = self.get_cereb_connectivity(
                s, atlas_map, runs=runs, type=type, cortical_atlas_parcels=surf_parcel, networks=networks)

            if type_info[-1] == 'All':  # Average across runs
                coef = np.nanmean(coef, axis=0)

                # Make info structure
                reg_ids = np.arange(len(seed_names)) + 1
                info = pd.DataFrame({'sn': [s] * coef.shape[0],
                                    'sess': [ses_id] * coef.shape[0],
                                     'half': [1] * coef.shape[0],
                                     'reg_id': reg_ids,
                                     'region_name': seed_names,
                                     'names': seed_names})

            elif type_info[-1] == 'Run':  # Concatenate over runs
                coef = np.concatenate(coef, axis=0)

                # Make info structure
                run_ids = np.repeat(runs, int(coef.shape[0] / len(runs)))
                reg_ids = np.tile(np.arange(len(seed_names)), 2) + 1
                names = ["{}_run-{}".format(reg_name, run_id)
                         for reg_name, run_id in zip(list(seed_names) * 2, run_ids)]
                info = pd.DataFrame({'sn': [s] * coef.shape[0],
                                    'sess': [ses_id] * coef.shape[0],
                                     'run': run_ids,
                                     'half': 2 - (run_ids < run_ids[-1]),
                                     'reg_id': reg_ids,
                                     'region_name': list(seed_names) * 2,
                                     'names': names})

                # update brain parcel axis (repeat names)
                # bpa = bpa + bpa

            # --- Save cerebellar data as dscalar CIFTI-file and write info to tsv ---
            C = suit_atlas.data_to_cifti(coef, info.names)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(C, dest_dir +
                    f'/{s}_space-{atlas}_{ses_id}_{type}.dscalar.nii')
            info.to_csv(
                dest_dir + f'/{s}_{ses_id}_info-{type}.tsv', sep='\t', index=False)

    def extract_all_fs32k(self, ses_id='ses-s1', type='IcoAll', res=162):
        """ MDTB extraction of atlasmap locations
        from nii files - and filterting or averaring
        as specified.

        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'IcoAll': Single estimate per session,
                          Correlation with cortical icosahedron parcels
                'IcoRun': Seperate estimates per run,
                          Correlation with cortical icosahedron parcels
                'NetAll': Single estimate per session,
                          Correlation with cortico-cerebellar resting-state networks
                'NetRun': Seperate estimates per run,
                          Correlation with cortico-cerebellar resting-state networks
                Defaults to 'IcoAll'.
            res: the resolution of underlying icosahedron. Default 162
        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row
        """
        # Make the atlas object
        mask_L = self.atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-L_mask.label.gii'
        mask_R = self.atlas_dir + '/tpl-fs32k/tpl-fs32k_hemi-R_mask.label.gii'
        atlas = am.AtlasSurface('fs32k', mask_gii=[mask_L, mask_R],
                                structure=['CORTEX_LEFT', 'CORTEX_RIGHT'])
        cortex_mask = [mask_L, mask_R]
        bmc = atlas.get_brain_model_axis()
        seed_names = []

        if type[0:3] == 'Ico':
            networks = None
            # Get the parcelation
            surf_parcel = []
            for i, h in enumerate(['L', 'R']):
                dir = self.atlas_dir + '/tpl-fs32k'
                gifti = dir + f'/Icosahedron-{res}.32k.{h}.label.gii'
                mask = dir + f'/tpl-fs32k_hemi-{h}_mask.label.gii'
                surf_parcel.append(am.AtlasSurfaceParcel(
                    self.hem_name[i], gifti, mask_gii=mask))
            bpa = surf_parcel[0].get_parcel_axis(
            ) + surf_parcel[1].get_parcel_axis()
            seed_names = list(bpa.name)
        elif type[0:3] == 'Net':
            surf_parcel = None
            # Get the networks
            networkdir = self.base_dir + '/group_ica/dim_25/'
            networkimg = nb.load(networkdir + 'melodic_IC.nii.gz')
            networks = networkimg.get_fdata()
            net_selected = pd.read_csv(networkdir + 'classified_components.txt',
                                       sep=', ', skiprows=[0], skipfooter=1,
                                       engine='python', header=None,
                                       names=['Network',
                                              'Classification', 'IsNoise'],
                                       dtype="category")
            networks = networks[:, :, :,
                                net_selected.Classification == 'Signal']
            seed_names = [
                f'Network-{n+1:02}' for n in np.arange(networks.shape[-1])]
        else:
            raise NameError("type must start with either 'Ico' or 'Net'!")
        # Making cifti2 axis for these network name
        bpa = nb.cifti2.ScalarAxis(seed_names)

        T = self.get_participants()
        for s in T.participant_id:
            print(f'Extract {s}, type {type}')
            if ses_id == 'ses-s1':
                runs = [0, 1]
            elif ses_id == 'ses-s2':
                runs = [2, 3]
            else:
                raise ValueError('Unknown session id.')

            coef = self.get_cortical_connectivity(s, cortex_mask=cortex_mask, runs=runs, type=type,
                                                  cortical_atlas_parcels=surf_parcel,
                                                  networks=networks)

            if type[3:7] == 'All':  # Average across runs
                coef = [np.nanmean(c, axis=0) for c in coef]

                # Make info structure
                info = []
                for i, d in enumerate(coef):
                    reg_ids = np.arange(len(seed_names)) + 1
                    this_info = pd.DataFrame({'sn': [s] * d.shape[0],
                                              'hemis': i + 1,
                                              'sess': [ses_id] * d.shape[0],
                                              'half': [1] * d.shape[0],
                                              'reg_id': reg_ids,
                                              'region_name': seed_names,
                                              'names': seed_names})
                    info.append(this_info)

            elif type[3:7] == 'Run':  # Concatenate over runs
                coef = [np.concatenate(c, axis=0) for c in coef]

                # Make info structure
                info = []
                for i, d in enumerate(coef):
                    run_ids = np.repeat(runs, int(d.shape[0] / len(runs)))
                    reg_ids = np.tile(np.arange(len(seed_names)), 2) + 1
                    names = ["{}_run-{}".format(reg_name, run_id)
                             for reg_name, run_id in zip(list(seed_names) * 2, run_ids)]
                    this_info = pd.DataFrame({'sn': [s] * d.shape[0],
                                              'hemis': i + 1,
                                              'sess': [ses_id] * d.shape[0],
                                              'run': run_ids,
                                              'half': 2 - (run_ids < run_ids[-1]),
                                              'reg_id': reg_ids,
                                              'region_name': list(seed_names) * 2,
                                              'names': names})
                    info.append(this_info)
                # update brain parcel axis (repeat names)
                bpa = bpa + bpa

            info = pd.concat([info[0], info[1]])

            # --- Build a connectivity CIFTI-file and save ---
            print(f'Writing {s}, type {type} ...')
            header = nb.Cifti2Header.from_axes((bpa, bmc))
            cifti_img = nb.Cifti2Image(
                dataobj=np.c_[coef[0], coef[1]], header=header)
            dest_dir = self.data_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            nb.save(cifti_img, dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_'
                                          f'{res}.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_space-fs32k_{ses_id}_info-{type}_{res}.tsv',
                        sep='\t', index=False)

    def extract_ts_volume(self,
                          participant_id,
                          atlas_map,
                          ses_id=[0, 1, 2, 3],
                          type='Run'):
        """ Returns the time series data for an atlas map
                runs=[0,1,2,3]):

        Args:
            participant_id (_type_): _description_
            atlas_map (_type_): _description_
        """
        # get the file name for the cifti time series
        fnames, info = self.get_data_fnames(participant_id)
        ts_volume = []
        for r in runs:
            # load the cifti
            ts_cifti = nb.load(fnames[r])

            # get the ts in volume for subcorticals
            ts_vol = nt.volume_from_cifti(ts_cifti)
            # transfer to suit space applying the deformation
            ts_vol = am.get_data4D(ts_vol, [atlas_map])
            ts_volume.append(ts_vol[0])
        return ts_volume

    def extract_ts_surface(self,
                           participant_id,
                           atlas_parcels,
                           runs=[0, 1, 2, 3]):
        """Returns the information from the CIFTI file
        in the 32K surface for left and right hemisphere.

        Args:
            participant_id (_type_): _description_
            atlas_parcel (_type_): _description_
        """
        hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT',
                    'CIFTI_STRUCTURE_CORTEX_RIGHT']
        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef = None
        ts_cortex = []
        for r in runs:
            # load the cifti
            ts_cifti = nb.load(fnames[r])

            # get the ts in surface for corticals
            ts_32k = util.surf_from_cifti(ts_cifti, hem_name)
            ts_parcel = []
            for hem in range(2):

                # get the average within parcels
                ts_parcel.append(
                    atlas_parcels[hem].agg_data(ts_32k[hem])
                )

            # concatenate them into a single array for correlation calculation
            ts_cortex.append(ts_parcel)
        return ts_cortex  # shape (n_tessl,P)

    def get_network_timecourse(self, networks, ts):
        """Regresses the group spatial map into the fMRI run.
        Returns the run-specific network timecourse.

        Args:
            networks (np.arry): 4D Network data of the signal components
                (default input networks are in MNI Space: 91 x 109 x 91 x nComponents )
            ts_vol (<nibabel CIFTI image object>): fMRI timeseries in volume
                Has to be in the same space as networks (91 x 109 x 91 x nTimepoints )
        Returns:
            ts_networks (np.ndarray):
                A numpy array (nTimepoints x nNetworks) with the fMRI timecourse for
                each resting-state network
        """
        X = networks.reshape(-1, networks.shape[3])
        Y = ts.reshape(-1, ts.shape[3])
        ts_networks = np.matmul(np.linalg.pinv(X), Y)

        return ts_networks  # shape (n_tessl, time_course)

    def get_cereb_connectivity(self, participant_id, cereb_atlas_map, runs=[0, 1, 2, 3],
                               type='IcoAll', cortical_atlas_parcels=None, networks=None):
        """Uses the original CIFTI files to produce cerebellar connectivity
           file

        """
        # Split type information on capital letters
        type_info = re.findall('[A-Z][^A-Z]*', type)

        hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT',
                    'CIFTI_STRUCTURE_CORTEX_RIGHT']

        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef = None
        for r, run in enumerate(runs):

            ts_cerebellum = am.get_data_cifti([fnames[run]], [cereb_atlas_map])
            ts_cerebellum = ts_cerebellum[0]

            # load the cifti
            ts_cifti = nb.load(fnames[run])

            if type_info[0] == 'Ico':
                ts_32k = util.surf_from_cifti(ts_cifti, hem_name)
                # get the ts in surface for corticals
                ts_parcel = []
                for hem in range(2):

                    # get the average within parcels
                    ts_parcel.append(
                        cortical_atlas_parcels[hem].agg_data(ts_32k[hem])
                    )

                # concatenate them into a single array for correlation calculation
                ts_seed = np.concatenate(ts_parcel, axis=1)
            elif type_info[0] == 'Net':
                # Regress network spatial map into the run's cortical data
                # (returns run-specific timecourse for each network)
                ts_vol = util.volume_from_cifti(ts_cifti)
                ts_seed = self.get_network_timecourse(
                    networks, ts_vol.get_fdata())
                ts_seed = ts_seed.T

            # Standardize the time series for easier calculation
            ts_cerebellum = util.zstandarize_ts(ts_cerebellum)
            ts_seed = util.zstandarize_ts(ts_seed)

            # Correlation calculation
            if coef is None:
                coef = np.empty((len(runs), ts_seed.shape[1],
                                 ts_cerebellum.shape[1]))
            N = ts_cerebellum.shape[0]
            coef[r, :, :] = ts_seed.T @ ts_cerebellum / N

        return coef

    def get_cortical_connectivity(self, participant_id, cortex_mask=None,
                                  runs=[0, 1, 2, 3], type='IcoAll',
                                  cortical_atlas_parcels=None,
                                  networks=None):
        """Uses the original CIFTI files to produce cortical connectivity
           file
        Args:
            participant_id (int): participant id
            cortex_mask (list): the list of L and R hemis cortex mask
            runs: index of runs
            type: the underlying cortical parcellation and type of extraction
                'IcoAll': Single estimate per session,
                          Correlation with cortical icosahedron parcels
                'IcoRun': Seperate estimates per run,
                          Correlation with cortical icosahedron parcels
                'NetAll': Single estimate per session,
                          Correlation with cortico-cerebellar resting-state networks
                'NetRun': Seperate estimates per run,
                          Correlation with cortico-cerebellar resting-state networks
                Defaults to 'IcoAll'.
            cortical_atlas_parcels: cortical random tessellation parcel object
            networks: group ICA networks
        Returns:
            [coef_1,coef_2]: List of cortical functional connectivity.
                             [left hemisphere, right hemisphere]
        """
        hem_name = ['CIFTI_STRUCTURE_CORTEX_LEFT',
                    'CIFTI_STRUCTURE_CORTEX_RIGHT']

        # get the file name for the cifti time series
        fnames = self.get_data_fnames(participant_id)
        coef_1, coef_2 = None, None
        for r, run in enumerate(runs):
            # load the cifti
            ts_cifti = nb.load(fnames[run])

            # get the ts in surface for corticals
            ts_32k = util.surf_from_cifti(
                ts_cifti, hem_name, mask_gii=cortex_mask)

            if type[0:3] == 'Ico':
                assert cortical_atlas_parcels is not None, \
                    "cortical_atlas_parcels must be given if extraction type is `Ico`!"
                ts_parcel = []
                for hem in range(2):
                    # get the average within parcels
                    ts_parcel.append(
                        cortical_atlas_parcels[hem].agg_data(ts_32k[hem]))

                # concatenate them into a single array for correlation calculation
                ts_parcel = np.concatenate(ts_parcel, axis=1)
            elif type[0:3] == 'Net':
                assert networks is not None, \
                    "networks must be given if extraction type is `Net`!"
                # Regress network spatial map into the run's wholebrain data
                # (returns run-specific timecourse for each network)
                ts_vol = util.volume_from_cifti(ts_cifti)
                ts = ts_vol.get_fdata()
                ts_parcel = self.get_network_timecourse(networks, ts)
                ts_parcel = ts_parcel.T

            # Standardize the time series for easier calculation
            ts_cortex = [util.zstandarize_ts(timeseries)
                         for timeseries in ts_32k]
            ts_parcel = util.zstandarize_ts(ts_parcel)

            # Correlation calculation
            if coef_1 is None:
                coef_1 = np.empty((len(runs), ts_parcel.shape[1],
                                   ts_cortex[0].shape[1]))  # (runs, parcels, vertices)
            if coef_2 is None:
                coef_2 = np.empty((len(runs), ts_parcel.shape[1],
                                   ts_cortex[1].shape[1]))  # (runs, parcels, vertices)

            N1 = ts_cortex[0].shape[0]
            N2 = ts_cortex[1].shape[0]
            coef_1[r, :, :] = ts_parcel.T @ ts_cortex[0] / N1
            coef_2[r, :, :] = ts_parcel.T @ ts_cortex[1] / N2

        return [coef_1, coef_2]  # shape (n_tessl,P)

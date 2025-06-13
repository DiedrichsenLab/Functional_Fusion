
class DataSetDmcc(DataSetMNIVol):
    def __init__(self, dir):
        super().__init__(dir)
        self.space = 'MNI152NLin2009cAsym'
        self.sessions = ['ses-axcpt-bas-mixed', 'ses-cuedts-bas-mixed', 'ses-stern-bas-mixed', 'ses-stroop-bas-mixed']
        self.default_type = 'CondHalf'
        self.cond_ind = 'cond_num'
        self.cond_name = 'cond_name'
        self.part_ind = 'knot_num'

    def get_data_fnames(self, participant_id, session_id=None, type='Cond'):
        """ Gets all raw data files
        Args:
            participant_id (str): Subject
            session_id (str): Session ID. Defaults to None.
            type (str): Type of data. Defaults to 'Cond' for task-based data. For rest data use 'Tseries'.
        Returns:
            fnames (list): List of fnames, last one is the resMS image
            T (pd.DataFrame): Info structure for regressors (reginfo)
        """
        dirw = self.estimates_dir.format(participant_id) + f'/{session_id}'
        # handle subjects with missing pro or rea sessions
        T = pd.read_csv(
            dirw + f'/{participant_id}_{session_id}_reginfo.tsv', sep='\t')
        
        if type == 'Contrast': # if you wanna look/work with at contrasts
            T = pd.read_csv(
                dirw + f'/{participant_id}_{session_id}_coninfo.tsv', sep='\t')
        
        
        if type[:4] == 'Cond' or type[:4] == 'Task' or type[:4] == 'Blck':
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.reg_id:02}_beta.nii' for i, t in T.iterrows()]
            # fnames.append(f'{dirw}/{participant_id}_{session_id}_resms.nii')
        elif type == 'Contrast':
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{t.run:02}_reg-{t.con_id:02}_con.nii' for i, t in T.iterrows()]
        elif type == 'Tseries' or type == 'FixTseries':
            fnames = [f'{dirw}/{participant_id}_{session_id}_run-{r:02}.nii' for r in T.run.unique().tolist()]
        return fnames, T

    def condense_data(self, data, info,
                      type='CondHalf',
                      participant_id=None,
                      ses_id=None):
        """ Extract data in a specific atlas space
        Args:
            participant_id (str): ID of participant
            atlas_maps (list): List of atlasmaps
            ses_id (str): Name of session
            type (str): Type of extraction:
                'CondHalf': Conditions with seperate estimates for first and second half of experient (Default)
                'CondRun': Conditions with seperate estimates per run
                    Defaults to 'CondHalf'.

        Returns:
            Y (list of np.ndarray):
                A list (len = numatlas) with N x P_i numpy array of prewhitened data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
            names: Names for CIFTI-file per row

        N.B.: Because some runs are missing for session 1-3, CondRun can only be run for session 04 (which has all runs for all subjects).
        Missing runs are: S3_sess03_MOTOR6, S3_sess01_MOTOR3, S3_sess01_MOTOR4, S3_sess01_MOTOR5, S4_sess01_MOTOR6, S4_sess02_MOTOR6 & S6_sess02_MOTOR2
        """
        # Depending on the type, make a new contrast
        info['half'] = (info.run % 2) + 1
        # n_cond = np.max(info.reg_id)
        
        if type == 'CondAll':
            data_info, C = agg_data(info, ['cond_num', 'cond_name'], ['knot_num', 'run'])
            data_info['names'] = [
                f'{d.cond_name}' for i, d in data_info.iterrows()]
        elif type == 'Contrast':
            data_info, C = agg_data(info, ['contrast_num', 'contrast_name'], ['knot_num', 'run'])
            data_info['names'] = [
                f'{d.contrast_name}' for i, d in data_info.iterrows()]
            
        

        # Prewhiten the data
        # data_n = prewhiten_data(data)
        # NOTE: I am currently using betas estimated using AFNI TentZero
        # It does not output ResMS and based on the documentation, it prewhitens the data
        # so the betas produced are already prewhitened.
        # data_n = prewhiten_data(data)
        data_n = data

        # Load the designmatrix and perform optimal contrast
        if type != 'Contrast':
            dir = self.estimates_dir.format(participant_id) + f'/{ses_id}'
            X = np.load(dir + f'/{participant_id}_{ses_id}_designmatrix.npy')
            reg_in = np.arange(C.shape[1], dtype=int)
            data_new = optimal_contrast(data_n, C, X, reg_in, baseline=None)
        else:
            data_new = data_n
            for i in range(len(data_n)):
                data_new[i] = pinv(C) @ data_n[i]

        
        return data_new, data_info

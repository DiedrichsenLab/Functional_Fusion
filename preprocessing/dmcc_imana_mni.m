function [ output_args ] = dmcc_imana_mni( what, varargin )

    % %========================================================================================================================
    %% ----- Initialize suit toolbox -----
    % Add dependencies to path
    if isdir('/Volumes/diedrichsen_data$/data')
        workdir='/Volumes/diedrichsen_data$/data';
    elseif isdir('/srv/diedrichsen/data')
        workdir='/srv/diedrichsen/data';
    else
        fprintf('Workdir not found. Mount or connect to server and try again.');
    end
    % check for SUIT installation
    if isempty(which('suit_isolate_seg')) % this function is only visible while SPM is actually "running" (not just on the path). This needs to happen for SUIT to run.
        warning('Cannot find SUIT, starting SPM12.'); % this should not happen since checked for in the beginning. Still leaving this snippet for robustnes (e.g. if someone closes SPM between starting Lead and pressing run).
        spm fmri
        if isempty(which('suit_isolate_seg')) % still not found.
        error('SUIT toolbox not found. Please install SUIT toolbox for SPM12 first (http://www.diedrichsenlab.org/imaging/suit.htm).');
        end
    end
    suit_defaults;
    
    
    %========================================================================================================================
    global base_dir
    global base_dir_mni
    global subj_str
    global ses_str
    
    base_dir = "/srv/diedrichsen/data/Cerebellum/DMCC"; % this is where tsv files are saved
    base_dir_mni = "/srv/diedrichsen/data/Cerebellum/DMCC/derivatives/fmriprep-1.3.2"; % parent directory of MNI space 
    
    %%% Freesurfer stuff
    path1 = getenv('PATH');
    path1 = [path1, sprintf('%s/freesurfer/6.0.0/bin', workdir)];
    setenv('PATH', path1);
    path1 = [path1, sprintf('%s/freesurfer/6.0.0/bin/surfreg', workdir)];
    setenv('PATH', path1);
    path1 = [path1, sprintf('%s/freesurfer/6.0.0/bin/xhemireg', workdir)];
    setenv('PATH', path1);
    path1 = [path1, ':/srv/software/freesurfer/6.0.0/fsfast/bin'];
    setenv('PATH', path1);
    path1 = [path1, ':/srv/software/freesurfer/6.0.0/mni/bin'];
    setenv('PATH', path1);
    setenv('FREESURFER_HOME','/srv/software/freesurfer/6.0.0');
    setenv(fullfile(base_dir_mni, 'surfaceFreesurfer'))
    setenv('SUBJECTS_DIR',fullfile(base_dir_mni, 'surfaceFreesurfer'));
    % setenv('PERL5LIB','/Applications/freesurfer/mni/Library/Perl/Updates/5.10.0');
    % setenv('PERL5LIB', '/Applications/freesurfer/mni/System/Library/Perl/5.8.6');
    
    % path1 = [path1 '/Applications/workbench/bin_macosx64'];
    % setenv('PATH', path1);
    
    % defining the names of other directories
    func_dir  = 'func';
    anat_dir  = 'anat';
    est_dir   = 'estimates';
    fs_dir    = 'surfaceFreeSurfer';
    wb_dir    = 'surfaceWB';
    mni_space = 'MNI152NLin2009cAsym'; 
    
    % names of the sessions: baseline: bas, reactivation: rea, probe: pro
    %%% TODO: I currently have only the baseline session, other sessions will be added later
    ses_str = {'bas', 'rea', 'pro'};
    
    % names of the tasks 
    task_str = {'Axcpt', 'Cuedts', 'Stern', 'Stroop'};
    % list of runs for each task
    %%% TODO: these are the runs for the baseline session
    run_list = {[1, 2], [1, 2], [1, 2], [1, 2]};
    run_enc_dir = {'mb4AP', 'mb4PA'}; % gradient encoding for each functional run
    
    % list of subjects
    subj_str = {'sub-f1027ao', 'sub-f1031ax', 'sub-f1342ku', 'sub-f1550bc',...
                'sub-f1552xo', 'sub-f1659oa', 'sub-f1670rz', 'sub-f1828ko',...
                'sub-f2157me', 'sub-f2499cq', 'sub-f2593wi', 'sub-f2648qw',...
                'sub-f2709ul', 'sub-f2968sp', 'sub-f3300jh', 'sub-f3387yq',...
                'sub-f3469wa', 'sub-f3526dz', 'sub-f3680fb', 'sub-f3681wf',...
                'sub-f3996sp', 'sub-f4138ge', 'sub-f4310gw', 'sub-f4354bs',...
                'sub-f4467ur', 'sub-f4796rs', 'sub-f4831tn', 'sub-f5001ob',...
                'sub-f5004cr', 'sub-f5094na', 'sub-f5094ya', 'sub-f5386yx',...
                'sub-f5416zj', 'sub-f5445nh', 'sub-f5635rv', 'sub-f5650zm',...
                'sub-f5930vp', 'sub-f6188io', 'sub-f6318if', 'sub-f6464bf',...
                'sub-f6950qp', 'sub-f7227ag', 'sub-f7688lh', 'sub-f7951pz',...
                'sub-f8113do', 'sub-f8194sp', 'sub-f8270up', 'sub-f8294bu',...
                'sub-f8298ds', 'sub-f8570ui', 'sub-f8710qa', 'sub-f8979ai',...
                'sub-f9057kp', 'sub-f9206gd', 'sub-f9271ex'};
    subj_id  = 1:length(subj_str);
    
    % number of TRs per task per run per session?
    % TODO: update with actual numbers
    % TODO: maybe fix this and hardcode number of TRs?
    % axcpt: 610, cuedts: 650, stern: 600, stroop: 540
    num_tr = [610, 650, 600, 540];
    %
    
    % number of dummy scans to be removed
    %%% TODO: update with actual numbers
    %%% do you even need to remove dummies? they might have already been removed
    num_dummy = 0; 
    % =========================================================================
    
    % =========================================================================
    % SUMMARY of the task design from the paper
    
    % =========================================================================
    
    % NOTES:
    % ANATOMICAL IMAGES:
    
    switch what
        case 'ANAT:unzip'  % unzipping gz files
            % simply unzips the anatomical files both native and MNI space
            % files
            % Example usage:dmcc_imana_mni('ANAT:unzip', 'sn', 21:26)
            sn = subj_id;
    
            vararginoptions(varargin, {'sn'});
            for s = sn
                fprintf('-%02d unziping %s anatomical\n', s, subj_str{s});
    
                % Get the directory of subjects anatomical
                subj_dir = fullfile(base_dir_mni, subj_str{s}, anat_dir);
                
                cd(subj_dir)
                gunzip('*.gz')
            end % sn (subjects)
        case 'ANAT:get_ants_def' % gets the deformation field nii image from ants
            % Example usage: dmcc_imana_mni('ANAT:get_ants_def', 'sn', 1)
            % NOTE: this cannot be done on the server due to this:
            % glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found 
            % on the server, use bash scripts to run the same command
            
            sn = subj_id;
    
            vararginoptions(varargin, {'sn'});
            for s = sn
                fprintf('-%02d %s creating deformation field image using ANTs\n', s, subj_str{s});
    
                
                % Get the directory of subjects anatomical
                subj_dir = fullfile(base_dir_mni, subj_str{s}, anat_dir);
                
                % get the deformation file name
                def_name = sprintf('%s_from-T1w_to-%s_mode-image_xfm.h5', subj_str{s}, mni_space);
                
                % get def file name
                def_filename = fullfile(subj_dir, def_name);
    
                % make the command string for moving from subject to mni
                command = sprintf('CompositeTransformUtil --disassemble %s y', def_filename);
                [status,cmdout] = system(command);
                fprintf(cmdout)
            end % sn (subjects)
            
            
            
        case 'FUNC:unzip'          % unzipping the files in mni space
            % TODO: figure out how to handle dummies. Do you want to remove
            % them here?
            % Example usage: dmcc_imana_mni('FUNC:unzip', 'sn', [2])
            sn = subj_id;
            
            vararginoptions(varargin, 'sn')
            for s = sn
                fprintf('-%02d unziping %s functionals\n', s, subj_str{s});
                % go to subject's directory
                func_subj_dir = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), func_dir);
                cd(func_subj_dir)
                % unzip the files
                gunzip('*.gz');
            end % sn (subjects)
        case 'FUNC:make_maskImage' % get the gray matter mask in functional space
            % Make maskImage in functional space
            % run this step after GLM!
            % Example usage: dmcc_imana_mni('FUNC:make_maskImage', 'sn', [1, 3])
            
            sn     = subj_id; % list of subjects
            vararginoptions(varargin, {'sn'});
            
            
            for s = sn
                % Get the directory of subjects anatomical and functional
                subj_anat_dir = fullfile(base_dir_mni, subj_str{s}, anat_dir);
                subj_func_dir = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), func_dir);
    
                fprintf('- make mask for %s\n', subj_str{s});
                cd(subj_func_dir);
                
                % loop over tasks
                for task_id = 1:length(task_str)
                    % get the task name
                    task_name = task_str{task_id};
                    
                    % loop over runs
                    for run = [1, 2]
                        nam{1}  = char(fullfile(subj_func_dir, sprintf('%s_ses-wave1bas_task-%s_acq-%s_run-%d_space-MNI152NLin2009cAsym_desc-aparcaseg_dseg.nii', subj_str{s}, task_str{task_id}, run_enc_dir{run}, run)));
                        nam{2}  = char(fullfile(subj_anat_dir, sprintf('%s_space-MNI152NLin2009cAsym_label-GM_probseg.nii', subj_str{s}))); % gray matter
                        % going with the voxels that have a gray matter probability
                        % of 0.3 and higher (not too conservative)
                        spm_imcalc(nam, sprintf('%s_ses-wave1bas_task-%s_acq-%s_run-%d_space-MNI152NLin2009cAsym_desc-gray_mask.nii', subj_str{s}, task_name, run_enc_dir{run}, run), 'i1 & (i2>0.3)')
    
                    end % runs
                end % task_id
            end % s (sn)            fprintf('- runs realigned for %s task %s\n',subj_str{s}, task_name);
        
        case 'GLM:design'        % glm design per task (uses canonical hrf) in mni space
            % this is kept as a separate case from native for my own sanity
            % the only difference between this case and the one above is the
            % images that will be used for GLM estimation
            % models each condition as a separate regressors
            % For conditions with multiple repetitions, one regressor
            % represents all the instances
            % dmcc_imana_mni('GLM:design', 'sn', [1])
            
            sn = subj_id;
            hrf_cutoff = Inf;
            glm = 1;
            vararginoptions(varargin, {'sn', 'hrf_cutoff', 'glm'});
            
            for s = sn
                func_subj_dir = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), func_dir);
                icond_uni = 1;
                Tu = []; % task/condition for all the tasks
                % loop over tasks
                for task_id = 1:length(task_str) 
                    % create a directory to save the design
                    subj_est_dir = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), est_dir, sprintf('glm%02d', glm), task_str{task_id});
                    dircheck(subj_est_dir)
                    
                    T = []; % task/condition + session + run info
                    J = []; % structure with SPM fields to make the design
                    
                    J.dir            = {char(subj_est_dir)};
                    J.timing.units   = 'secs';
                    J.timing.RT      = 1.2; % TR = 1.2 sec
                    J.timing.fmri_t  = 16;
                    J.timing.fmri_t0 = 1;
                    
                    % get the list of runs for the current task
                    runs = run_list{task_id};
                    
                    % loop through runs within the current sessions
                    for run = 1:length(runs)
                        
                        % fill in nifti image names for the current run
                        N = cell(num_tr(task_id) - num_dummy, 1); % preallocating!
                        for i = 1:(num_tr(task_id) - num_dummy)
    %                         sub-f1027ao_ses-wave1bas_task-Axcpt_acq-mb4AP_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii
                            N{i} = char(fullfile(func_subj_dir, sprintf('%s_ses-wave1bas_task-%s_acq-%s_run-%d_space-MNI152NLin2009cAsym_desc-preproc_bold.nii, %d', subj_str{s}, task_str{task_id}, run_enc_dir{run}, runs(run), i)));
                        end % i (image numbers)
                        J.sess(run).scans = N; % scans in the current runs
                        
                        % get the path to the tsv file
                        % NOTE make sure that the tsv files are downloaded into
                        % derivatives directory not the raw one
                        tsv_path = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), func_dir);
                        % get the tsvfile for the current run
    %                     A = dload(fullfile(tsv_path, 'sub-f1027ao_ses-wave1bas_task-Stroop_acq-mb4PA_run-2_events.tsv'))
                        D = dload(fullfile(tsv_path, sprintf('%s_ses-wave1bas_task-%s_acq-%s_run-%d_events.tsv', subj_str{s}, task_str{task_id}, run_enc_dir{run}, runs(run))));
                        
                        % first get the unique trial types
                        % I will use sorted so that the order is consistent
                        % across subjects
                        unique_trials = unique(D.trial_type, 'sorted');
                        % loop over trials within the current run and build up
                        % the design matrix
                        for ic = 1:length(unique_trials)
                            
                            % get the indices corresponding to the current
                            % condition.
                            idx = strcmp(D.trial_type, unique_trials{ic});
                            fprintf('* %d instances found for condition %s in run %02d\n', sum(idx), unique_trials{ic}, runs(run))
                            
                            
                            % filling in "reginfo"
                            TT.sn_id       = s;
                            TT.sn          = subj_str(s);
                            TT.task        = task_id;
                            TT.task_name   = task_str(task_id); 
                            TT.run         = runs(run);
                            TT.run_enc_dir = run_enc_dir(run);
                            TT.trial_type  = unique_trials(ic);
                            TT.cond_num    = ic;
                            TT.reg_id      = icond_uni;
                            TT.n_rep       = sum(idx);
                            
                            % filling in fields of J (SPM Job)
                            J.sess(run).cond(ic).name = unique_trials{ic};
                            J.sess(run).cond(ic).tmod = 0;
                            J.sess(run).cond(ic).orth = 0;
                            J.sess(run).cond(ic).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            
                            % get onset and duration (should be in seconds)
                            onset    = D.onset(idx) - (J.timing.RT*num_dummy);
                            if onset < 0
                                warning("negative onset found")
                            end
                            duration = D.duration(idx);
                            
                            J.sess(run).cond(ic).onset    = onset;
                            J.sess(run).cond(ic).duration = duration;
                            fprintf('**icond_uni is %d\n', icond_uni);
                            
                            icond_uni = icond_uni + 1;
                            
                            % add the condition info to the reginfo structure
                            T = addstruct(T, TT);
                            Tu = addstruct(Tu, TT);
                        end % ic (conditions)
                        
                        J.sess(run).multi     = {''};
                        J.sess(run).regress   = struct('name', {}, 'val', {});
                        J.sess(run).multi_reg = {''};
                        J.sess(run).hpf       = hrf_cutoff; % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
                    end % run (runs of current session)
                    
                    J.fact             = struct('name', {}, 'levels', {});
                    J.bases.hrf.derivs = [0 0];
                    J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
                    J.volt             = 1;
                    J.global           = 'None';
                    J.mask             = {char(fullfile(func_subj_dir,sprintf('%s_ses-wave1bas_task-%s_acq-%s_run-%d_space-MNI152NLin2009cAsym_desc-brain_mask.nii', subj_str{s}, task_str{task_id}, run_enc_dir{run}, run)))};
                    J.mthresh          = 0.05;
                    J.cvi_mask         = {char(fullfile(func_subj_dir,sprintf('%s_ses-wave1bas_task-%s_acq-%s_run-%d_space-MNI152NLin2009cAsym_desc-gray_mask.nii', subj_str{s}, task_str{task_id}, run_enc_dir{run}, run)))};
                    J.cvi              =  'fast';
    
                    spm_rwls_run_fmri_spec(J);
                    save(fullfile(J.dir{1},'SPM_info.mat'), '-struct', 'T');
                    dsave(fullfile(J.dir{1},sprintf('%s_ses-wave1bas_%s_reginfo.tsv', subj_str{s}, task_str{task_id})), T);
                    fprintf('- estimates for glm_%d task %s has been saved for %s \n', glm, task_str{task_id}, subj_str{s});            
                end % task_id
                Tu_path = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), est_dir, sprintf('glm%02d', glm));
                dsave(fullfile(Tu_path,sprintf('SPM_info.mat')), Tu);
                dsave(fullfile(Tu_path,sprintf('%s_ses-wave1bas_reginfo.tsv', subj_str{s})), Tu);
                fprintf(' FINISHED %s\n', subj_str{s})
            end % sn (subject) 
        case 'GLM:check_design'  % checking the design matrix
            % run GLM:make_design, GLM:estimate, and GLM:contrast before this step
            % Example usage:nishimoto_imana('GLM:check_design', 'sn', 1, 'ses', 1, 'glm', 1)
            
            sn       = subj_id; % list of subjects you want to inspect
            ses = 1;
            glm      = 1;           % glm number
            runs     = 1:12;         % list of runs you want to inspect
            
            vararginoptions(varargin, {'sn', 'ses', 'glm', 'runs'});
    
            for s = sn 
                
                % glm subject directory
                glm_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), ses_str{ses});
                % load SPM.mat file
                load(fullfile(glm_dir, 'SPM.mat'));
                
                for r = runs
                    % get the design matrix for run r
                    X = SPM.xX.X(SPM.Sess(r).row, SPM.Sess(r).col);
                    
                    % get the variance of beta estimates
                    indx = SPM.Sess(r).col;
                    Bcov = SPM.xX.Bcov(indx, indx);
                    
                    % create visualizations
                    h = figure('Name', sprintf('%s run %02d', subj_str{s}, r));
                    subplot(1, 2, 1)
                    imagesc(X); colorbar;
                    title(sprintf('Design Matrix for run %02d GLM %02d', r, glm));
                    subplot(1, 2, 2)
                    imagesc(Bcov); axis square; colorbar;
                    title(sprintf('Variance of Beta estimates %s for run %02d GLM %02d', r, glm));
                    
                    
                end % r (runs)
            end % s (sn)
        case 'GLM:estimate'      % estimate beta values
            % Example usage: dmcc_imana_mni('GLM:estimate', 'glm', 1, 'sn', [1])
            
            sn       = subj_id; % subject list
            glm      = 1;       % glm number
            
            vararginoptions(varargin, {'sn', 'glm'})
            
            for s = sn
                % loop over tasks
                for task_id = 1:length(task_str)
                    fprintf('- Doing glm estimation for glm %02d %s task %s\n', glm, subj_str{s}, task_str{task_id});
                    subj_est_dir = char(fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), est_dir, sprintf('glm%02d', glm), task_str{task_id}));
                    
                    load(char(fullfile(subj_est_dir,'SPM.mat')));
                    SPM.swd = subj_est_dir;
                    
                    spm_rwls_spm(SPM);
                end % task_id
            end % s (sn),
        case 'GLM:save_design'   % save design matrix
            % Example usage: dmcc_imana_mni('GLM:save_design', 'sn', 1)
            
            sn       = subj_id; % subject list
            glm      = 1;       % glm number
            
            vararginoptions(varargin, {'sn', 'glm'})
            for s = sn
                % loop over tasks
                for task_id = 1:length(task_str)
                    fprintf('- Getting design matrix for %s %s\n', subj_str{s}, task_str{task_id});
                    subj_est_dir = char(fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), est_dir, sprintf('glm%02d', glm), task_str{task_id}));
                    
                    load(char(fullfile(subj_est_dir,'SPM.mat')));
                    X = SPM.xX.xKXs.X;
                        save(fullfile(subj_est_dir, ...
                            sprintf('design_matrix.mat')), 'X', '-v7');
                        
    %                 T = load(fullfile(subj_est_dir, 'SPM_info.mat'));
                    
    %                 spm_rwls_spm(SPM);
                end % task_id
            end % s (sn),
        case 'GLM:T_contrast'    % make T contrasts for each condition
            %%% Calculating contrast images.
            % Example usage: dmcc_imana_mni('GLM:T_contrast', 'sn', 2:20, 'glm', 1)
            
            sn             = subj_id;    % subjects list
            glm            = 1;          % glm number
            baseline       = 'rest';     % contrast will be calculated against base (available options: 'rest')
            
            vararginoptions(varargin, {'sn', 'glm', 'baseline'})
            
            for s = sn
                
                % loop through tasks
                for task_id = 1:length(task_str)
                    % get the subject id folder name
                    fprintf('Contrasts for task %s %s\n', task_str{task_id}, subj_str{s})
                    glm_dir = fullfile(base_dir_mni, subj_str{s}, sprintf('ses-wave1bas'), est_dir, sprintf('glm%02d', glm), task_str{task_id});
                    
                    cd(glm_dir);
                    
                    % load the SPM.mat file
                    load(fullfile(glm_dir, 'SPM.mat'))
                    
                    SPM  = rmfield(SPM,'xCon');
                    T    = dload(fullfile(glm_dir,sprintf('%s_ses-wave1bas_%s_reginfo.tsv', subj_str{s}, task_str{task_id})));
                    
                    % t contrast for each condition type
                    utrial = unique(T.trial_type)';
                    idx = 1;
                    for ic = utrial
                        switch baseline
                            case 'myBase' % contrast vs future baseline :)))
                                % put your new contrasts here!
                            case 'rest' % contrast against rest
                                con                          = zeros(1,size(SPM.xX.X,2));
                                con(:,logical(strcmp(T.trial_type, ic)& (T.n_rep>0))) = 1;
                                con                          = con/abs(sum(con));
                        end % switch base
                        
                        % set the name of the contrast
                        contrast_name = sprintf('%s-%s', char(unique(T.trial_type(strcmp(T.trial_type, ic)))), baseline);
                        SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', con', SPM.xX.xKXs);
                        
                        idx = idx + 1;
                    end % ic (conditions)
                    
                    SPM = spm_contrasts(SPM,1:length(SPM.xCon));
                    save('SPM.mat', 'SPM','-v7.3');
                    SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
                    save(fullfile(glm_dir, 'SPM_light.mat'), 'SPM')
                    
                    % rename contrast images and spmT images
                    conName = {'con','spmT'};
                    for i = 1:length(SPM.xCon)
                        for n = 1:numel(conName)
                            oldName = fullfile(glm_dir, sprintf('%s_%2.4d.nii',conName{n},i));
                            newName = fullfile(glm_dir, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                            movefile(oldName, newName);
                        end % conditions (n, conName: con and spmT)
                    end % i (contrasts)
                end % task_id
            end % sn
        case 'GLM:F_contrast'    % make F contrast
            %%% Calculating contrast images.
            % 'SPM_light' is created in this step (xVi is removed as it slows
            % down code for FAST GLM).
            % Example1: nishimoto_imana('GLM:F_contrast', 'sn', 1, 'glm', 1, 'ses', 2)
            
            sn       = subj_id;   %% list of subjects
            ses      = 2;
            glm      = 1;             %% The glm number
            
            vararginoptions(varargin, {'sn', 'glm', 'ses'})
            
            for s = sn
                % subject name
                fprintf('- Calculating F contrast for %s %s \n', ses_str{ses}, subj_str{s});
                glm_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), ses_str{ses});
                load(fullfile(glm_dir, 'SPM.mat'))
                
                SPM  = rmfield(SPM,'xCon');
                cd(fullfile(glm_dir))
                T    = dload(fullfile(glm_dir, sprintf('%s_%s_reginfo.tsv', subj_str{s}, ses_str{ses})));
                
                % F contrast
                numConds = max(T.task); 
                con = zeros(numConds,size(SPM.xX.X,2));
                for i=1:numConds
                    con(i,T.task==i)=1-1/numConds;
                    con(i,T.task>0 & T.task~=i)=-1/numConds;
                end
                
                SPM.xCon(1) = spm_FcUtil('Set', 'Fcon', 'F', 'c',con',SPM.xX.xKXs);
                SPM = spm_contrasts(SPM,1:length(SPM.xCon));
                save('SPM.mat', 'SPM','-v7.3');
                SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
                save(fullfile(glmDir, subj_name, 'SPM_light.mat'), 'SPM');
    
            end % sn 
        case 'GLM:check'         % visually inspect design matrix
            % run GLM:make_design, GLM:estimate, and GLM:contrast before this step
            % Example usage:nishimoto_imana('GLM:check', 'sn', 1, 'ses', 1, 'glm', 1)
            
            sn       = subj_id; % list of subjects you want to inspect
            ses      = 1;       % session number
            glm      = 1;       % glm number
            runs     = 1:5;     % list of runs you want to inspect or the number of the run you want to expect
            
            vararginoptions(varargin, {'sn', 'ses', 'glm', 'runs'});
    
            for s = sn 
                % glm subject directory
                glm_dir = fullfile(base_dir,subj_str{s}, est_dir, sprintf('glm%02d', glm), ses_str{ses});
                
                % load SPM.mat file
                load(fullfile(glm_dir, 'SPM.mat'));
                
                for r = runs
                    % get the design matrix for run r
                    X = SPM.xX.X(SPM.Sess(r).row, SPM.Sess(r).col);
                    
                    % get the variance of beta estimates
                    indx = SPM.Sess(r).col;
                    Bcov = SPM.xX.Bcov(indx, indx);
                    
                    % create visualizations
                    h = figure('Name', sprintf('%s run %02d', subj_str{s}, r));
                    subplot(1, 2, 1)
                    imagesc(X); colorbar;
                    title(sprintf('Design Matrix session %02d for run %02d GLM %02d', ses, r, glm));
                    subplot(1, 2, 2)
                    imagesc(Bcov); axis square; colorbar;
                    title(sprintf('Variance of Beta estimates session %02d for run %02d GLM %02d', ses, r, glm));
                    
                    keyboard;
                    
                end % r (runs)
            end % s (sn)     
                
        case 'SURF:reconall'       % Freesurfer reconall routine
            % Calls recon-all, which performs, all of the
            % FreeSurfer cortical reconstruction process
            % Example usage: dmcc_imana_mni('SURF:reconall', 'sn', 1)
            % NOTE: if you get textScalar error, open freesurfer_reconall and
            % change the system comman line to:
            % system(sprintf('recon-all -s %s  -i  %s  -all -cw256', subj_name, anatomical));
            
            sn   = subj_id; % subject list
            
            vararginoptions(varargin, {'sn'});
            % set freesurfer directory
            subj_fs_dir = fullfile(base_dir_mni, fs_dir);
            
            for s = sn
                fprintf('- recon-all %s\n', subj_str{s});
                subj_dir = fullfile(base_dir_mni, subj_str{s}, anat_dir);
                freesurfer_reconall(subj_fs_dir, subj_str{s}, ...
                                    fullfile(subj_dir,sprintf('%s_desc-preproc_T1w.nii', subj_str{s})));
            end % s (sn)        
        case 'SURF:xhemireg'       % Cross-register surfaces left / right hem
            % surface-based interhemispheric registration
            % example: dmcc_imana_mni('SURF:xhemireg', 'sn', [1])
            
            sn   = subj_id; % list of subjects
    
            vararginoptions(varargin, {'sn'})
            
            % set freesurfer directory
            fs_dir = fullfile(base_dir_mni, 'surfaceFreeSurfer');
            
            for s = sn
                fprintf('- xhemiregl %s\n', subj_str{s});
                freesurfer_registerXhem(subj_str(s), fs_dir,'hemisphere', [1 2]); % For debug... [1 2] orig
            end % s (sn)
        case 'SURF:fs2wb'          % Resampling subject from freesurfer fsaverage to fs_LR
            % Example usage: dmcc_imana_mni('SURF:fs2wb', 'sn', [1])
            
            sn   = subj_id; % list of subjects
            hemi = [1, 2];      % list of hemispheres
            
            vararginoptions(varargin, {'sn', 'hemi'});
            
            % set freesurfer directory
            fs_dir = fullfile(base_dir_mni, 'surfaceFreeSurfer');
            wb_dir  = char(fullfile(base_dir_mni, wb_dir, 'data'));
            for s = sn 
                fprintf('- fs2wb %s\n', subj_str{s});
                surf_resliceFS2WB(subj_str{s},...
                                  fs_dir,...
                                  wb_dir, ...
                                  'hemisphere',[1 2], ...
                                  'resolution', '32k')            
            end % s (sn)
        case 'SURF:run_all'        % Pipeline running all of surface preprocessing routines
            % Example usage: dmcc_imana_mni('SURF:run_all', 'sn', 1:20)
            
            sn = subj_id;
            
            vararginoptions(varargin, {'sn'});
            dmcc_imana_mni('SURF:reconall', 'sn', sn)
            dmcc_imana_mni('SURF:fs2wb', 'sn', sn);
                        
        case 'SUIT:isolate_segment'    % Segment cerebellum into grey and white matter
            % Example usage: dmcc_imana_mni('SUIT:isolate_segment', 'sn', 1);
            
            sn = subj_id;
            % results of both will be saved under suit directory
            
            vararginoptions(varargin, {'sn'});
            
            for s = sn
                fprintf('- Isolate and segment the cerebellum for %s\n', subj_str{s})
                spm_jobman('initcfg')
                
                % Get the directory of subjects anatomical
                subj_dir = fullfile(base_dir_mni, subj_str{s}, anat_dir);
                % Get the name of the anatpmical image
                % sub-f1027ao_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii
                anat_name = sprintf('%s_space-%s_desc-preproc_T1w.nii', subj_str{s}, mni_space);
                
                
                suit_subj_dir = fullfile(base_dir_mni, subj_str{s}, 'suit');
                dest   = char(fullfile(suit_subj_dir, sprintf('%s_space-%s_desc-preproc_T1w.nii', subj_str{s}, mni_space)));
    
                dircheck(suit_subj_dir);
                
                source = char(fullfile(subj_dir, anat_name));
                
                copyfile(source,dest);
                
                % go to subject directory for suit and isolate segment
                cd(fullfile(suit_subj_dir));
                suit_isolate_seg(fullfile(suit_subj_dir, sprintf('%s_space-%s_desc-preproc_T1w.nii', subj_str{s}, mni_space)), 'keeptempfiles', 1);
            end % s (sn)
        case 'SUIT:normalise_dartel'   % SUIT normalization using dartel
            % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
            % example usage: dmcc_imana_mni('SUIT:normalise_dartel')
            sn = subj_id; %subjNum
            vararginoptions(varargin, 'sn');
            
            for s = sn
                suit_subj_dir = fullfile(base_dir_mni, subj_str{s}, 'suit');
                cd(suit_subj_dir)
                job.subjND.gray       = {char(fullfile(suit_subj_dir, sprintf('c_%s_space-%s_desc-preproc_T1w_seg1.nii', subj_str{s}, mni_space)))};
                job.subjND.white      = {char(fullfile(suit_subj_dir, sprintf('c_%s_space-%s_desc-preproc_T1w_seg2.nii', subj_str{s}, mni_space)))};
                job.subjND.isolation  = {char(fullfile(suit_subj_dir, sprintf('c_%s_space-%s_desc-preproc_T1w_pcereb.nii', subj_str{s}, mni_space)))};
    
                suit_normalize_dartel(job);
            end % s (subjects)    
        case 'SUIT:save_dartel_def'    % save the deformation field image
            % Saves the dartel flow field as a deformation file. 
            % example usage: dmcc_imana_mni('SUIT:save_dartel_def', 'sn', [1, 2, 3])
            sn = subj_id; %subjNum
            vararginoptions(varargin, 'sn');
    
            for s = sn
                fprintf('- saving suit deformation for %s\n', subj_str{s})
                suit_subj_dir = fullfile(base_dir_mni, subj_str{s}, 'suit');
    
                cd(suit_subj_dir);
                anat_name = sprintf('%s_space-MNI152NLin2009cAsym_desc-preproc_T1w', subj_str{s});
                suit_save_darteldef(anat_name);
            end % s (subjects)
        case 'SUIT:make_func_mask'     % make a cerebellar mask for functional images
            % Example usage: dmcc_imana_mni('SUIT:make_func_mask', 'glm', 1, 'sn', 1)
            
            sn       = subj_id; % list of subjects
            glm      = 1;           % glm number
            
            vararginoptions(varargin, {'sn', 'glm'})
    
            
            for s = sn
                suit_dir = fullfile(base_dir_mni, subj_str{s}, 'suit');
                for task_id = 1 % no need to loop through tasks and have separate masks per tasks. They are the same!
                    % get the task name
                    task_name = task_str{task_id};
                    glm_dir = fullfile(base_dir_mni, subj_str{s}, 'ses-wave1bas', 'estimates', sprintf('glm%02d', glm), task_name);
                    
                    mask  = char(fullfile(glm_dir, 'mask.nii')); % mask for functional image
                    %
                    pcereb  = char(fullfile(suit_dir, sprintf('c_%s_space-%s_desc-preproc_T1w_pcereb.nii', subj_str{s}, mni_space)));
                    cereb_gray = char(fullfile(suit_dir, sprintf('c_%s_space-%s_desc-preproc_T1w_seg1.nii', subj_str{s}, mni_space))); %
                    suit_glm_dir = char(fullfile(base_dir_mni, subj_str{s}, 'suit')); dircheck(suit_glm_dir);
                    omask = char(fullfile(suit_glm_dir, 'maskbrainSUITGrey.nii'));
    %                 omask = char(fullfile(suit_glm_dir, sprintf('%s_%s_desc_preproc_suit_gray.nii', subj_str{s}, mni_space))); % output mask image - grey matter
                    
                    cd(suit_dir);
                    spm_imcalc({mask,pcereb,cereb_gray}, omask, 'i1>0 & i2>0 & i3>0.01', {});
    
                end % task_id
                
            end % s (sn)           
    
    end
    
    end
    
    
    % Checking for directories and creating if not exist
    function dircheck(dir)
    if ~exist(dir,'dir')
        warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
        mkdir(dir);
    end
    end
    
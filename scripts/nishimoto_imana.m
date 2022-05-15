function [ output_args ] = nishimoto_imana( what, varargin )

% %========================================================================================================================
% PATH DEFINITIONS

% Add dependencies to path
% if isdir('/Volumes/diedrichsen_data$/data')
%     workdir='/Volumes/diedrichsen_data$/data';
% elseif isdir('/srv/diedrichsen/data')''
%     workdir='/srv/diedrichsen/data';
% else
%     fprintf('Workdir not found. Mount or connect to server and try again.');
% end
% addpath(sprintf('%s/../matlab/spm12',workdir));
% addpath(sprintf('%s/../matlab/spm12/toolbox/suit/',workdir));
% addpath(sprintf('%s/../matlab/dataframe',workdir));
% addpath(sprintf('%s/../matlab/imaging/tools/',workdir));

%% ----- Initialize suit toolbox -----
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

% base_dir = sprintf('%s/FunctionalFusion/Nishimoto_103Task/',workdir);
base_dir = '/Users/ladan/Documents/DATA/nishimoto';

%%% Freesurfer stuff
path1 = getenv('PATH');
path1 = [path1, ':/Applications/freesurfer/bin'];
setenv('PATH', path1);
path1 = [path1, ':/Applications/freesurfer/fsfast/bin'];
setenv('PATH', path1);
path1 = [path1, ':/Applications/freesurfer/mni/bin'];
setenv('PATH', path1);
setenv('FREESURFER_HOME','/Applications/freeSurfer');
setenv(fullfile(base_dir, 'surfaceFreesurfer'))
setenv('SUBJECTS_DIR',fullfile(base_dir, 'surfaceFreesurfer'));
% setenv('PERL5LIB','/Applications/freesurfer/mni/Library/Perl/Updates/5.10.0');
% setenv('PERL5LIB', '/Applications/freesurfer/mni/System/Library/Perl/5.8.6');

path1 = [path1 '/Applications/workbench/bin_macosx64'];
setenv('PATH', path1);

% defining the names of other directories
func_dir = 'func';
anat_dir = 'anat';
est_dir  = 'est';
fs_dir   = 'surfaceFreeSurfer';
wb_dir   = 'surfaceWB';

% list of subjects
subj_str = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06'};
subj_id  = [1, 2, 3, 4, 5, 6];

% AC coordinates
loc_AC = {[-103, -140, -140],...       %sub-01
          [-103, -133, -142],...       %sub-02
          [-101, -138, -130],...       %sub-03
          [-104, -147, -139],...       %sub-04
          [-100, -134, -134],...       %sub-05
          [-102, -134, -145],...       %sub-06
        };
numDumm = 3; % based on the paper and comments on https://openneuro.org/datasets/ds002306/versions/1.1.0

sess = [1, 2, 3];
run_list  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18];
switch what
    case 'ANAT:reslice_lpi' % reslice anatomical to LPI
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});

        for s = sn
            fprintf('- Reslicing %s anatomical to LPI\n', subj_str{s});
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            
            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w', subj_str{s});
            
            % Reslice anatomical image to set it within LPI co-ordinate frames
            source  = fullfile(subj_dir, sprintf('%s.nii', anat_name));
            dest    = fullfile(subj_dir, sprintf('%s_lpi.nii', anat_name));
            if ~isfile(source) && isfile(sprintf('%s.gz', source))  % unzip file
                gunzip(sprintf('%s.gz', source));
            end
            spmj_reslice_LPI(source,'name', dest);
            
            % In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
        end % sn (subjects)
    case 'ANAT:center_ac'   % recenter to AC (manually retrieve coordinates)
        % run spm display to get the AC coordinates
        fprintf('MANUALLY RETRIEVE AC COORDINATES')
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Centre AC for %s\n', subj_str{s});
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            
            % Get the name of the anatomical image
            anat_name = sprintf('%s_T1w_lpi.nii', subj_str{s});
            
            img             = fullfile(subj_dir, anat_name);
            V               = spm_vol(img);
            dat             = spm_read_vols(V);
            %%oldOrig         = V.mat(1:3,4);
            %%V.mat(1:3,4)    = oldOrig-loc_AC{sn(s)};
            V.mat(1:3,4)    = loc_AC{s};
            spm_write_vol(V,dat);
        end % s (subjects)
    case 'ANAT:segment'     % segment the anatomical image
        % also saves the bias corrected anatomical
        % check results when done
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        SPMhome = fileparts(which('spm.m'));
        J       = []; % spm jobman
        for s = sn
            fprintf('- Anatomical segmentation for %s\n', subj_str{s});
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            
            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w_lpi.nii', subj_str{s}); 
            J.channel.vols     = {fullfile(subj_dir,sprintf('%s,1', anat_name))};
            J.channel.biasreg  = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write    = [1 0];
            J.tissue(1).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus  = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus  = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus  = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus  = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus  = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus  = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];
            
            J.warp.mrf     = 1;
            J.warp.cleanup = 1;
            J.warp.reg     = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg  = 'mni';
            J.warp.fwhm    = 0;
            J.warp.samp    = 3;
            J.warp.write   = [1 1];
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
        end % s (subject)
    
    case 'FUNC:rm_dumm' % removing dummies and renaming
        % removes dummies from the beginning of the functional images and
        % save the new images with new names so that we don't lose the
        % original images. For code efficiency, it will also rename the tsv
        % files
        % Example usage: nishimoto_imana('FUNC:rm_dumm', 'sn', [1])
        sn = subj_id;
        
        vararginoptions(varargin, 'sn')
        
        for s = sn
            % go to subject's directory
            func_subj_dir = fullfile(base_dir, subj_str{s}, func_dir);
            cd(func_subj_dir)
%             funScans = dir('*bold*.nii');
            funScans = temp_list;
            for i = 1:length(funScans)
                srcfilename = sprintf('%s%sbold.nii', subj_str{s}, funScans{i});
                desfilename = sprintf('%s_run-%02d.nii', subj_str{s}, i);
                mkdir temp;
                unzip(srcfilename)
                
                spm_file_split(funScans{i}.name,'temp');
                spm_file_split(sprintf('%s%sbold.nii', subj_str{s}, funScans{i}),'temp');
                
                cd temp;
                list = dir('sub-*.nii');
                list = list(numDummys+1:end);  % Remove dummies
                V    = {list(:).name};
                spm_file_merge(V,srcfilename);
                movefile(srcfilename,fullfile(func_subj_dir, desfilename))
                cd ..
                rmdir('temp','s');
                % change the names of the tsv files
                tsv_source = sprintf('%s%sevents.tsv', subj_str{s}, funScans{i});
                tsv_dest = sprintf('%s_run-%02d_events.tsv', subj_str{s}, i);
                movefile(fullfile(func_subj_dir, tsv_source), fullfile(func_subj_dir, tsv_dest))
                fprintf('-Dummies removed for run-%02d done for %s \n', i, subj_str{s});
            end % i (number of runs: funScans)
        end % sn (subjects)    
    case 'FUNC:realign'
    case 'FUNC:move_data' 
    case 'FUNC:coreg' 
    case 'FUNC:make_samealign'
    case 'FUNC:make_maskImage'
        
    case 'GLM:design1' % make the design matrix for the glm
        % models each condition as a separate regressors
        % For conditions with multiple repetitions, one regressor
        % represents all the instances
        % nishimoto_imana('GLM:design1', 'sn', [1])
        
        sn = subj_id;
        hrf_cutoff = Inf;
        vararginoptions(varargin, {'sn', 'hrf_cutoff'});
        
        prefix = 'r'; % prefix of the preprocessed epi we want to use
        glm = 1;
        for s = sn
            func_subj_dir = fullfile(base_dir, subj_str{s}, func_dir);
            % loop over runs
            for ss = [1, 2]
                % create a directory to save the design
                subj_est_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('ses-%d', ss));
                dircheck(subj_est_dir)
                
                T = []; % task/condition + session + run info
                J = []; % structure with SPM fields to make the design
                
                J.dir            = {subj_est_dir};
                J.timing.units   = 'secs';
                J.timing.RT      = 2.0;
                J.timing.fmri_t  = 16;
                J.timing.fmri_t0 = 1;
                
                % get the list of runs for the current session
                runs_list = sess{ss};
                
                % loop through runs within the current sessions
                for run = 1:length(runs_list)
                    
%                   % fill in nifti image names for the current run
                    N = cell(numTRs - numDummys, 1); % preallocating!
                    for i = 1:(numTRs-numDummys)
                        N{i} = fullfile(func_subj_dir, sprintf('%s%s_run-%02d.nii', prefix, subj_str{s}, i));
                    end % i (image numbers)
                    J.sess(run).scans = N; % scans in the current runs
                    
                    % get the path to the tsv file
                    tsv_path = fullfile(base_dir, subj_str{s}, func_dir);
                    % get the tsvfile for the current run
                    D = dload(fullfile(tsv_path, sprintf('%s_run-%02d_events.tsv', subj_str{s}, run)));
                    
                    % get unique conditions
                    unique_conds = unique(D.trial_type, 'stable');
                    
                    % loop over trials within the current run and build up
                    % the design matrix
                    for ic = 1:length(unique_conds)
                        % get the indices corresponding to the current
                        % condition.
                        % this line is necessary as there are some
                        % conditions with more than 1 repetition
                        idx = strcmp(D.trial_type, unique_conds{ic});
                        fprintf('* %d instances found for condition %s in run %02d\n', sum(idx), unique_conds{ic}, run)
                        
                        % filling in "reginfo"
                        TT.sn    = s;
                        TT.sess  = ss;
                        TT.run   = run;
                        TT.CN    = unique_conds(ic);
                        TT.cond  = ic;
                        TT.n_rep = sum(idx);
                        
                        % filling in fields of J (SPM Job)
                        J.sess(run).cond(ic).name = unique_conds{ic};
                        J.sess(run).cond(ic).tmod = 0;
                        J.sess(run).cond(ic).orth = 0;
                        J.sess(run).cond(ic).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        
                        % get onset and duration (should be in seconds)
                        onset    = D.onset(idx) - (J.timing.RT*numDummys);
                        duration = D.duration(idx); 
                        
                        J.sess(run).cond(ic).onset    = onset;
                        J.sess(run).cond(ic).duration = duration;
                        
                        % add the condition info to the reginfo structure
                        T = addstruct(T, TT);
                    end % ic (conditions)
                    
                    J.sess(run).multi     = {''};
                    J.sess(run).regress   = struct('name', {}, 'val', {});
                    J.sess(run).multi_reg = {''};
                    J.sess(run).hpf       = hrf_cutoff; % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
                end % run (runs of current session)
                
%                 J.fact             = struct('name', {}, 'levels', {});
%                 J.bases.hrf.derivs = [0 0];
%                 J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
%                 J.volt             = 1;
%                 J.global           = 'None';
%                 J.mask             = {fullfile(imagingDir, task_names{task_num}, subj_name,'rmask_noskull.nii,1')};
%                 J.mthresh          = 0.05;
%                 J.cvi_mask         = {fullfile(imagingDir, task_names{task_num}, subj_name,'rmask_gray.nii')};
%                 J.cvi              =  'fast';
%                 
%                 spm_rwls_run_fmri_spec(J);
                
                save(fullfile(J.dir{1},sprintf('%s_ses-%d_reginfo.tsv', subj_str{s}, ss)), '-struct', 'T');
                fprintf('- estimates for glm_%d session %d has been saved for %s \n', glm, ss, subj_str{s});
            end % ss (session)
            
            
        end % sn (subject)
        
    case 'GLM:estimate'
    case 'GLM:contrast' 
        
        
    case 'SURF:reconall'       % Freesurfer reconall routine
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        % Example usage: nishimoto_imana('SURF:reconall', 'sn', 1)
        
        sn   = subj_id; % subject list
        
        vararginoptions(varargin, {'sn'});
        % set freesurfer directory
        subj_fs_dir = fullfile(base_dir, fs_dir);
        
        for s = sn
            fprintf('- recon-all %s\n', subj_str{s});
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            freesurfer_reconall(subj_fs_dir, subj_str{s}, ...
                                fullfile(subj_dir,sprintf('%s_T1w_lpi.nii', subj_str{s})));
        end % s (sn)
    case 'SURF:xhemireg'       % Cross-register surfaces left / right hem
        % surface-based interhemispheric registration
        % example: nishimoto_imana('SURF:xhemireg', 'sn', 95)
        
        sn   = subj_id; % list of subjects

        vararginoptions(varargin, {'sn'})
        
        % set freesurfer directory
        fs_dir = fullfile(base_dir, 'FreeSurfer');
        
        for s = sn
            fprintf('- xhemiregl %s\n', subj_str{s});
            freesurfer_registerXhem(subj_str(s), fs_dir,'hemisphere', [1 2]); % For debug... [1 2] orig
        end % s (sn)
    case 'SURF:map_ico'        % Align to the new atlas surface (map icosahedron)
        % Resampels a registered subject surface to a regular isocahedron
        % This allows things to happen in atlas space - each vertex number
        % corresponds exactly to an anatomical location
        % Makes a new folder, called ['x' subj] that contains the remapped subject
        % Uses function mri_surf2surf
        % mri_surf2surf: resamples one cortical surface onto another
        % Example usage: nishimoto_imana('SURF:map_ico', 'sn', 95)
        
        sn = subj_id; % list of subjects
        
        vararginoptions(varargin, {'sn'});
        
        % set freesurfer directory
        fs_dir = fullfile(base_dir, 'FreeSurfer');
        for s = sn
            fprintf('- map_ico %s\n', subj_str{s});
            freesurfer_mapicosahedron_xhem(subj_str{s}, fs_dir ,'smoothing',1,'hemisphere',[1, 2]);
        end % s (sn)
    case 'SURF:fs2wb'          % Resampling subject from freesurfer fsaverage to fs_LR
        % Example usage: nishimoto_imana('SURF:fs2wb', 'sn', [1], 'res', 32)
        
        sn   = subj_id; % list of subjects
        res  = 32;          % resolution of the atlas. options are: 32, 164
        hemi = [1, 2];      % list of hemispheres
        
        vararginoptions(varargin, {'sn', 'res', 'hemi'});
        
        % set freesurfer directory
        fs_dir = fullfile(base_dir, 'FreeSurfer');
        
        for s = sn 
            fprintf('- fs2wb %s\n', subj_str{s});
            wb_subj_dir  = fullfile(base_dir, subj_str{s}, wb_dir);
            dircheck(wb_subj_dir)
            surf_resliceFS2WB(subj_str{s}, fs_dir, wb_subj_dir, 'hemisphere', hemi, 'resolution', sprintf('%dk', res))
        end % s (sn)
    case 'SURF:run_all'
        % Example usage: nishimoto_imana('SURF:run_all')
        nishimoto_imana('SURF:reconall')
        nishimoto_imana('SURF:xhemireg')
        nishimoto_imana('SURF:map_ico')
        nishimoto_imana('SURF:fs2wb')
        
        
    case 'SUIT:isolate_segment' % Segment cerebellum into grey and white matter
        % Example usage: nishimoto_bids_imana('SUIT:isolate_segment', 'sn', 1);
        
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Isolate and segment the cerebellum for %s\n', subj_str{s})
            spm_jobman('initcfg')
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w_lpi.nii', subj_str{s});

            suit_subj_dir = fullfile(base_dir, subj_str{s}, 'suit');
            dircheck(suit_subj_dir);
            
            source = fullfile(subj_dir, anat_name);
            dest   = fullfile(suit_subj_dir, sprintf('%s_T1w.nii',subj_str{s}));
            
            copyfile(source,dest);
            
            % go to subject directory for suit and isolate segment
            cd(fullfile(suit_subj_dir));
            suit_isolate_seg({fullfile(suit_subj_dir, sprintf('%s_T1w.nii', subj_str{s}))}, 'keeptempfiles', 1);
%             suit_isolate_seg({fullfile(suit_subj_dir, sprintf('%s_T1w.nii', subj_str{s}))}, 'keeptempfiles', 1);
        end % s (sn)
    case 'SUIT:correct_cereb_mask' % Corrected cerebellum cortex mask
        % uses the results from SUIT:isolate_segment to create corrected
        % masks for the cerebellum and cortex. It removes buffer voxels
        % from both masks and create cereb_prob_corr_grey.nii and
        % cortical_mask_grey_corr.nii
        % Example usage: nishimoto_bids_imana('SUIT:correct_cereb_mask', 'sn', 1)
        
        sn   = returnSubjs; % subject list        
        vararginoptions(varargin, {'sn'});
        
        for s = sn 
            
            % Get the directory of subjects suit (where the result of segmentation is stored)

            suit_subj_dir = fullfile(base_dir, subj_str{s}, 'suit');

            % cortexGrey : cortex grey matter mask
            cortexGrey = fullfile(suit_subj_dir, sprintf('c7%s_T1w.nii', subj_str{s}));
            % cerebGrey  : cerebellar grey matter mask
            cerebGrey  = fullfile(suit_subj_dir, sprintf('c1%s_T1w.nii', subj_str{s}));
            % bufferVox  : buffer voxels
            bufferVox = fullfile(suitDir, 'anatomicals', subj_name, sprintf('sess-%02d', ss), 'buffer_voxels.nii');
            
            % isolate overlapping voxels
            spm_imcalc({cortexGrey,cerebGrey}, bufferVox,'(i1.*i2)')
            
            % mask buffer
            spm_imcalc({bufferVox},bufferVox,'i1>0')
            
            cerebGreyC  = fullfile(suitDir, 'anatomicals', subj_name, sprintf('sess-%02d', ss), 'cereb_prob_corr_grey.nii');
            cortexGreyC = fullfile(suitDir, 'anatomicals', subj_name, sprintf('sess-%02d', ss), 'cortical_mask_grey_corr.nii');
            
            % remove buffer from the cerebellum
            spm_imcalc({cerebGrey,bufferVox}, cerebGreyC,'i1-i2');
            
            % remove buffer from cortex
            spm_imcalc({cortexGrey,bufferVox}, cortexGreyC,'i1-i2');
            
        end % s (sn)
    case 'SUIT:normalize_darte'   
    case 'SUIT:reslice'                      %reslice cerebellum into suit space to check normalization
        % 'suit_normalise_dartel'.
        % example: nishimoto_bids_imana('SUIT:reslice','anatomical','pcereb')
        % make sure that you reslice into 2mm^3 resolution
        vararginoptions(varargin, {'sn'});
        type=varargin{1}; % 'betas' or 'contrast' or 'ResMS' or 'cerebellarGrey'
        mask=varargin{2}; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        subjs=length(subj_str);
        for s=1:subjs
            switch type
                case 'anatomical'
                    subj_dir = fullfile(base_dir, subj_str{s}, 'suit');
                    % Get the name of the anatpmical image
                    anat_name = sprintf('%s_T1w.nii', subj_str{s});
                    J.channel.vols     = {fullfile(subj_dir,sprintf('%s,1', anat_name))};
                    source=dir(fullfile(subj_dir,sprintf('%s', anat_name))); % images to be resliced
                    cd(subj_dir);
            end
            job.subj.affineTr = {fullfile(subj_dir,sprintf('%c_s_T1w_seg8.mat', subj_str{s}))};
            job.subj.flowfield= {fullfile(subj_dir,sprintf('y_%s_T1w.nii', subj_str{s}))};
            job.subj.resample = {source.name};
            job.subj.mask     = {fullfile(subj_dir,sprintf('c_%s_T1w_pcereb.nii', subj_str{s}))};
            job.vox           = [2 2 2];
            suit_reslice_dartel(job);
            fprintf('- Resliced %s anatomical to SUIT\n', subj_str{s});
        end
end

end


% Checking for directories and creating if not exist
function dircheck(dir)
if ~exist(dir,'dir')
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end
end

function unzip()
if ~isfile(source) && isfile(sprintf('%s.gz', source))  % unzip file
    gunzip(sprintf('%s.gz', source));
end
end




function [ output_args ] = nishimoto_imana( what, varargin )

% %========================================================================================================================
% PATH DEFINITIONS

% Add dependencies to path
if isdir('/Volumes/diedrichsen_data$/data')
    workdir='/Volumes/diedrichsen_data$/data';
elseif isdir('/srv/diedrichsen/data')
    workdir='/srv/diedrichsen/data';
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end
addpath(sprintf('%s/../matlab/spm12',workdir));
addpath(sprintf('%s/../matlab/spm12/toolbox/suit/',workdir));
addpath(sprintf('%s/../matlab/dataframe',workdir));
addpath(sprintf('%s/../matlab/imaging/tools/',workdir));
addpath(sprintf('%s/../matlab/imaging/coregtool/',workdir));

% Some suit functions are not released yet, hence clone the SUIT develop
% branch. Run:
% git clone https://github.com/jdiedrichsen/suit.git;
% git checkout Develop
addpath('~/Matlab/suit')


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

base_dir = sprintf('%s/FunctionalFusion/Nishimoto_103Task/raw/',workdir);
% base_dir = '/Users/ladan/Documents/DATA/nishimoto';

%%% Freesurfer stuff
path1 = getenv('PATH');
path1 = [path1, ':/srv/software/freesurfer/6.0.0/bin'];
setenv('PATH', path1);
path1 = [path1, ':/srv/software/freesurfer/6.0.0/fsfast/bin'];
setenv('PATH', path1);
path1 = [path1, ':/srv/software/freesurfer/6.0.0/mni/bin'];
setenv('PATH', path1);
setenv('FREESURFER_HOME','/srv/software/freesurfer/6.0.0');
setenv(fullfile(base_dir, 'surfaceFreesurfer'))
setenv('SUBJECTS_DIR',fullfile(base_dir, 'surfaceFreesurfer'));
% setenv('PERL5LIB','/Applications/freesurfer/mni/Library/Perl/Updates/5.10.0');
% setenv('PERL5LIB', '/Applications/freesurfer/mni/System/Library/Perl/5.8.6');

path1 = [path1 '/Applications/workbench/bin_macosx64'];
setenv('PATH', path1);

% defining the names of other directories
func_dir = 'func';
anat_dir = 'anat';
est_dir  = 'estimates';
fs_dir   = 'surfaceFreeSurfer';
wb_dir   = 'surfaceWB';

% list of subjects
subj_str = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06'};
subj_id  = [1, 2, 3, 4, 5, 6];
% list of runs within each session
%%% run_list{1} for session 1 and run_list{2} for session 2
run_list = {[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [13, 14, 15, 16, 17, 18]};

% AC coordinates
loc_AC = {[-103, -140, -140],...       %sub-01
          [-103, -133, -142],...       %sub-02
          [-101, -138, -130],...       %sub-03
          [-104, -147, -139],...       %sub-04
          [-100, -134, -134],...       %sub-05
          [-102, -134, -145],...       %sub-06
        };
    

sess = {'training', 'test'}; % training runs are considered to be ses-01 and testing runs are ses-02
ses_str = {'ses-01', 'ses-02'};
numTRs = 281;

hem     = {'L', 'R', 'cereb'};
hemName = {'CortexLeft', 'CortexRight', 'Cerebellum'};
% =========================================================================
numDummys = 0; % we need to make sure that this is correct
% based on comments here there should be 6 dummies: https://openneuro.org/datasets/ds002306/versions/1.1.0 
% based on the paper, there are 6 dummies: 3 dummies in the beginning of
% the run and 3 dummies at the end (281 - 6 = 275, the paper says 275 time points)
% BUT, using the times recorded in the tsv files, if we subtract the 3
% dummies in the beginning, the onsets become negative, no dummies are
% removed for now, as we think the times reported in the tsv files are
% including the dummies
% =========================================================================

% =========================================================================
% SUMMARY of the task design from the paper
% Each task has 12 instances: 8 (training) + 4 (testing)
% 3 days, 6 runs per day: 12 training runs + 6 test runs
% each run contains 77-83 trials lasting for 6-12 seconds
% Task order is pseudorandomized in the training runs! =
%("some tasks depended on each other and were therefore presented close to each other in time")
% In the test runs, 103 tasks were presented four times in the same order across all six runs
% Task order is fixed in the testing runs.
% Looking at the tsv files, the order during =training runs is the same
% across all participants.
% =========================================================================

switch what
    case 'ANAT:reslice_lpi'  % reslice anatomical to LPI
        % Example usage:nishimoto_imana('ANAT:reslice_lpi')
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
    case 'ANAT:center_ac'    % recenter to AC (manually retrieve coordinates)
        % Example usage: nishimoto_imana('ANAT:center_ac')
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
    case 'ANAT:segment'      % segment the anatomical image
        % also saves the bias field estimated in SPM
        % ********IF YOU WANT TO APPLY SPM BIAS CORRECTION TO, USE
        % ANAT:T1w_bcorrect CASE***********************
        % Example usage: nishimoto_imana('ANAT:segment')
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
    case 'ANAT:T1w_bcorrect' % bias correction for anatomical T1w (optional)
        % nishimoto_imana('ANAT:T1w_bcorrect')
        sn = subj_id;
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Bias correcting the anatomica image for %s\n', subj_str{s});
            % Get the directory of subjects anatomical and functional
            subj_anat_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            % make copy of original mean epi, and work on that
            % Get the name of the anatomical image
            source  = fullfile(subj_anat_dir, sprintf('%s_T1w_lpi.nii', subj_str{s}));
            dest    = fullfile(subj_anat_dir, sprintf('b%s_T1w_lpi.nii', subj_str{s}));
            copyfile(source, dest);
            
            % bias correct mean image for grey/white signal intensities
            P{1}    = dest;
            spmj_bias_correct(P);
        end % s (sn)
    case 'ANAT:bet'          % brain extraction for the anatomical image
        % Run bash script /srv/diedrichsen/shell/optiBET.sh
        % Edit command variable to set path to optiBET.sh script
        % Example usage: nishimoto_imana('ANAT:bet', 'sn', [1, 2, 3, 4, 5, 6])
        sn   = subj_id; % list of subjects
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);

            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w', subj_str{s});

            img    = fullfile(subj_dir, sprintf('%s_lpi.nii', anat_name));
            command = sprintf('bash /srv/diedrichsen/shell/optiBET.sh -i %s', img);
            system(command)
            
            in = fullfile(subj_dir, sprintf('%s_lpi_optiBET_brain.nii.gz', anat_name));
            out = fullfile(subj_dir, sprintf('%s_lpi_brain.nii.gz', anat_name));
            move_command = sprintf('mv %s %s', in, out);
            system(move_command);
            
            fprintf('optiBET completed for %s \n',subj_str{s})
            fprintf('Check the output of optiBET using FSLeyes or some other visualization software.')
           
        end   

    case 'FUNC:rename'           % removing dummies and renaming 
        % "WARNING" NO NEED TO REMOVE DUMMIES (WE THINK SO),SO THE REMOVING
        % PART IS COMMENTED OUT AND THIS CASE WILL ONLY BE USED FOR
        % RENAMING THE FILES.
        % Files are renamed based on an excel file which you can find in
        % the repository(MultTask_RunOrderInfo.xlsx)
        % removes dummies from the beginning of the functional images and
        % save the new images with new names so that we don't lose the
        % original images. For code efficiency, it will also rename the tsv
        % files. Runs done on the same day are considered as one session. 
        % this case calls FUNC:get_in_info to get the run information and
        % assign runs to sessions
        % Example usage: nishimoto_imana('FUNC:rename', 'sn', 1)
        sn = subj_id;
        
        vararginoptions(varargin, 'sn')
        
        for s = sn
            % go to subject's directory
            func_subj_dir = fullfile(base_dir, subj_str{s}, func_dir);
            cd(func_subj_dir)
            % unzip the files
            gunzip('*.gz');
            % loop over sessions
            for ss = [1, 2]
                % get the session name
%                 ses_name = sess{ss};
                % get the functional scans for the session
%                 fun_scans = dir(sprintf('%s*-%s*bold.nii', subj_str{s}, ses_name));
                fun_scans = dir(sprintf('%s_ses-02_run*.nii', subj_str{s}));
                
                
                % loop over functional scans and rename
                for i = 1:length(fun_scans)
                    fprintf('- Doing %s %s run %02d\n', subj_str{s}, ses_name, i);
                    src_filename = fun_scans(i).name;
                    des_filename = sprintf('%s_ses-%02d_run-%02d.nii', subj_str{s}, ss, run_list{ss}(i));
                    movefile(src_filename,fullfile(func_subj_dir, des_filename))
                    
                    % change the names of the tsv files
                    %%% get the run id
                    tsv_source = sprintf('%s_task-%s_run-%02d_events.tsv', subj_str{s}, ses_name, i);
                    tsv_dest = sprintf('%s_ses-%02d_run-%02d_events.tsv', subj_str{s}, ss, run_list{ss}(i));
                    movefile(fullfile(func_subj_dir, tsv_source), fullfile(func_subj_dir, tsv_dest))
                end % i (runs within a scan)
            end % ss (session)
        end % sn (subjects)    
    case 'FUNC:realign'          % realign functional images
        % SPM realigns all volumes to the first volume of first run
        % example usage: nishimoto_imana('FUNC:realign', 'sn', 1)
        
        sn   = subj_id; % list of subjects
        
        vararginoptions(varargin, {'sn'});
                
        for s = sn
            func_subj_dir = fullfile(base_dir, subj_str{s}, func_dir);
            % cd to the folder with raw functional data
            cd(func_subj_dir)
            spm_jobman('initcfg')
            
            data = {}; % initialize data cell array which will contain file names for runs/TR images
            for ses = [1, 2]
                runs = run_list{ses};
                for r = runs
                    for j = 1:numTRs-numDummys
                        data{r}{j,1} = sprintf('%s_ses-%02d_run-%02d.nii,%d', subj_str{s},ses, r,j);
                    end % j (TRs/images)
                end % r (runs)
            end
            spmj_realign(data);
            fprintf('- runs realigned for %s  ses %02d\n',subj_str{s}, ses);
        end % s (sn)
    case 'FUNC:meanepi_bcorrect' % bias correction for the mean image before coreg (optional)
        % uses the bias field estimated in SPM segmenttion
        % Example usage: nishimoto_imana('FUNC:meanepi_bcorrect')
        sn = subj_id;
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Bias correcting mean epi for %s\n', subj_str{s});
            % Get the directory of subjects anatomical and functional
            subj_func_dir = fullfile(base_dir, subj_str{s}, func_dir);
            % make copy of original mean epi, and work on that
            source  = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01.nii', subj_str{s}));
            dest    = fullfile(subj_func_dir, sprintf('bmean%s_ses-01_run-01.nii', subj_str{s}));
            copyfile(source, dest);
            
            % bias correct mean image for grey/white signal intensities
            P{1}    = dest;
            spmj_bias_correct(P);
        end % s (sn)
    case 'FUNC:coreg'            % coregistration with the anatomicals using spm
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and mean functional image to overlay
        % - Manually adjust mean functional image and save the results ("r" will be added as a prefix)
        % Example usage: nishimoto_imana('FUNC:coreg', 'sn', [2], 'prefix', 'r')
        
        sn       = subj_id;   % list of subjects        
        step     = 'manual';  % first 'manual' then 'auto'
        prefix   = 'rb';      % to use the bias corrected version, set it to 'rbb'
        % ===================
        % After the manual registration, the mean functional image will be
        % saved with r as the prefix which will then be used in the
        % automatic registration
        vararginoptions(varargin, {'sn', 'step', 'prefix'});
        spm_jobman('initcfg')
        for s = sn
            % Get the directory of subjects anatomical and functional
            subj_anat_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            subj_func_dir = fullfile(base_dir, subj_str{s}, func_dir);
            
            cd(subj_anat_dir); % goes to subjects anatomical dir so that coreg tool starts from that directory (just for convenience)
            
            switch step
                case 'manual'
                    coregtool;
                    keyboard;
                case 'auto'
                    % do nothing
            end % switch step
            
            % (2) Automatically co-register functional and anatomical images
            J.ref = {fullfile(subj_anat_dir, sprintf('%s_T1w_lpi.nii', subj_str{s}))}; % just one anatomical or more than one?
            
            J.source = {fullfile(subj_func_dir, sprintf('%smean%s_ses-01_run-01.nii', prefix, subj_str{s}))};
            
            J.other             = {''};
            J.eoptions.cost_fun = 'nmi';
            J.eoptions.sep      = [4 2];
            J.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm     = [7 7];
            
            % 
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);
            
            % (3) Manually check again
%             coregtool;
%             keyboard();
            
            % NOTE:
            % Overwrites meanepi, unless you update in step one, which saves it
            % as rmeanepi.
            % Each time you click "update" in coregtool, it saves current
            % alignment by appending the prefix 'r' to the current file
            % So if you continually update rmeanepi, you'll end up with a file
            % called r...rrrmeanepi.
        end % s (sn) 
    case 'FUNC:make_samealign'   % align all the functionals
        % Aligns all functional images to rmean functional image
        % Example usage: nishimoto_imana('FUNC:make_samealign', 'prefix', 'r', 'sn', [6])
        
        sn     = subj_id;     % subject list
        prefix = 'r';         % prefix for the meanepi: r or rbb if bias corrected
        
        vararginoptions(varargin, {'sn', 'prefix'});
                
        for s = sn
            % Get the directory of subjects functional
            subj_func_dir = fullfile(base_dir, subj_str{s}, func_dir);
            
            % get the images from both sessions
            for ss = [1, 2]
                % get the list of runs for the session
                runs = run_list{ss};
                fprintf('- make_samealign  %s \n', subj_str{s})
                cd(subj_func_dir);
                
                % Select image for reference 
                %%% note that functional images are aligned with the first
                %%% run from first session hence, the ref is always rmean<subj>_ses-01_run-01
                P{1} = fullfile(subj_func_dir, sprintf('%smean%s_ses-01_run-01.nii', prefix,subj_str{s}));
%                 ref_image = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01_func2highres.nii',subj_str{s}));
%                 P{1} = ref_image
                % Select images to be realigned
                Q = {};
                for r = runs
                    for i = 1:numTRs - numDummys
                        Q{end+1}    = fullfile(subj_func_dir,...
                                               sprintf('%s%s_ses-%02d_run-%02d.nii,%d', prefix, subj_str{s}, ss, r, i));
                    end
                end % r(runs)
                
                spmj_makesamealign_nifti(char(P),char(Q));
            end % ss (sess)
        end % s (sn)
    case 'FUNC:make_maskImage'   % make mask images (noskull and grey_only)
        % Make maskImage in functional space
        % run this step after GLM!
        % Example usage: nishimoto_imana('FUNC:make_maskImage', 'prefix', 'r', 'sn', [6])
        
        sn     = subj_id; % list of subjects
        prefix = 'r';     % prefix for the meanepi: r or rbb if bias corrected
        vararginoptions(varargin, {'sn', 'prefix'});
        
        
        for s = sn
            % Get the directory of subjects anatomical and functional
            subj_anat_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            subj_func_dir = fullfile(base_dir, subj_str{s}, func_dir);

            fprintf('- make mask for %s\n', subj_str{s});
            cd(subj_func_dir);
            
            nam{1}  = fullfile(subj_func_dir, sprintf('%smean%s_ses-01_run-01.nii', prefix, subj_str{s}));
            nam{2}  = fullfile(subj_anat_dir, sprintf('c1%s_T1w_lpi.nii', subj_str{s}));
            nam{3}  = fullfile(subj_anat_dir, sprintf('c2%s_T1w_lpi.nii', subj_str{s}));
            nam{4}  = fullfile(subj_anat_dir, sprintf('c3%s_T1w_lpi.nii', subj_str{s}));
            spm_imcalc(nam, 'rmask_noskull.nii', 'i1>0 & (i2+i3+i4)>0.1')
            
            nam     = {};
            nam{1}  = fullfile(subj_func_dir, sprintf('%smean%s_ses-01_run-01.nii', prefix, subj_str{s}));
            nam{2}  = fullfile(subj_anat_dir, sprintf('c1%s_T1w_lpi.nii', subj_str{s}));
            spm_imcalc(nam, 'rmask_gray.nii', 'i1>0 & i2>0.1')
            
            nam     = {};
            nam{1}  = fullfile(subj_func_dir, sprintf('%smean%s_ses-01_run-01.nii', prefix, subj_str{s}));
            nam{2}  = fullfile(subj_anat_dir, sprintf('c1%s_T1w_lpi.nii', subj_str{s}));
            nam{3}  = fullfile(subj_anat_dir, sprintf('c5%s_T1w_lpi.nii', subj_str{s}));
            spm_imcalc(nam, 'rmask_grayEyes.nii', 'i1>2400 & i2+i3>0.4')
            
            nam     = {};
            nam{1}  = fullfile(subj_func_dir, sprintf('%smean%s_ses-01_run-01.nii', prefix, subj_str{s}));
            nam{2}  = fullfile(subj_anat_dir, sprintf('c5%s_T1w_lpi.nii', subj_str{s}));
            nam{3}  = fullfile(subj_anat_dir, sprintf('c1%s_T1w_lpi.nii', subj_str{s}));
            nam{4}  = fullfile(subj_anat_dir, sprintf('c2%s_T1w_lpi.nii', subj_str{s}));
            nam{5}  = fullfile(subj_anat_dir, sprintf('c3%s_T1w_lpi.nii', subj_str{s}));
            spm_imcalc(nam, 'rmask_noskullEyes.nii', 'i1>2000 & (i2+i3+i4+i5)>0.2')

        end % s (sn)
    case 'FUNC:run'              % add functional pipelines here
        % Example usage: nishimoto_imana('FUNC:run', 'sn', [2, 3, 4, 5, 6])
        
        sn  = subj_id;        
        vararginoptions(varargin, {'sn'});
        
        nishimoto_imana('FUNC:rename', 'sn', sn)
        nishimoto_imana('FUNC:realign', 'sn', sn);
%         nishimoto_imana('FUNC:coreg_fsl', 'sn', sn, 'prefix', 'r')
%         nishimoto_imana('FUNC:make_samealign', 'prefix', 'r', 'sn', sn);
%         nishimoto_imana('FUNC:make_maskImage', 'prefix', 'r', 'sn', sn);            
        
    case 'GLM:task_info'    % creates a text file and assign numbers to the tasks/conditions
        % Example usage: nishimoto_imana('GLM:task_info')
        
        info_dir = fullfile(base_dir, 'sub-02', func_dir);
        
        % loop over sessions/runs and load the info tsv file for each run
        TN_cell = {};
        info_struct = [];
        
        for r = run_list{1}
            % load the tsv file
            tsv_file = fullfile(info_dir, sprintf('sub-02_ses-01_run-%02d_events.tsv', r));
            T = dload(tsv_file);

            % get task names:
            TN_cell = [TN_cell; T.trial_type];
        end % r (runs)

        % get the unique task names
        task_name = unique(TN_cell)';

        % assign numbers to each task
        task = 1:length(task_name);

        % create a structure with the info you want
        info_tmp.task_name = task_name';
        info_tmp.task = task';

        info_struct = addstruct(info_struct, info_tmp);
        
        % save as a tsv file 
        dsave(fullfile(base_dir, 'tasks_info.tsv'), info_struct);        
    case 'GLM:design1_ses'  % make the design matrix for the glm
        % models each condition as a separate regressors
        % For conditions with multiple repetitions, one regressor
        % represents all the instances
        % nishimoto_imana('GLM:design1', 'sn', [6])
        
        sn = subj_id;
        hrf_cutoff = Inf;
        ses = 1;
        prefix = 'r'; % prefix of the preprocessed epi we want to use
        glm = 1;
        vararginoptions(varargin, {'sn', 'hrf_cutoff', 'ses'});
        
        
        % get the info file that specifies the order of the tasks
        Dd = dload(fullfile(base_dir, 'tasks_info.tsv'));
        
        for s = sn
            func_subj_dir = fullfile(base_dir, subj_str{s}, func_dir);
            % loop over runs
%             for ss = [1, 2]
                % create a directory to save the design
                subj_est_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), ses_str{ses});
                dircheck(subj_est_dir)
                
                T = []; % task/condition + session + run info
                J = []; % structure with SPM fields to make the design
                
                J.dir            = {subj_est_dir};
                J.timing.units   = 'secs';
                J.timing.RT      = 2.0;
                J.timing.fmri_t  = 16;
                J.timing.fmri_t0 = 1;
                
                % get the list of runs for the current session
                runs = run_list{ses};
                
                % loop through runs within the current sessions
                itaskUni = 0;
                for run = 1:length(runs)
                    
%                   % fill in nifti image names for the current run
                    N = cell(numTRs - numDummys, 1); % preallocating!
                    for i = 1:(numTRs-numDummys)
                        N{i} = fullfile(func_subj_dir, sprintf('%s%s_ses-%02d_run-%02d.nii, %d', prefix, subj_str{s}, ses, runs(run), i));
                    end % i (image numbers)
                    J.sess(run).scans = N; % scans in the current runs
                    
                    % get the path to the tsv file
                    tsv_path = fullfile(base_dir, subj_str{s}, func_dir);
                    % get the tsvfile for the current run
                    D = dload(fullfile(tsv_path, sprintf('%s_ses-%02d_run-%02d_events.tsv', subj_str{s}, ses, runs(run))));
                                        
                    % loop over trials within the current run and build up
                    % the design matrix
                    for ic = 1:length(Dd.task_name)
                        itaskUni = itaskUni+1;
                        % get the indices corresponding to the current
                        % condition.
                        % this line is necessary as there are some
                        % conditions with more than 1 repetition
                        idx = strcmp(D.trial_type, Dd.task_name{ic});
                        fprintf('* %d instances found for condition %s in run %02d\n', sum(idx), Dd.task_name{ic}, runs(run))
                        
                        %
                        % filling in "reginfo"
                        TT.sn        = s;
                        TT.sess      = ses;
                        TT.run       = runs(run);
                        TT.task_name = Dd.task_name(ic);
                        TT.task      = ic;
                        TT.taskUni   = itaskUni;
                        TT.n_rep     = sum(idx);
                        
                        % filling in fields of J (SPM Job)
                        J.sess(run).cond(ic).name = Dd.task_name{ic};
                        J.sess(run).cond(ic).tmod = 0;
                        J.sess(run).cond(ic).orth = 0;
                        J.sess(run).cond(ic).pmod = struct('name', {}, 'param', {}, 'poly', {});
                        
                        % get onset and duration (should be in seconds)
                        onset    = D.onset(idx) - (J.timing.RT*numDummys);
                        fprintf("The onset is %f\n", onset)
                        if onset < 0
                            warning("negative onset found")
                        end
                        duration = D.duration(idx);
                        fprintf("The duration is %f\n", duration);
                        
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
                
                J.fact             = struct('name', {}, 'levels', {});
                J.bases.hrf.derivs = [0 0];
                J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
                J.volt             = 1;
                J.global           = 'None';
                J.mask             = {fullfile(func_subj_dir,'rmask_noskull.nii,1')};
                J.mthresh          = 0.05;
                J.cvi_mask         = {fullfile(func_subj_dir,'rmask_gray.nii')};
                J.cvi              =  'fast';
                
%                 spm_rwls_run_fmri_spec(J);
                
                dsave(fullfile(J.dir{1},sprintf('%s_ses-%02d_reginfo.tsv', subj_str{s}, ses)), T);
                fprintf('- estimates for glm_%d session %d has been saved for %s \n', glm, ses, subj_str{s});
%             end % ss (session)
            
            
        end % sn (subject)      
    case 'GLM:design2_all'  % make the design matrix for the glm
        % models each condition as a separate regressors
        % For conditions with multiple repetitions, one regressor
        % represents all the instances
        % nishimoto_imana('GLM:design1', 'sn', [6])
        
        sn = subj_id;
        hrf_cutoff = Inf;
        prefix = 'r'; % prefix of the preprocessed epi we want to use
        glm = 1;
        vararginoptions(varargin, {'sn', 'hrf_cutoff', 'ses'});
        
        
        % get the info file that specifies the order of the tasks
        Dd = dload(fullfile(base_dir, 'tasks_info.tsv'));
        
        for s = sn
                func_subj_dir = fullfile(base_dir, subj_str{s}, func_dir);
                % loop over runs
                % create a directory to save the design
                subj_est_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), ses_str{ses});
                dircheck(subj_est_dir)
                
                T = []; % task/condition + session + run info
                J = []; % structure with SPM fields to make the design
                
                J.dir            = {subj_est_dir};
                J.timing.units   = 'secs';
                J.timing.RT      = 2.0;
                J.timing.fmri_t  = 16;
                J.timing.fmri_t0 = 1;
                
                % loop through runs within the current sessions
                itaskUni = 0;
                for ses = [1, 2]
                    % get the list of runs for the current session
                    runs = run_list{ses};
                    for run = 1:length(runs)
                        
                        % fill in nifti image names for the current run
                        N = cell(numTRs - numDummys, 1); % preallocating!
                        for i = 1:(numTRs-numDummys)
                            N{i} = fullfile(func_subj_dir, sprintf('%s%s_ses-%02d_run-%02d.nii, %d', prefix, subj_str{s}, ses, run, i));
                        end % i (image numbers)
                        J.sess(run).scans = N; % scans in the current runs
                        
                        % get the path to the tsv file
                        tsv_path = fullfile(base_dir, subj_str{s}, func_dir);
                        % get the tsvfile for the current run
                        D = dload(fullfile(tsv_path, sprintf('%s_ses-%02d_run-%02d_events.tsv', subj_str{s}, ses, run)));
                        
                        % loop over trials within the current run and build up
                        % the design matrix
                        for ic = 1:length(Dd.task_name)
                            itaskUni = itaskUni+1;
                            % get the indices corresponding to the current
                            % condition.
                            % this line is necessary as there are some
                            % conditions with more than 1 repetition
                            idx = strcmp(D.trial_type, Dd.task_name{ic});
                            fprintf('* %d instances found for condition %s in run %02d\n', sum(idx), Dd.task_name{ic}, run)
                            
                            %
                            % filling in "reginfo"
                            TT.sn        = s;
                            TT.sess      = ses;
                            TT.run       = run;
                            TT.task_name = Dd.task_name(ic);
                            TT.task      = ic;
                            TT.taskUni   = itaskUni;
                            TT.n_rep     = sum(idx);
                            
                            % filling in fields of J (SPM Job)
                            J.sess(run).cond(ic).name = Dd.task_name{ic};
                            J.sess(run).cond(ic).tmod = 0;
                            J.sess(run).cond(ic).orth = 0;
                            J.sess(run).cond(ic).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            
                            % get onset and duration (should be in seconds)
                            onset    = D.onset(idx) - (J.timing.RT*numDummys);
                            fprintf("The onset is %f\n", onset)
                            if onset < 0
                                warning("negative onset found")
                            end
                            duration = D.duration(idx);
                            fprintf("The duration is %f\n", duration);
                            
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
                end % session
                
                
                J.fact             = struct('name', {}, 'levels', {});
                J.bases.hrf.derivs = [0 0];
                J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
                J.volt             = 1;
                J.global           = 'None';
                J.mask             = {fullfile(func_subj_dir,'rmask_noskull.nii,1')};
                J.mthresh          = 0.05;
                J.cvi_mask         = {fullfile(func_subj_dir,'rmask_gray.nii')};
                J.cvi              =  'fast';
                
                spm_rwls_run_fmri_spec(J);
                
                dsave(fullfile(J.dir{1},sprintf('%s_ses-%02d_reginfo.tsv', subj_str{s}, ses)), T);
                fprintf('- estimates for glm_%d session %d has been saved for %s \n', glm, ses, subj_str{s});
%             end % ss (session)
            
            
        end % sn (subject)      
    case 'GLM:check_design' % checking the design matrix
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
    case 'GLM:estimate'     % estimate beta values
        % Example usage: nishimoto_imana('GLM:estimate', 'glm', 1, 'ses', 1, 'sn', 6)
        
        sn       = subj_id; % subject list
        glm      = 1;       % glm number
        ses      = 1;       % session number
        
        vararginoptions(varargin, {'sn', 'glm', 'ses'})
        
        for s = sn
            fprintf('- Doing glm estimation for glm %02d %s\n', glm, subj_str{s});
            subj_est_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), sprintf('ses-%02d', ses));         
            
            load(fullfile(subj_est_dir,'SPM.mat'));
            SPM.swd = subj_est_dir;
            
            spm_rwls_spm(SPM);
        end % s (sn),
    case 'GLM:T_contrast'   % make T contrasts for each condition
        %%% Calculating contrast images.
        % Example usage: nishimoto_imana('GLM:T_contrast', 'sn', 2, 'glm', 1, 'ses', 1, 'baseline', 'rest')
        
        sn             = subj_id;    % subjects list
        ses            = 1;              % task number
        glm            = 1;              % glm number
        baseline       = 'rest';         % contrast will be calculated against base (available options: 'rest')
        
        vararginoptions(varargin, {'sn', 'glm', 'ses', 'baseline'})
        
        for s = sn
            
            % get the subject id folder name
            fprintf('Contrasts for session %02d %s\n', ses, subj_str{s})
            glm_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), ses_str{ses}); 
            
            cd(glm_dir);
            
            % load the SPM.mat file
            load(fullfile(glm_dir, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = dload(fullfile(glm_dir, sprintf('%s_%s_reginfo.tsv', subj_str{s}, ses_str{ses})));
            
            % t contrast for each condition type
            utask = unique(T.task)';
            idx = 1;
            for ic = utask
                switch baseline
                    case 'myBase' % contrast vs future baseline :)))
                        % put your new contrasts here!
                    case 'rest' % contrast against rest
                        con                          = zeros(1,size(SPM.xX.X,2));
                        con(:,logical((T.task == ic)& (T.n_rep>0))) = 1;
%                         n_rep = length(T.run(T.task == ic));
%                         n_rep_t = T.n_rep(T.task == ic);
%                         name = unique(T.task_name(T.task == ic));
%                         fprintf('- task is %s: \n', name{1});
%                         fprintf('number of reps in all runs = %d\n', n_rep);
%                         fprintf('numberof reps recorded in tsv = %d\n', n_rep_t);
                        con                          = con/abs(sum(con));            
                end % switch base

                % set the name of the contrast
                contrast_name = sprintf('%s-%s', char(unique(T.task_name(T.task == ic))), baseline);
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
        end % sn
    case 'GLM:F_contrast'   % make F contrast
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
    case 'GLM:check'        % visually inspect design matrix
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
    case 'GLM:run'          % add glm routines you want to run as pipeline
        % Example usage: nishimoto_imana('GLM:run', 'sn', [3, 4, 5, 6], 'glm', 1, 'ses', 1)
        
        sn  = subj_id; % subject id
        ses = 1;       % which dataset? ses1 or ses2
        glm = 1;       % number assigned to the glm
        
        vararginoptions(varargin, {'sn', 'ses', 'glm'});
        
        nishimoto_imana('GLM:design1', 'sn', sn, 'ses', ses);
        nishimoto_imana('GLM:estimate', 'sn', sn, 'glm', glm, 'ses', ses);
%         nishimoto_imana('GLM:F_contrast', 'sn', sn, 'glm', glm, 'ses', ses)
        nishimoto_imana('GLM:T_contrast', 'sn', sn, 'glm', glm, 'ses', ses, 'baseline', 'rest')
            
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
        % example: nishimoto_imana('SURF:xhemireg', 'sn', [1, 2, 3, 4, 5])
        
        sn   = subj_id; % list of subjects

        vararginoptions(varargin, {'sn'})
        
        % set freesurfer directory
        fs_dir = fullfile(base_dir, 'surfaceFreeSurfer');
        
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
        % Example usage: nishimoto_imana('SURF:map_ico', 'sn', [1, 2, 3, 4, 5, 6])
        
        sn = subj_id; % list of subjects
        
        vararginoptions(varargin, {'sn'});
        
        % set freesurfer directory
        fs_dir = fullfile(base_dir, 'surfaceFreeSurfer');
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
        fs_dir = fullfile(base_dir, 'surfaceFreeSurfer');
        
        for s = sn 
            fprintf('- fs2wb %s\n', subj_str{s});
            wb_subj_dir  = fullfile(base_dir, wb_dir, 'data', subj_str{s});
            surf_resliceFS2WB(subj_str{s}, fs_dir, wb_subj_dir, 'hemisphere', hemi, 'resolution', sprintf('%dk', res))
        end % s (sn)
    case 'SURF:run_all'        % Pipeline running all of surface preprocessing routines
        % Example usage: nishimoto_imana('SURF:run_all')
        
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
%         nishimoto_imana('SURF:reconall')
        nishimoto_imana('SURF:xhemireg', 'sn', sn);
        nishimoto_imana('SURF:map_ico', 'sn', sn);
        nishimoto_imana('SURF:fs2wb', 'sn', sn);
    case 'SURF:vol2surf'       % Mapping volumetric data to surface
        % first univariately whiten data and then map to surface
        % Example usage: nishimoto_imana('SURF:vol2surf', 'sn', 1, 'group', 0)
        
        sn             = subj_id;        % subjects list
        ses            = 1;              % task number
        glm            = 1;              % glm number
        type           = 'con';          % type of data to be mapped. Options are: 'con', 'beta'
        baseline       = 'rest';         % contrast will be calculated against base (available options: 'rest')
        kernel         = 1;              % smoothing kernel 
        group          = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'ses', 'baseline', 'type', 'kernel', 'group'})
        
        % loop over hemispheres
        for h = 1:2
            for s = sn
                
                % get directory where subject surfaces are saved
                surf_dir = fullfile(base_dir, wb_dir, sprintf('glm%02d', glm), subj_str{s});
                dircheck(surf_dir);
                glm_dir = fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), sprintf('ses-%02d', ses));
                
                % file containing the giftis for the current subject
                filename{s} = fullfile(surf_dir, sprintf('%s.w%s-%s.ses-%02d.%s.func.gii', subj_str{s}, type, baseline, ses, hem{h}));
                
                white   = fullfile(base_dir, wb_dir, 'data', subj_str{s}, sprintf('%s.%s.white.32k.surf.gii', subj_str{s}, hem{h}));
                pial    = fullfile(base_dir, wb_dir, 'data', subj_str{s}, sprintf('%s.%s.pial.32k.surf.gii', subj_str{s}, hem{h}));
                C1      = gifti(white);
                C2      = gifti(pial);% map the files
                
                switch type
                    case 'con'
                        % get all the contrasts
                        files2map = dir(fullfile(glm_dir, sprintf('con_*-%s.nii', baseline)));
                    case 'beta'
                        files = dir(fullfile(base_dir, subj_str{s}, est_dir, sprintf('glm%02d', glm), sprintf('ses-%02d', ses), 'beta*.nii'));
                        files2map = cat(1, {files(:).name});
                        names = string(1:length(files2map));
                end % switch type
                %%% get file paths and names of contrasts
                for i= 1:length(files2map)
                    name{i} = fullfile(glm_dir, files2map(i).name);
                    column_names{i} = files2map(i).name(5:end-4);
                end % i (contrast names)
                
                
                % % if the subj specific gifti files have not been created, then
                % create them
                if group ~=1
                    fprintf('- Transforming %s for %s hemi %s\n', type, subj_str{s}, hem{h});
                    maps = surf_vol2surf(C1.vertices, C2.vertices, name, 'column_names', column_names, ...
                        'anatomicalStruct', hemName{h});
                    
                    % map ResMS (will be used for univariate prewhitening)
                    Gres = surf_vol2surf(C1.vertices, C2.vertices, {fullfile(glm_dir, 'ResMS.nii')}, ...
                        'anatomicalStruct', hemName{h});
                    
                    % do univariate prewhitening
                    data    = bsxfun(@rdivide, maps.cdata, Gres.cdata);
                    
                    % create one single gifti file containing all the contrasts
                    data(:, length(files2map)+1) = Gres.cdata;
                    column_names{length(files2map)+1} = 'ResMS';
                    G = surf_makeFuncGifti(data,'anatomicalStruct', hemName{h}, 'columnNames', column_names);
                    
                    % save the single gifti file
                    %%% a cell array containing the filenames.
                    %%% will be used in creating group maps and summary
                    save(G, filename{s});
                    fprintf('- Done %s %s map2surf\n', subj_str{s}, type);
                    % smooth the single gifti
                    atlas_dir = fullfile(base_dir, 'fs_LR_32');
                    atlas_file = fullfile(atlas_dir, sprintf('fs_LR.32k.%s.inflated.surf.gii',hem{h}));
                    surf_smooth(filename{s},'surf',atlas_file,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                end % if not group
            end % s (subjects)
            
            % create group and group summary
            if group == 1
                wb_group_dir = fullfile(base_dir, wb_dir, sprintf('glm%02d', glm), 'group');dircheck(wb_group_dir);
                cd(wb_group_dir)
                summaryname     = fullfile(wb_group_dir,sprintf('wgroup.%s-%s.ses-%02d.glm%02d.%s.func.gii', type, baseline, ses, glm, hem{h}));
                surf_groupGiftis(filename, 'groupsummary', summaryname, 'replaceNaNs', 1, 'outcolnames', subj_str, 'outfilenamePattern', ['%s.', hem{h}, '.func.gii']);
            end % if group
        end % hemisphere
        
               
    case 'SUIT:isolate_segment'    % Segment cerebellum into grey and white matter
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
    case 'SUIT:normalise_dartel'   % SUIT normalization using dartel
        % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
        % example usage: nishimoto_imana('SUIT:normalise_dartel')
        sn = subj_id; %subjNum
        vararginoptions(varargin, 'sn');
        
        for s = sn
            suit_subj_dir = fullfile(base_dir, subj_str{s}, 'suit');
            mkdir(suit_subj_dir)
            
            cd(suit_subj_dir)
            job.subjND.gray       = {fullfile(suit_subj_dir, sprintf('c_%s_T1w_lpi_seg1.nii', subj_str{s}))};
            job.subjND.white      = {fullfile(suit_subj_dir, sprintf('c_%s_T1w_lpi_seg2.nii', subj_str{s}))};
            job.subjND.isolation  = {fullfile(suit_subj_dir, sprintf('c_%s_T1w_lpi_pcereb_corr.nii', subj_str{s}))};
            suit_normalize_dartel(job);
        end % s (subjects)    

    case 'SUIT:save_dartel_def'    
        % Saves the dartel flow field as a deformation file. 
        % example usage: nishimoto_imana('SUIT:save_dartel_def')
        sn = subj_id; %subjNum
        vararginoptions(varargin, 'sn');

        for s = sn
            suit_subj_dir = fullfile(base_dir, subj_str{s}, 'suit', 'anat');

            cd(suit_subj_dir);
            anat_name = sprintf('%s_T1w_lpi', subj_str{s});
            suit_save_darteldef(anat_name);
        end % s (subjects)

    case 'SUIT:reslice'            % Reslice stuff into suit space 
        % run the case with 'anatomical' to check the suit normalization
        % Example usage: nishimoto_imana('SUIT:reslice','type','ResMS', 'mask', 'pcereb_corr')
        % make sure that you reslice into 2mm^3 resolution
        
        sn   = subj_id;
        type = 'anatomical';  % 'betas' or 'con' or 'ResMS' or 'cerebellarGrey' or 'anatomical'
        mask = 'pcereb_corr'; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask' or 'pcereb'
        glm  = 1;             % glm number. Used for reslicing betas and contrasts 
        
        vararginoptions(varargin, {'sn', 'type', 'mask', 'glm'})
        
        for s = sn
            suit_dir = fullfile(base_dir, subj_str{s}, 'suit', 'anat');
            switch type
                case 'anatomical'
                    subj_dir = suit_dir;
                    % Get the name of the anatpmical image
                    files2map = sprintf('%s_T1w_lpi.nii', subj_str{s});
                    
                    job.subj.resample = {sprintf('%s,1', files2map)};                 
                case 'con'
                    subj_dir     = fullfile(base_dir, subj_str{s}, 'estimates', sprintf('glm%02d', glm), 'ses-01');
                    out_dir      = fullfile(base_dir, subj_str{s}, 'suit',sprintf('glm%02d',glm));
                    files2map    = dir(fullfile(subj_dir,sprintf('*con*'))); % find images to be resliced
                    
                    job.subj.resample     = {files2map.name};
                case 'ResMS'
                    subj_dir     = fullfile(base_dir, subj_str{s}, 'estimates', sprintf('glm%02d', glm), 'ses-01');
                    out_dir      = fullfile(base_dir, subj_str{s}, 'suit',sprintf('glm%02d',glm));
                    files2map    = dir(fullfile(subj_dir,sprintf('*ResMS*'))); % find images to be resliced
                    
                    job.subj.resample     = {files2map.name};
                    
            end
            
            cd(subj_dir);
            job.subj.affineTr  = {fullfile(suit_dir,sprintf('Affine_c_%s_T1w_lpi_seg1.mat', subj_str{s}))};
            job.subj.flowfield = {fullfile(suit_dir,sprintf('u_a_c_%s_T1w_lpi_seg1.nii', subj_str{s}))};
            job.subj.mask      = {fullfile(suit_dir, sprintf('c_%s_T1w_lpi_%s.nii', subj_str{s}, mask))};
            job.vox            = [2 2 2];
            suit_reslice_dartel(job);
            
            if ~strcmp(type,'anatomical')
                source=fullfile(subj_dir, '*wd*');
                dircheck(fullfile(out_dir));
                movefile(source,out_dir);
            end
            fprintf('- Resliced %s %s to SUIT\n', subj_str{s}, type);
        end % s (subject)
    case 'SUIT:map2flat'           % Creates flatmaps
        % this case also creates group average for each task
        % First RUN nishimoto_imana('SUIT:map2flat', 'group', 0) to get the
        % gifti file for all the contrasts for each subject
        % Then RUN nishimoto_imana('SUIT:map2flat', 'group', 1) to get the
        % group summaries
        % Example usage: nishimoto_imana('SUIT:map2flat')
        
        sn    = subj_id;
        glm   = 1;
        type  = 'con'; % type of the image to be mapped to flatmap
        baseline = 'rest';
        group = 1;     % if this flag is set to 1, it bypasses the step that creates gifti files for each subject
        
        vararginoptions(varargin, {'sn', 'glm', 'type', 'group', 'baseline'});
        
        for s = sn
            suit_dir = fullfile(base_dir, subj_str{s}, 'suit', sprintf('glm%02d', glm));
            % file containing the giftis for the current subject
            filename{s} = fullfile(suit_dir, sprintf('%s.w%s.cerebellum.func.gii', subj_str{s}, type));
            
            switch type
                case 'con'
                    % get all the contrasts
                    files2map = dir(fullfile(suit_dir, sprintf('wd%s*%s', type, baseline)));
            end
            
            % map the files
            %%% get file paths and names of contrasts
            for i= 1:length(files2map)
                name{i} = fullfile(suit_dir, files2map(i).name);
                column_names{i} = files2map(i).name(7:end-4);
            end % i (contrast names)
            
            % if the subj specific gifti files have not been created, then
            % create them
            if group ~=1
                maps = suit_map2surf(name, 'stats', 'nanmean' );
                        
                % map ResMS (will be used for univariate prewhitening)
                mapResMS = suit_map2surf(fullfile(suit_dir, 'wdResMS.nii'), 'stats', 'nanmean');
                Gres = surf_makeFuncGifti(mapResMS,'anatomicalStruct', 'Cerebellum', 'columnNames', {'ResMS'});
                
                % do univariate prewhitening
                data    = bsxfun(@rdivide, maps, mapResMS);
                
                % create one single gifti file containing all the contrasts
                data(:, length(files2map)+1) = Gres.cdata;
                column_names{length(files2map)+1} = 'ResMS';
                G = surf_makeFuncGifti(data,'anatomicalStruct', 'Cerebellum', 'columnNames', column_names);
                
                % save the single gifti file
                %%% a cell array containing the filenames.
                %%% will be used in creating group maps and summary
                save(G, filename{s});
                fprintf('- Done %s %s map2flat\n', subj_str{s}, type);                
            end % if not group  
        end % s (subject)
        
        % create group and group summary
        if group == 1
            suit_group_dir = fullfile(base_dir, 'suit', sprintf('glm%02d', glm), 'group');dircheck(suit_group_dir);
            cd(suit_group_dir)
            summaryname     = fullfile(suit_group_dir,sprintf('wgroup.%s.glm%02d.func.gii', type, glm));
            surf_groupGiftis(filename, 'groupsummary', summaryname, 'replaceNaNs', 1, 'outcolnames', subj_str);
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

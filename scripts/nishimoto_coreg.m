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

base_dir = sprintf('%s/FunctionalFusion/Nishimoto_103Task/',workdir);
% base_dir = '/Users/ladan/Documents/DATA/nishimoto';

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
numTRs = 281;
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
% Task order is pseudorandomized in the training runs! 
%("some tasks depended on each other and were therefore presented close to each other in time")
% In the test runs, 103 tasks were presented four times in the same order across all six runs
% Task order is fixed in the testing runs.
% Looking at the tsv files, the order during training runs is the same
% across all participants.
% =========================================================================

switch what
    case 'ANAT:bet'             % Brain extraction for anatomical.nii
        % Run bash script /srv/diedrichsen/shell/optiBET.sh
        % Edit command variable to set path to optiBET.sh script
        % example: bsp_imana('ANAT:bet',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs,
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);

            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w', subj_str{s});

            img    = fullfile(subj_dir, sprintf('%s_lpi.nii', anat_name))
            command = sprintf('bash /srv/diedrichsen/shell/optiBET.sh -i %s', img)
            system(command)
            
            in = fullfile(subj_dir, sprintf('%s_lpi_optiBET_brain.nii.gz', anat_name));
            out = fullfile(subj_dir, sprintf('%s_lpi_brain.nii.gz', anat_name));
            copy_command = sprintf('cp %s %s', in, out)
            system(copy_command)
            
            fprintf('optiBET completed for %s \n',subj_name{sn(s)})
            fprintf('Check the output of optiBET using FSLeyes or some other visualization software.')
           
        end
    case 'FUNC:coreg_fsl'        % Coregister meanepi to anatomical with fsl
        % Need meanrun_01 and brain-extracted anatomical
        % example: bsp_imana('FUNC:coreg_meanepi_fsl',1,8)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs,
            
            % Get the directory of subjects anatomical and functional
            subj_anat_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            subj_func_dir = fullfile(base_dir, subj_str{s}, func_dir);

            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w', subj_str{s});

            
            t1 = fullfile(subj_dir, sprintf('%s_lpi.nii', anat_name));
            t1_brain = fullfile(subj_dir, sprintf('%s_lpi_brain.nii.gz', anat_name));
            epi = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01.nii',subj_str{s}));
            out = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01_func2highres',subj_str{s}));
            reg_command = sprintf('epi_reg --epi=%s --t1=%s --t1brain=%s --out=%s', epi, t1, t1_brain, out);
            system(reg_command)
            
            
            fprintf('epi_reg completed for %s \n',subj_name{sn(s)})
            fprintf('Check the output of epi_reg using FSLeyes or some other visualization software.')
           
        end

     case 'FUNC:coreg_epi2epi'        % Coregister meanrun_01 to meanrun_01_func2struct
        % Need meanrun_01 in epi resolution coregistered to anatomical
        % example: bsp_imana('FUNC:coreg_meanepi_fsl',1,8)
        sn=varargin{1}; % subjNum
        runnum=varargin{2} %runNum
        
        subjs=length(sn);
        
        J = [];
        for s=1:subjs,
            
            % Get the directory of subjects anatomical and functional
            subj_anat_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            subj_func_dir = fullfile(base_dir, subj_str{s}, func_dir);

            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w', subj_str{s});

            
            t1 = fullfile(subj_dir, sprintf('%s_lpi.nii', anat_name));
            t1_brain = fullfile(subj_dir, sprintf('%s_lpi_brain.nii.gz', anat_name));
            epi = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01.nii',subj_str{s}));
            out = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01_func2highres',subj_str{s}));
            
            orig = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01.nii',subj_str{s}));
            orig_copied = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01_reg.nii',subj_str{s}));
            command = sprintf('cp %s %s',orig,orig_copied)
            system(command)

            ref_image = fullfile(subj_func_dir, sprintf('mean%s_ses-01_run-01_func2highres.nii.gz',subj_str{s}));
            command_gunzip = sprintf('gunzip %s', ref_image)
            system(command_gunzip)
            fprintf('gunzip completed for run %d \n',runs(r))
            [filedir,filestem, ext] = fileparts(ref_image);
            
            J.ref = {fullfile(filedir,filestem)};
            J.source = {orig_copied};
            J.other = {''};
            J.eoptions.cos_fun = 'nmi';
            J.eoptions.sep = [4 2];
            J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);
            fprintf('mean epi coregistered for %s \n',subj_name{sn(s)})
            
        end    
        
   case 'FUNC:make_samealign_fsl'        % Align functional images to rmeanepi of run 1, session 1
        % Aligns all functional images from both sessions
        % to rmeanepi of run 1 of session 1
        % example: bsp_imana('FUNC:make_samealign_fsl',1,8,[1:8])
        sn=varargin{1}; % subjNum
        runnum=varargin{2}; % runNum used for coregistration
        runs=varargin{3}; % runNum to coregister
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(baseDir,imagingDir,subj_name{sn(s)}));
            
            % Select image for reference
            % For ants-registered data: TSE 
            % P{1} = fullfile(fullfile(baseDir,anatomicalDir,subj_name{sn},'tse.nii'));
            % for tradition way: rmeanepi 
            P{1} = fullfile(fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('rmeanrun_%2.2d.nii',runnum)));
            
            % Select images to be realigned
            Q={};
            for r=1:numel(runs)
              Q{end+1}    = fullfile(baseDir,imagingDir,subj_name{sn(s)},sprintf('run_%2.2d.nii',runs(r)));
            end;
            
            % Run spmj_makesamealign_nifti
            spmj_makesamealign_nifti(char(P),char(Q));
            fprintf('functional images realigned for %s \n',subj_name{sn(s)})
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
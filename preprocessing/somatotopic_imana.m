function [ output_args ] = somatotopic_imana( what, varargin )

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
addpath(sprintf('%s/../matlab/imaging/freesurfer/',workdir));
addpath(sprintf('%s/../matlab/imaging/surfAnalysis/',workdir));
addpath(sprintf('%s/../matlab/imaging/coregtool/',workdir));

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

base_dir = sprintf('%s/Cerebellum/Somatotopic/raw',workdir);

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
subj_str = {'sub-02'};
subj_id_old = {'S2'};
subj_id  = [1];

ses_str = {'ses-01'};

% =========================================================================

switch what
    case 'ANAT:rename'    % Rename anatomical image to match script convention
        % Example usage: somatotopic_imana('ANAT:rename', 'sn', 1);
        
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            % Get the name of the anatpmical image
            anat_name = sprintf('%s_anat_mni_underlay_defaced.nii.gz', subj_id_old{s});

            
            source = fullfile(subj_dir, anat_name);
            dest   = fullfile(subj_dir, sprintf('%s_T1w.nii',subj_str{s}));

            [filepath,name,ext] = fileparts(source)
            if strcmp(ext,'.gz')
                gunzip(source)
                source = fullfile(filepath,name)
            end
            
            copyfile(source,dest);
            
        end % s (sn)
    case 'SUIT:isolate_segment'    % Segment cerebellum into grey and white matter
        % Example usage: somatotopic_imana('SUIT:isolate_segment', 'sn', 1);
        
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Isolate and segment the cerebellum for %s\n', subj_str{s})
            spm_jobman('initcfg')
            
            % Get the directory of subjects anatomical
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            % Get the name of the anatpmical image
            anat_name = sprintf('%s_T1w.nii', subj_str{s});

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
        % example usage: somatotopic_imana('SUIT:normalise_dartel')
        sn = subj_id; %subjNum
        vararginoptions(varargin, 'sn');
        
        for s = sn
            suit_subj_dir = fullfile(base_dir, subj_str{s}, 'suit');
            mkdir(suit_subj_dir)
            
            cd(suit_subj_dir)
            job.subjND.gray       = {fullfile(suit_subj_dir, sprintf('c_%s_T1w_seg1.nii', subj_str{s}))};
            job.subjND.white      = {fullfile(suit_subj_dir, sprintf('c_%s_T1w_seg2.nii', subj_str{s}))};
            job.subjND.isolation  = {fullfile(suit_subj_dir, sprintf('c_%s_T1w_pcereb_corr.nii', subj_str{s}))};
            suit_normalize_dartel(job);
        end % s (subjects)    

    case 'SUIT:save_dartel_def'    
        % Saves the dartel flow field as a deformation file. 
        % example usage: somatotopic_imana('SUIT:save_dartel_def')
        sn = subj_id; %subjNum
        vararginoptions(varargin, 'sn');

        for s = sn
            suit_subj_dir = fullfile(base_dir, subj_str{s}, 'suit', 'anat');

            cd(suit_subj_dir);
            anat_name = sprintf('%s_T1w', subj_str{s});
            suit_save_darteldef(anat_name);
        end % s (subjects)

    case 'SUIT:reslice'            % Reslice stuff into suit space 
        % run the case with 'anatomical' to check the suit normalization
        % Example usage: somatotopic_imana('SUIT:reslice','type','ResMS', 'mask', 'pcereb_corr')
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
                    files2map = sprintf('%s_T1w.nii', subj_str{s});
                    
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
            job.subj.affineTr  = {fullfile(suit_dir,sprintf('Affine_c_%s_T1w_seg1.mat', subj_str{s}))};
            job.subj.flowfield = {fullfile(suit_dir,sprintf('u_a_c_%s_T1w_seg1.nii', subj_str{s}))};
            job.subj.mask      = {fullfile(suit_dir, sprintf('c_%s_T1w_%s.nii', subj_str{s}, mask))};
            job.vox            = [2 2 2];
            suit_reslice_dartel(job);
            
            if ~strcmp(type,'anatomical')
                source=fullfile(subj_dir, '*wd*');
                dircheck(fullfile(out_dir));
                movefile(source,out_dir);
            end
            fprintf('- Resliced %s %s to SUIT\n', subj_str{s}, type);
        end % s (subject)

    case 'SURF:reconall'       % Freesurfer reconall routine
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        % Example usage: somatotopic_imana('SURF:reconall', 'sn', 1)
        
        sn   = subj_id; % subject list
        
        vararginoptions(varargin, {'sn'});
        % set freesurfer directory
        subj_fs_dir = fullfile(base_dir, fs_dir);
        
        for s = sn
            fprintf('- recon-all %s\n', subj_str{s});
            subj_dir = fullfile(base_dir, subj_str{s}, anat_dir);
            freesurfer_reconall(subj_fs_dir, subj_str{s}, ...
                                fullfile(subj_dir,sprintf('%s_anat_mni_underlay_defaced.nii.gz', subj_id_old{s})));
        end % s (sn)
    case 'SURF:fs2wb'          % Resampling subject from freesurfer fsaverage to fs_LR
        % Example usage: somatotopic_imana('SURF:fs2wb', 'sn', [1], 'res', 32)
        
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
        % Example usage: somatotopic_imana('SURF:run_all')
        
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        somatotopic_imana('SURF:reconall')
        somatotopic_imana('SURF:fs2wb', 'sn', sn);  
       
end

end


% Checking for directories and creating if not exist
function dircheck(dir)
if ~exist(dir,'dir')
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end
end
function [ output_args ] = hcptask_imana( what, varargin )

% %========================================================================================================================
% PATH DEFINITIONS

% Add dependencies to path
if isdir('/Volumes/diedrichsen_data$/data')
    workdir='/Volumes/diedrichsen_data$/data';
elseif isdir('/cifs/diedrichsen/data')
    workdir='/cifs/diedrichsen/data';
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end


%========================================================================================================================
global base_dir

base_dir = sprintf('%s/ExternalOpenData/HCP_UR100_tfMRI_full',workdir);

% defining the names of other directories
wb_dir   = 'surfaceWB';

D = dload(fullfile(base_dir,'participants.tsv'));

switch what
    case 'SUIT:isolate_segment'    % Segment cerebellum into grey and white matter
        % Example usage: nishimoto_bids_imana('SUIT:isolate_segment', 'sn', 1);
        
        sn = [1:50];
        vararginoptions(varargin, {'sn'});
               
        for s = sn
            subj = D.participant_id{s};
            fprintf('- Isolate and segment the cerebellum for %s\n', subj)
            
            % Get the name of the anatpmical image
            anat_name = fullfile(base_dir, subj, 'anat','T1w.nii');

            suit_subj_dir = fullfile(base_dir, subj, 'suit');
            dircheck(suit_subj_dir);
            
            n_anat_name   = fullfile(suit_subj_dir,'T1w.nii');            
            copyfile(anat_name,n_anat_name);
            
            % go to subject directory for suit and isolate segment
            suit_isolate_seg({n_anat_name}, 'keeptempfiles', 0);
        end % s (sn)
    case 'SUIT:normalise_dartel'   % SUIT normalization using dartel
        % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
        % example usage: nishimoto_imana('SUIT:normalise_dartel')
        sn = [1:50];
        vararginoptions(varargin, {'sn'});
               
        for s = sn
            subj = D.participant_id{s};
            suit_subj_dir = fullfile(base_dir, subj, 'suit');
            
            job.subjND.gray       = {fullfile(suit_subj_dir, 'c_T1w_seg1.nii')};
            job.subjND.white      = {fullfile(suit_subj_dir, 'c_T1w_seg2.nii')};
            job.subjND.isolation  = {fullfile(suit_subj_dir, 'c_T1w_pcereb_corr.nii')};
            suit_normalize_dartel(job);
        end % s (subjects)    
    case 'SUIT:save_dartel_def'    
        % Saves the dartel flow field as a deformation file. 
        % example usage: nishimoto_imana('SUIT:save_dartel_def')
        sn = [1:50];
        vararginoptions(varargin, {'sn'});
               
        for s = sn
            subj = D.participant_id{s};
            suit_subj_dir = fullfile(base_dir, subj, 'suit');

            cd(suit_subj_dir);
            anat_name = sprintf('', subj_str{s});
            suit_save_darteldef(anat_name);
        end % s (subjects)
    case 'SUIT:mask_cereb'         % Make cerebellar mask using SUIT
        % Example usage: nishimoto_imana('SUIT:mask_cereb', 'glm', 1, 'ses', 1)
        
        sn       = subj_id; % list of subjects
        glm      = 1;           % glm number
        ses = 1;
        
        vararginoptions(varargin, {'sn', 'glm', 'ses'})

        
        for s = sn
            suit_dir = fullfile(base_dir, subj_str{s}, 'suit', 'anat');
            glm_dir = fullfile(base_dir, subj_str{s}, 'estimates', sprintf('glm%02d', glm), sprintf('ses-%02d', ses));
            
            mask  = fullfile(glm_dir, 'mask.nii'); % mask for functional image
            
            suit  = fullfile(suit_dir, 'cereb_prob_corr_grey');
%             suit  = fullfile(suit_dir, sprintf('c1%s_T1w_lpi.nii', subj_str{s})); % cerebellar mask grey (corrected)
            suit_glm_dir = fullfile(base_dir, subj_str{s}, 'suit', sprintf('glm%02d', glm), sprintf('ses-%02d', ses)); dircheck(suit_glm_dir);
            omask = fullfile(suit_glm_dir, 'maskbrainSUITGrey2.nii'); % output mask image - grey matter
            
            cd(suit_dir);
            spm_imcalc({mask,suit}, omask, 'i1>0 & i2>0.', {});
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
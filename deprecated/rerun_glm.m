%% Use designmatrix.npy to rerun GLM
% % Import NumPy from Python
% np = py.importlib.import_module('numpy');
% 
% % Load the .npy file
% X = np.load('/Volumes/diedrichsen_data$/data/FunctionalFusion/MDTB/derivatives/sub-02/estimates/ses-s1/sub-02_ses-s1_designmatrix.npy');
% 
% % Convert the Python array to MATLAB array
% X = double(X); % SPM.xX.xKXs.X

%% Use SPM.mat to rerun GLM
% Load the SPM.mat file
source_dir = [workdir '/Cerebellum/super_cerebellum/sc1/']
design_file = [source_dir 'GLM_firstlevel_7/s02/SPM.mat'];
load(design_file);

% Import the data
% Loop through run 01 to 08 for both s1 and s2 and concatenate the timeseries
for run = 1:8
    for ses = 1:2
        data_file = [source_dir '/imaging_data_fix/sub-02_ses-s' num2str(ses) '_run-0' num2str(run) '_fix.nii'];
        data = niftiread(data_file);
        if run == 1 && ses == 1
            y_raw = data;
        else
            y_raw = cat(4, y_raw, data);
        end
    end
end
y_raw = double(y_raw);


X = SPM.xX.xKXs.X;
y_filt = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
B = SPM.xX.pKX*y_filt;               %-Parameter estimates
y_res  = spm_sp(‘r’,SPM.xX.xKXs,y_filt);       %-Residuals
y_hat = X(:,reg_interest)*B(reg_interest,:); %- predicted values


y_adj = y_hat + y_res;




save('new.mat','-struct','S','-v7.3','-nocompression')
save('test.mat','-struct','SPM')

/cifs/diedrichsen/data/Cerebellum/super_cerebellum/sc1/GLM_firstlevel_7/s31/
matlab -nodesktop -nosplash -r "load('SPM.mat'); save('SPM.mat', '-struct',  'SPM', '-v7'); exit"

for i in /cifs/diedrichsen/data/Cerebellum/super_cerebellum/sc1/GLM_firstlevel_7/s*/; do echo $i;
cd $i;
matlab -nodesktop -nosplash -r "load('SPM.mat'); save('SPM.mat', '-struct',  'SPM', '-v7' ); exit"
 done

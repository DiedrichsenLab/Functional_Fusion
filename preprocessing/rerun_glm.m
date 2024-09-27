% Import NumPy from Python
np = py.importlib.import_module('numpy');

% Load the .npy file
X = np.load('/Volumes/diedrichsen_data$/data/FunctionalFusion/MDTB/derivatives/sub-02/estimates/ses-s1/sub-02_ses-s1_designmatrix.npy');

% Convert the Python array to MATLAB array
X = double(X);




if (~isfield(SPM.xX,‘pKX’)) % SPM not run -
  y_filt = spm_filter(SPM.xX.K,y_raw);
  SPM.xX.xKXs.X = spm_filter(SPM.xX.K,SPM.xX.X);
  SPM.xX.pKX = pinv(SPM.xX.xKXs.X);
  B = SPM.xX.pKX*y_filt;               %-Parameter estimates
  y_res  = y_filt - SPM.xX.xKXs.X*B;       %-Residuals
  y_hat = SPM.xX.xKXs.X(:,reg_interest)*B(reg_interest,:); %- predicted values
else
  y_filt = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
  B = SPM.xX.pKX*y_filt;               %-Parameter estimates
  y_res  = spm_sp(‘r’,SPM.xX.xKXs,y_filt);       %-Residuals
  y_hat = SPM.xX.xKXs.X(:,reg_interest)*B(reg_interest,:); %- predicted values
end;
y_adj = y_hat + y_res;
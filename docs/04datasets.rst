Datasets and Data Import
========================

Each dataset is stored in a separate directory in the `basedir` directory. Depending to the dataset, the preprocessed time series or individual effect-size estimates are stored.
The Data Set class `DataSet` is designed to be the entry of getting the data in standard format. To be able to reuse a lot of the code across data sets, it is useful if the way the data is
stored is homogenous across datasets

Directory structure
-------------------
The folder structure of derivatives (example for DataSetNative)

:: 

    basedir/
    |--- participant.tsv: # TSV file with subject information/
    |--- README.md:       # General information file for data/
    |--- derivatives/
    | |--- ffimport
    | | |--- sub-<id>/
    | | | |---anat/
    | | | | |---- sub-<id>_T1w.nii                             # Native space T1w (space defining)
    | | | | |---- sub-<id>_label-CSF_probseg.nii               # probabilistic segmentation (CSF)
    | | | | |---- sub-<id>_label-GM_probseg.nii                # probabilistic segmentation (GM)
    | | | | |---- sub-<id>_label-WM_probseg.nii                # probabilistic segmentation (WM)
    | | | | |---- sub-<id>_space-32k_hemi-L_white.surf.gii     # 32K white matter surface
    | | | | |---- sub-<id>_space-32k_hemi-L_pial.surf.gii      # 32K pial surfaceces
    | | | | |---- sub-<id>_desc-brain_mask.nii                 # Mask of within brain tissue
    | | | | |---- sub-<id>_label-GMc_probseg.nii                # probabilistic segmentation (GM-cereb)
    | | | | |---- sub-<id>_label-WMc_probseg.nii                # probabilistic segmentation (WM-cereb)
    | | | | |---- sub-<id>_label-GMb_probseg.nii                # probabilistic segmentation (GM-rest)
    | | | | |---- sub-<id>_label-WMb_probseg.nii                # probabilistic segmentation (WM-rest)
    | | | | |---- sub-<id>_desc-cereb_mask.nii                 # hand corrected cerebellar mask in functional space
    | | | | |---- sub-<id>_space-SUIT_xfm.nii                  # coordinate transformation file into SUIT
    | | | | |---- sub-<id>_space-MNI_xfm.nii                   # coordinate transformation file into MNI152Nonlin
    | | | |---func/ 
    | | | | |---- ses-<id>
    | | | | | |--- sub-<id>_ses-<id>_designmatrix.npy   # Design matrix used for estimation (optional)
    | | | | | |--- sub-<id>_ses-<id>_mask.nii           # Brain mask in functional space
    | | | | | |--- sub-<id>_ses-<id>_reginfo.tsv        # Information on regression estimate values gressor in design matrix
    | | | | | |--- sub-<id>_ses-<id>_resms.nii          # Model Variance (ResMS.nii in SPM, sigmasquareds.nii.gz in FSL)
    | | | | | |--- sub-<id>_ses-<id>_run-<id>_reg-<id>_beta.nii    # Parameter estimates (beta_0001.nii in SPM, pe1.nii.gz in FSL)
    | | | | | |--- sub-<id>_ses-<id>_run-<id>.nii       # Preprocessed time series file (4D, optional)
    | | | | | |--- ....                                   
    | | | | |---- ses-<id>
    | | | | | |--- ....                                 # Second session 
    | |--- ffextract
    | | |--- group/                                         # Group averaged data 
    | | | |--- group_ses-<id>_info-<type>.tsv               # Information file for specific session / data type
    | | | |--- group_space-<id>_ses-<id>_<type>.dscalar.nii # Extracted data for specific space / session / data type
    | | | |--- ...
    | | |--- sub-<id>/
    | | | |--- sub_<id>_ses-<id>_info-<type>.tsv               # Information file for specific session / data type
    | | | |--- sub_<id>_space-<id>_ses-<id>_<type>.dscalar.nii # Extracted data for specific space / session / data type
    | | | |--- ...



participant.tsv file
-----------------------

`participant.tsv` is a tab-delimted text file that contains the subject information with the following recommended columns:

* participant_id: Subject ID with prepended sub- (e.g. sub-01)
* sex: M or F 
* age: Age of the subject in years since birth

Optional columns include: 

* handedness: L or R
* group: Group of the subject 
* exclude: 0 or 1 (1 = exclude from analysis)


reginfo.tsv file
-------------------

`reginfo.tsv`` is a tab-delimted text file that contains the information for a specific functional session for the regression estimates. By default these regression estimates are save 
under the name `sub-<id>_run-<id>_reg-<id>_beta.nii`, where numerical ids are printed with 2 places (i.e 01, 02, 03). 

The file should contain the following columns:

* run: Run number (numerical): Run number 
* reg_id: Regressor ID (numerical) 

Optional columns can include: 

* 
* 


Import Anatomical and MNI normalization parameters from SPM (Segement)
----------------------------------------------------------------------

If you run the SPM Segmentation algorithm in a source directory, the anatomical, segmentations, and normalization parameters to MNI152Nonlin can be imported by:

.. code-block:: python

    import import_data as id
    source_dir = <Directory where you ran segementation (outside functional fusion) >
    dest_dir = '<base_dir/derivates/sub-xx/anat'
    anat_name = '<something>.nii'
    id.import_anat(source_dir,dest_dir,anat_name,'sub-xx')


Import Cortical surfaces from Freesurfer reconstruction
--------------------------------------------------------

Import SUIT normalization
-------------------------
Run SUIT isolation, and normalization outside of the Functional Fusion framework.
To produce the cerebellar mask in functional space, you need to combine the functional mask from the GLM (mask.nii), the cerebellar mask from suit (c_anatimical_pcerebe(_corr).nii) and the gray matter segmentation (c_anatomical_seg1.nii)

.. code-block:: matlab

    mask  = fullfile(glm_dir, 'mask.nii'); % mask for functional image
    suitm  = fullfile(suit_dir, 'c_anatomical_pcereb_corr.nii');
    gray  = fullfile(suit_dir, c_anatomical_seg1.nii)); %
    omask = fullfile(suit_glm_dir, 'maskbrainSUITGrey.nii'); %
    spm_imcalc({mask,suitm,gray}, omask, 'i1>0 & i2>0 & i3>0.01', {});


Additionally, you need to save the non-linear transformation between SUIT and individual subject space as a deformation file.

.. code-block:: matlab

    suit_save_darteldef(<c_anat_name>,'wdir',workingdirectory)

Then you can run ,,import_suit`` in Python to copy and rename.

Import functional estimates and design matrix from SPM
------------------------------------------------------

Import task-specific beta files (ex: beta_0001.nii) for each subject, and rename them according to subject, session, run, and condition/ regressor (ex: sub-01_ses-01_run-01_reg-00_beta.nii). 

Import the SPM_info.tsv file for each subject and rename according to subject and session (ex: sub-01_ses-01_reginfo.tsv).

Save the prewhitened design matrix (SPM.xX.nKX) as a numpy array (ex: sub-01_ses-01_designmatrix.npy). 
To do this, run this sequence of code in Matlab:

.. code-block:: matlab

    load('SPM.mat')
    nKX = SPM.xX.nKX;
    save('/directory_of_your_choice/nKX_data.mat','nKX')

and this sequence of code in Python: 

.. code-block:: matlab
    
        import numpy as np
        import scipy.io as sio
        nKX_data = sio.loadmat('/directory_of_your_choice/nKX_data.mat')
        np.save('/directory_of_your_choice/nKX.npy',nKX_data)

---------------

Add the information to dataset_description.
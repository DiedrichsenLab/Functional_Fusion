# Functional_Fusion
Diedrichsen Lab, Western University

This repository is the structure and preprocessing code for the multi-data set project in the Diedrichsenlab.

## Installation and dependencies
This project depends on several third party libraries, including:

[numpy](https://numpy.org/) (version>=1.22.2)

nibabel []

[nilearn](https://nilearn.github.io/stable/index.html) (version>=0.9.0), ...

	pip install numpy nilearn ...

[nitools]
    pip install neuroimagingtools

Or you can install the package manually from the original binary source as above links.

Once you clone the functional fusion repository, you need to add it to your PYTHONPATH, so you can import the functionality. Add these lines to your .bash_profile, .bash_rc .zsh_profile file... 

```
PYTHONPATH=/Users/.../parentdir:${PYTHONPATH}
export PYTHONPATH
```


## Structures of the project
### Overall structure
![ScreenShot](docs/data_flow.png)



#### 1. Data Set class
As can be seen in above diagram, the integration across data sets is achieved through  `DataSet` objects, with each data set being an instantiation of that type. The `DataSet` has access to the subject list and the individual preprocessed imaging data. The main function of the Data set class is the  `get_data()` function, which provides the
processed data `Y` in desired format `(N * P)` where `N` is the number of measurements (tasks, etc) `P` is the number of brain locations for a specific subject.

### 2. Atlas map class

To read out different data sets in a consistent anatomical location, we need to
decide an `Atlas` to be used. For each subejct, we need to have an `AtlasMap`, which provides the mapping function from the raw data space (per subject) to the common atlas space. We will have subclasses `AtlasMapMNI`, `AtlasMapFS32K`, and `AtlasMapSUIT3`. For each subject (and hence Dataset), there will be a seperate Instantiation (Object) of this class. If two dataset share the same raw data space for the same subject, they can rely on the same AtlasMap.

### 3. Altas
Each group atlas also has some subject / study independent behaviors that will be implemented in the `Atlas` Class. Atlas are paired with a specific `Atlas Map` and have the information that it takes to map the `P` locations back into brain space. Subclasses indicates the brain areas to be studied with a structural template (i.e MNI 152).

## Directory structure for derivatives
=======
The Data Set class `DataSet` is designed to be the entry of getting the data in standard format. To be able to reuse a lot of the code across data sets, it is useful if the

### Derivatives structure

The folder structure of derivatives

    derivatives/
        │   README.md
        │
        └───group/
        │
        │       ...
        │
        └───sub-<label>/
        │       └───anat/
        │       │       sub-<id>_T1w.nii                # Native space T1w (space defining)
        │       │       sub-<id>_label-CSF_probseg.nii               # probabilistic segmentation (CSF)
        │       │       sub-<id>_label-GM_probseg.nii                # probabilistic segmentation (GM)
        │       │       sub-<id>_label-WM_probseg.nii                # probabilistic segmentation (WM)
        │       │       sub-<id>_space-32k_hemi-L_white.surf.gii     # 32K white matter surface
        │       │       sub-<id>_space-32k_hemi-L_pial.surf.gii      # 32K pial surfaceces
        |       |       sub-<id>_desc-brain_mask.nii                 # Mask of within brain tissue
        │       └───suit/
        │       │       sub-<id>_label-GMc_probseg.nii                # probabilistic segmentation (GM-cereb)
        │       │       sub-<id>_label-WMc_probseg.nii                # probabilistic segmentation (WM-cereb)
        │       │       sub-<id>_label-GMb_probseg.nii                # probabilistic segmentation (GM-rest)
        │       │       sub-<id>_label-WMb_probseg.nii                # probabilistic segmentation (WM-rest)
        │       │       sub-<id>_desc-cereb_mask.nii                  # hand corrected cerebellar mask in functional space
        |       | 		sub-<id>_space-SUIT_xfm.nii 				  #	coordinate transformation file into native
        │       └───func/
          								sess-s1/
        |				| 					Minimally preprocessed fMRI data, ideally in the subjects original space
        │       │       		sub-<label>_ses-<label>_run-<label>_bold.nii[.gz]
        |				|
        |				|						Information of different characteristics of runs (phase-encoding direction, etc)
        |				|						should be stored in a separate json or tsv file....
        │       │
        │       └───estimates/
          								sess-s1/
    			    │               beta_info.tsv: Information on regression estimate values structure
    			    									TSV-file with obligatory columns
    			    										run: run number (reflected in file name)
    			    										reg_id: regressor id (reflected in file name)
    			    										reg_num: column number of regressor in design matrix
    			    								sub-<label>_ses-<label>_matrix.npy: Design matrix used for estimation
        │                     sub-<label>_ses-<label>_run-<label>_reg-<label>_beta.nii
        │                     sub-<label>_ses-<label>_run-<label>_reg-<label>_beta.nii
        │                     sub-<label>_ses-<label>_mask.nii
        │                     sub-<label>_ses-<label>_resms.nii


### AtlasMap structure

Need to be discussed later.


## Import data to the Functional Fusion framework
### Import Anatomical and MNI normalization parameters from SPM (Segement)
If you run the SPM Segmentation algorithm in a source directory, the anatomical, segmentations, and normalization parameters to MNI152Nonlin can be imported by: 
```
    import import_data as id
    source_dir = <Directory where you ran segementation (outside functional fusion) >
    dest_dir = '<base_dir/derivates/sub-xx/anat'
    anat_name = '<something>.nii' 
    id.import_anat(source_dir,dest_dir,anat_name,'sub-xx') 
```
### Import Cortical surfaces from Freesurfer reconstruction 
### Import SUIT normalization
Run SUIT isolation, and normalization outside of the Funtional Fusion framework. Additionally, you need to save the non-linear transformation between SUIT and individual subject space as a deformation file. 

```
    suit_save_darteldef(<c_anat_name>,'wdir',workingdirectory)
```
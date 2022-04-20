Functional_Fusion
====
Diedrichsen Lab, Western University

This repository is the structure and preprocessing code for the multi-data set project in the Diedrichsenlab.

Installation and dependencies
------
This project depends on several third party libraries, including:

[numpy](https://numpy.org/) (version>=1.22.2)

nibabel []

[nilearn](https://nilearn.github.io/stable/index.html) (version>=0.9.0), ...

	pip install numpy nilearn ...

Or you can install the package manually from the original binary source as above links.

Structures of the project
------
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

    @Ana: IMHO, all different individual derivatives should be under a ses-<label> dir, under sub<label> dir.
          I added this suggestion to the diagram.
    @Joern: This is clearly useful for some experiments, but not for others. Certainly, we want only one anatomical directory per study - as we want the data in a common space. I have removed the session for now - this is meant to be the smallest common denominator across all studies. Extensions can always be made for longer, multi-session studies.

    derivatives/
        │   README.md
        │
        └───group/
        │
        │       ...
        │
        └───sub-<label>/
        │       └───anat/
        │       │       sub-<label>_desc-preproc_T1w.nii                # Native space T1w (space defining)
        │       │       sub-<label>_label-CSF_probseg.nii               # probabilistic segmentation (CSF)
        │       │       sub-<label>_label-GM_probseg.nii                # probabilistic segmentation (GM)
        │       │       sub-<label>_label-WM_probseg.nii                # probabilistic segmentation (WM)
        │       │       sub-<label>_space-32K_hemi-L_white.surf.gii     # 32K white matter surface
        │       │       sub-<label>_space-32K_hemi-L_pial.surf.gii      # 32K pial surfaceces
        |       |       sub-<label>_desc-brain_mask.nii                 # Mask of within brain tissue
        │       └───suit/
        │       │       sub-<label>_label-GMc_probseg.nii                # probabilistic segmentation (GM-cereb)
        │       │       sub-<label>_label-WMc_probseg.nii                # probabilistic segmentation (WM-cereb)
        │       │       sub-<label>_label-GMb_probseg.nii                # probabilistic segmentation (GM-rest)
        │       │       sub-<label>_label-WMb_probseg.nii                # probabilistic segmentation (WM-rest)
        │       │       sub-<label>_desc-cereb_mask.nii                  # hand corrected cerebellar mask
        │       │       sub-<label>_label-WMb_probseg.nii                # probabilistic segmentation (WM-rest)
        |       |       
        │       └───func/
        │       │       file naming (@Ana)
        │       │
        │       │       Template suggested by BIDS:
        │       │       sub-<label>[_ses-<label>]_task-<label>[_acq-<label>][_ce-<label>][_dir-<label>]
        | [_rec-<label>] \
        │       │           [_run-<index>][_echo-<index>]_<contrast_label>.nii[.gz]
        │       │
        │       │       [_acq-<label>][_ce-<label>][_dir-<label>][_rec-<label>] are optional keys/values.
        │       │       Multi-echo data MUST be split into one file per echo. We can skip this if not Multi-Echo.
        │       │
        │       │       Check:
        │       │       https://bids-specification.readthedocs.io/en/v1.2.0/04-modality-specific-files/ \
        │       │           01-magnetic-resonance-imaging-data.html
        │       │
        │       │       Example for raw data after conversion from Dicom to NIfTI, considering a task named
        │       │       Theory-of-Mind (TOM):
        │       │       sub-01_task-TOM_run-01_bold.nii.gz
        │       │
        │       │       For preprocessed volume data, we should probably take advantage of the SPM prefixes.
        │       │       We might be using data that need to be tagged with the following:
        │       │       a - slice timing correction
        │       │       r - resliced (this can be from coregistration or realignment)
        │       │       u - undistorted, (from Realign unwarp - which requires reslicing)
        │       │       w - warped - typically this is done by normalization
        │       │
        │       │       In addition, AFAIK, there's no BIDS convention for different system coordinates. So,
        │       │       as suggested for the surface, we can simply add another key/value for that.
        │       │
        │       │       Example:
        │       │       wurasub-01_ses-01_task-TOM_space-MNI152_dir-ap_run-01_bold.nii.gz
        │       │
        │       │       Motion files should also be included here.
        │       │       Example:
        │       │       rpsub-01_ses-01_task-TOM_dir-ap_run-01_bold.txt
        │       │
        │       │       Paradigm-descriptors files for all runs to build the design matrix should also be
        │       │       included here.
        │       │       Example:
        │       │       sub-01_ses-01_task-TOM_dir-ap_run-01_events.tsv
        │       │
        │       └───contrast/
        │               @Ana: This dir should be renamed as 'first-level_analysis' because we may not only want
        │                     to include contrast (stat) maps, but also effect size or effect variance maps.
        │                     Inside this dir, we should have sub-folders dedicated to each run and one for ffx.
        │                     Examples of sub-folders names:
        │                     results_volume_task-<label>_run<label>
        │                     results_volume_task-<label>_ffx
        │                     results_fsaverage7_task-<label>_run<label>
        │                     results_fsaverage7_task-<label>_ffx
        │
        │               beta_info.tsv file structure (@Ana / @Ladan)
        │       ...
        │
        └───sub002/
        │       ...
        │
        └───subxxx/
                ...

### AtlasMap structure

Need to be discussed later.

# Data Check
Diedrichsen Lab, Western University

Guide on how to check data quality.

## Surface checks
Surfaces are located in folder anat:

derivatives/
        └───sub-<label>/
        │       └───anat/
        │       │       sub-<id>_T1w.nii                # Native space T1w (space defining)
        │       │       sub-<id>_label-CSF_probseg.nii               # probabilistic segmentation (CSF)
        │       │       sub-<id>_label-GM_probseg.nii                # probabilistic segmentation (GM)
        │       │       sub-<id>_label-WM_probseg.nii                # probabilistic segmentation (WM)
        │       │       sub-<id>_space-32k_hemi-L_white.surf.gii     # 32K white matter surface
        │       │       sub-<id>_space-32k_hemi-L_pial.surf.gii      # 32K pial surfaceces
        |       |       sub-<id>_desc-brain_mask.nii                 # Mask of within brain tissue

To check that surfaces were reconstructed correctly and match the T1, pull up workbench:

```
wb_view
```

### Pial Surface

Then open the native space T1 (sub-01_T1w.nii) and the pial matter surfaces for the left and right hemisphere (sub-01_space-32k_hemi-L_white.surf.gii & sub-01_space-32k_hemi-R_pial.surf.gii ). In the montage tab, click on ```All``` (is set to ```Montage``` by default).

You should be seeing the T1 image and reconstructed pial surface now. Rotate the image to see whether the surface is aligned with the sulci in the T1. An extreme mismatch would look like this:

![PialSurface_largeOffset](docs/surface_check_1.png)

A subtler mismatch could look like this:

![PialSurface_slightOffset](docs/surface_check_2.png)

You can see a perfect match below.

![PialSurface_Match](docs/surface_check_3.png)

### White Matter Surface

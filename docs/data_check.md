# Data Check
Diedrichsen Lab, Western University

Guide on how to check data quality.

CRN, 2022

## Surfaces
Surfaces are located in folder anat:

    derivatives/
            └───sub-<label>/
            │       └───anat/
            │       │       sub-<id>_T1w.nii                # Native space T1w (space defining)
            │       │       sub-<id>_space-32k_hemi-L_white.surf.gii     # 32K white matter surface
            │       │       sub-<id>_space-32k_hemi-L_pial.surf.gii      # 32K pial surfaceces



To check that surfaces were reconstructed correctly and match the T1, pull up workbench:

```
wb_view
```

### Pial Surface

Open the native space T1 (sub-01_T1w.nii) and the pial surfaces for the left and right hemisphere (sub-01_space-32k_hemi-L_white.surf.gii & sub-01_space-32k_hemi-R_pial.surf.gii ). In the montage tab, click on ```All``` (is set to ```Montage``` by default). You might have to also click ```On``` for the T1 image in the overlay toolbox to show the T1.

You should now see the T1 image and reconstructed pial surface. Rotate the image to see whether the surface is aligned with the sulci in the T1. An extreme mismatch would look like this:

<img src="../docs/surface_check_1.png" alt="PialSurface_largeOffset" width="400"/>

A subtler mismatch of the right hemisphere could look like this:

<img src="../docs/surface_check_2.png" alt="PialSurface_slightOffset" width="400"/>

And here is a perfect match between pial surface and T1 image:

<img src="../docs/surface_check_3.png" alt="PialSurface_Match" width="400"/>

### White Matter Surface

Remove the pial surface and add the white matter surface reconstruction. Make sure you removed the pial surface from you viewer, otherwise you won't see the white matter surface under the pial surface. Follow the same process as for the pial surface check, click on ```All``` and rotate the image to see how T1 and the white matter surface align. Mismatches for the white matter surface should be easy to spot since the white matter surface at the cortex should neatly fit into the grey matter on the T1 like below.

<img src="../docs/surface_check_4.png" alt="WMSurface_Match" width="600"/>


Focus on the white matter surface going into the coronal T1 slice, not the saggital slice. Since the white matter connects the hemispheres only at the corpus callosum, which is hidden behind the rest of the white matter in this view, it's hard to spot mismatches to the saggital T1 slice. Instead, look at where the white matter surface enters the coronal slice and the axial slice.



## Functional Data

### Group Averaged Functional Data

To check whether your functional data makes sense, extract the group average (use function ```group_average()```) and open the group averaged data in the workbench viewer (```wb_view```). Depending on which space you the data is extracted in, you will first have to load the space defining file. For example, if you want to insepct data on the fs32k surfaces, you will first have to load the .spec file for the fs32k space. Afterwards, you can load in the group average cifti file.

Group averaged functional data is located in folder data:

    derivatives/
            └───group/
            │       group_ses-<label>_space-32k_<type>.dscalar.nii      # Group averaged functional estimates in 32K 

Look at some contrasts you are familiar with and where you know what to expect. For example, a left hand movement task should give you high activation in the right M1 hand area:

<img src="../docs/left_hand.png" alt="LeftHandActivation" width="300"/>

And the right hand movement task should give you left M1 hand area activation.

<img src="../docs/right_hand.png" alt="RightHandActivation" width="300"/>

You can also check visual tasks, where you would expect activation in the visual cortices:

<img src="../docs/visual_task.png" alt="VisualTaskActivation" width="300"/>

These inspections should confirm that on average, you are seeing activity in the areas that you would expect for the different tasks.

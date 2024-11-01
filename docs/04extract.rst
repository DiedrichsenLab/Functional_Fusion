Data Extraction
###############

Extraction reads out the effect-size estimates or timeseries data (the type) from a dataset out in a specific atlas space (the atlas). The result will be stored in a cifti-file in the subjects data directory, with a uniform Spatial and Data dimension.  This process is done can be done for all subject using the function :py:meth:`dataset.DataSet.extract_all`.

For example, extracting sessions 1 from the MDTB dataset, with a split half-estimate for the condition in MNISymC3,  would look like this: 

.. code-block:: python

    dataset = ds.DataSetMDTB(base_dir + '/MDTB')
    dataset.extract_all(ses_id='ses-s1',
                        type='CondAll',
                        atlas='MNISymC3',
                        smooth=2)

The resulting data for each subject and session is stored in a cifti-file in the ``basedir/derivatives/<subj_id>/data`` directory under the name ``sub-<xx>_space-<atlas>_ses-xx_<type>.dscalar.nii``. The description of the data-axis in the cifti-file is stored in ``sub-xx_ses-xx_<type>.tsv`` (note that there is of course only one of these files for all atlas spaces). 

That's all you need. The rest of documentation will explain the different steps in the extraction process in detail:

Spatial resampling
------------------
Along the spatial dimension, extraction reads out the data from the source space and brings them into the atlas space (``atlas.space``). Not that each atlas has a ``space`` and within this space a set of defined xyz coordinates ``atlas.world``.This mapping between different spaces is store in an ``AtlasMap`` object, and can rely either on a subject-specific deformation, a atlas-to-atlas deformation, or even a combination of both.

The specific mapping rules for the dataset are defined in the dataset method :py:meth:`dataset.DataSet.get_atlasmaps`. We have implemented a set of default mapping rules in the generic Dataset class, and some more specific ones in ``DataSetNative`` and ``DataSetMNIVol``. If you require more specialize mapping rules, you can overwrite this behavior by defining a ``get_atlasmaps`` function for your specific dataset.Each atlas map also uses a mask image, which defines the source voxel-space, and the voxels that are available.

* Defined in ``DataSet``

    * ``SUIT``: Here we use the individual deformation map into suit space, and the cerebellar mask in source space.
         .. code-block:: python

            deform = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'

    * ``MNI152NLin2009cSymC`` & ``MNI152NLin6AsymC``: Here we use the individual deformation map into SUIT space and then from SUIT space into the MNI-space. Mask is the cerebellar mask in source space. We used this type of mapping for our first papers (Nettekoven et al., 2024), but will replace it with a direct mapping from native to MNI space (see below).

         .. code-block:: python

            deform1 = am.get_deform(atlas.space, source='SUIT')
            deform2 = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'

    * ``MNI152NLin2009cSym``, ``MNI152NLin2009cAsym``: Here we use ideally direct mapping between native space and the corresponding space, found in the individual anat folder. If the corresponding deformation file is not found, a warning is given, We use then the individual deformation into MNI152NLin6Asym space and the functional mask image for that subject and session. Then we use the deformation from MNI152NLin6Asym to MNI152NLin2009cSym or MNI152NLin2009cAsym.

            .. code-block:: python

                deform = self.anatomical_dir.format(sub) + f'/{sub}_space-MNI_xfm.nii'
                mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'

    * ``MNI152NLin6Asym``: We use the individual deformation into MNI152NLin6Asym (SPM segmentation) space and the functional mask image for that subject and session.

            .. code-block:: python

                deform = self.anatomical_dir.format(sub) + f'/{sub}_space-MNI_xfm.nii'
                mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'

    * ``fs32k``: The atlasmap is defined by the individual pial and white surfaces. The surface coordinates need to be in the source space - so if the source data is in native space, the surfaces need to be defined in native space. If the source data is in MNI152 space, the individual surfaces need to be in MNI152 space. The mask is the functional mask image for that subject and session.

         .. code-block:: python

            pial = self.anatomical_dir.format(sub) + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
            white = self.anatomical_dir.format(sub) + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
            mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'

* Defined in ``DataSetNative``

    Same as in ``DataSet`` 

* Defined in ``DataSetMNIVol``

    * ``MNI152NLin6Asym,MNI152NLin2009cSym,MNI152NLin2009cAsym``: For any deformation into an MNI space, we either use no deformation (if the atlas.space and dataset.space match), or we use the deformation between the two MNI spaces. No individual deformation is expected. 

            .. code-block:: python

                deform = None
                mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'

Datafile specification
----------------------
Depending on the type and dataset, the filenames of the raw datafiles need to be correctly specified. This is done in the method  :py:meth:`dataset.DataSet.get_data_fnames`.

The default behavior is:

* For ``type == 'TSeries'``: ``derivaties/estimates/{participant_id}_{session_id}_run-01.nii'``
* For ``type == 'task / cond'``: ``derivaties/estimates/{participant_id}_{session_id}_reg_00_beta.nii'``

If the naming convention differs, your Dataset class needs to overwrite this function.


Data aggregation
----------------
After the data has been sampled into atlas space, it is (potentially) aggregated across different runs and conditions. This dataset-specific function is done in the function :py:meth:`dataset.DataSet.condense_data`.

Typically, there are different `type`s:

* ``'TSeries'``: No aggregation (or z-standardization).
* ``'CondAll'``: A single estimate per condition, averaged across all runs.
* ``'CondHalf'``: Two estimates per condition, one per half
* ``'CondRun'``: A separate estimate per condition and run.

The averaging is done in the function :py:meth:`dataset.optimal_contrast`, which can take into account the first-level design matrix. This procedure will result in the same estimate that you would have gotten if you had defined a design matrix with a regressor for each condition across runs.

After thius step (and depending how baseline has been modeled) the function :py:meth:`dataset.remove_baseline` can be called to remove the mean of the voxels across conditions within each run. This is totally optional.  

Finally, we are dividing the beta estimates by the estimate of the noise standard-deviation per voxel, using :py:meth:`dataset.prewhiten`, coming from the resms.nii.


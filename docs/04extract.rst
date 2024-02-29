Data Extraction
###############

Extraction is the process by which the effect-size estimates or timeseries data from a dataset are read out. The result will be stored in a cifti-file, with a uniform Spatial and Data dimension.

Spatial resampling
------------------
Along the spatial dimension, extraction reads out the data from the source space and brings them into the atlas space (``atlas.space``). Not that each atlas has a ``space`` and within this space a set of defined xyz coordinates ``atlas.world``.This mapping between different spaces is store in an ``AtlasMap`` object, and can rely either on a subject-specific deformation, a atlas-to-atlas deformation, or even a combination of both.

The specific mapping rules for the dataset are defined in the dataset function ``get_atlasmap``. We have implemented a set of default mapping rules in the generic Dataset class, and some more specific ones in ``DataSetNative`` and ``DataSetMNIVol``. If you require more specialize mapping rules, you can overwrite this behavior by defining a ``get_atlasmaps`` function for your specific dataset.Each atlas map also uses a mask image, which defines the source voxel-space, and the voxels that are available.

* Defined in ``DataSet``

    * ``SUIT``: Here we use the individual deformation map into suit space, and the cerebellar mask in source space.
         .. code-block:: python

            deform = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'

    * ``MNI152NLin2009cSymC`` & ``MNI152NLin6AsymC``: Here we use the individual deformation map into SUIT space and then from SUIT space into the MNI-space. Mask is the cerebellar mask in source space.

         .. code-block:: python

            deform1, _ = am.get_deform(self.atlas_dir, atlas.name, source='SUIT2')
            deform2 = self.suit_dir.format(sub) + f'/{sub}_space-SUIT_xfm.nii'
            mask = self.suit_dir.format(sub) + f'/{sub}_desc-cereb_mask.nii'

    * ``fs33k``: The atlasmap is defined by the individual pial and white surfaces. The surface coordinates need to be in the source space - so if the source data is in native space, the surfaces need to be defined in native space. If the source data is in MNI152 space, the individual surfaces need to be in MNI152 space. The mask is the functional mask image for that subject and session.

         .. code-block:: python

            pial = self.anatomical_dir.format(sub) + f'/{sub}_space-32k_hemi-{hem}_pial.surf.gii'
            white = self.anatomical_dir.format(sub) + f'/{sub}_space-32k_hemi-{hem}_white.surf.gii'
            mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'

* Defined in ``DataSetNative``

    * ``MNI152NLin6Asym``: We use the individual deformation into MNI space and the functional mask image for that subject and session.

            .. code-block:: python

                deform = self.anatomical_dir.format(sub) + f'/{sub}_space-MNI_xfm.nii'
                mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'

* Defined in ``DataSetMNIVol``

    * ``MNI152NLin6Asym``: We use no deformation (already in correct space) and the functional mask image for that subject and session.

            .. code-block:: python

                deform = None
                mask = self.estimates_dir.format(sub) + f'/{ses_id}/{sub}_{ses_id}_mask.nii'


Data aggregation
----------------


is implemented in the dataset module in

``dataset.extract_all``

The steps are to
* Define the altas with `am.get_atlas``
* Get the participant list with `mydataset.get_participants``
* Define the individual atlas maps (when using Native space data)
* Get the file names - this is automatically done from `_reginfo.tsv` files
* call `condense_data` to prewhiten and average the data in the desired format
* Save the extracted data to a cifti file as derivatives/subj/data/subj_space-S_condHalf.dscalar.nii







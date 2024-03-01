Data Extraction
###############

Extraction is the process by which the estimates or timeseries data from a dataset are read out. The result will be stored in a cifti-file, with a Spatial and Data dimension.

Spatial resampling
------------------
Along the spatial dimension, extraction reads out the data from the source space (stored in ``dataset.space``) and deforms and sample them into the atlas space (``atlas.space``). This mapping between different spaces is store in an ``AtlasMap`` object, and can rely either on a subject-specific deformation, a atlas-to-atlas deformation, or even a combination of both. 

The specific mapping rules for the dataset are defined in the dataset function ``get_atlasmap``. We have implemented a set of dafault mapping rules - however you can overwrite this behavior by inherting the dataset object. 

Default behavior: 






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







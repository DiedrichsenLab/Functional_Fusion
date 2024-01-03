Data Extraction
###############

Extraction is implemented in the dataset module in 

``dataset.extract_all``

The steps are to 
* Define the altas with `am.get_atlas``
* Get the participant list with `mydataset.get_participants``
* Define the individual atlas maps (when using Native space data)
* Get the file names - this is automatically done from `_reginfo.tsv` files 
* call `condense_data` to prewhiten and average the data in the desired format
* Save the extracted data to a cifti file as derivatives/subj/data/subj_space-S_condHalf.dscalar.nii







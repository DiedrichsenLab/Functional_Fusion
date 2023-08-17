Overview
========

Each dataset is stored in a separate directory in the `basedir` directory. Depending to the dataset, the preprocessed time series or individual effect-size estimates are stored.
These data are usually stored in the Native subject space, but could also be stored in a volumetric group space or even in CIFTI files (surface + subcortical ROIs).

To use this data there are two steps:

* Data Extraction: Pull the data in a specific atlas space, defined by an `atlas`. The resulting data is stored in the `data` directory for each participant.
* Data Analysis: After the data is extracted, you can simply load the data with `dataset.get_data()` which gives you a `(nsubj x  nfeatures x voxel/vertices)` tensor. You can then perform any analysis on this data.


Overview
========

The Functional_Fusion repository is designed to bring data from many different fMRI experiments into a common analysis space - flexibly and efficiently.
This allows the aggregation of data to build better and bigger models of brain function.

.. image:: _static/extraction.png
  :alt: overview

Each dataset may be stored in different locations and in different format. The user can communicate with each dataset over a ``DataSet`` class, which knows where to find the data and how to load it. This allows the user to interact with each dataset in a uniform way.
We provide a number of common DataSet structures that each assume a specific way of data organization, usually following a BIDS-derivative structure.

The Analysis spaces (or regions of interest) are defined by an ``Atlas``. While we have predefined some common Atlases or regions, you can define your own easily. The framework supports both Volume- and Surface-based atlases. The mapping between each Atlas and each subject in each dataset is determined by an ``AtlasMap``, which allows you to extract the data from each subject without reslicing it into a common space. This is especially usefule if you want to get the time-series data from a specific region in the native space of the subject. Extracted data can be storted in CIFTI-files for further use, so you do not have to redo the extraction step every tie you want to use the data.

There are three main steps to using this framework on new data:

* Data Import: Bringing the data into the common framework. This includes the import of the preprocessed time series, or the import of the contrast estimates.
* Data Extraction: Pull the data in a specific atlas space, defined by an ``atlas``. The resulting CIFTI-file can be stored in a ``data`` directory for each dataset for quick retrieval.
* Data Analysis: After the data is extracted, you can simply load the data with ``dataset.get_data()``. which gives you a ``(nsubj x  nfeatures x voxel/vertices)`` tensor. You can then perform any analysis on this data.

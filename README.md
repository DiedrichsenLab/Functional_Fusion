Functional_Fusion
====
Diedrichsen Lab, Western University

This repository is the structure and preprocessing code for the multi-data set project in the Diedrichsenlab. 

Installation and dependencies
------
This project depends on several third party libraries, including: 

[numpy](https://numpy.org/) (version>=1.22.2), [nilearn](https://nilearn.github.io/stable/index.html) (version>=0.9.0), ...

	pip install numpy nilearn ...

Or you can install the package manually from the original binary source as above links.	

Structures of the project
------
### Overall structure
![ScreenShot](docs/structure.png)

As can be seen in above diagram, we first get the individual raw imaging data from a specific anatomical 
structure, then pass the raw data to `dataset` class for further processing. In parallel, we also need to
decide an `atlas` to be used and this `atlas` will be the input for the `atlas map` function to find
the vertices -> voxels mapping list. Lastly, the `get_data()` function in `dataset` class outputs the 
processed data `Y` in desired format `(N * P)` where `N` is the number of tasks `P` is the number of brain
locations.

### 1. Anatomical structure

The anatomical structure indicates the brain areas to be studied with a structural template (i.e MNI 152).
The current interested brain areas include but not limited to cerebellum and cortex.

### 2. Atlas map

The mapping function from the vertices of surface-based template to the list of connecting voxels. This 
is the main function to integrate and connect multiple dataset into the same format, which contains many
lookup table (potentially python dictionary, but I think we can do better) and find the correct one 
given by the current atlas being used and imaging resolution.

### 3. Data Set class

The Data Set class `dataset` is designed to be the entry of getting the data in standard format. It 
reads the input of raw individual data and other parameters for further use. The class function 
`get_data()` is to get the final processed data matrix `Y` after minimally preprocessed data.
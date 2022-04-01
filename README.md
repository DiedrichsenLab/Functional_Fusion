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
![ScreenShot](docs/structure.png)

#### 1. Data Set class
As can be seen in above diagram, the integration across data sets is achieved through  `DataSet` objects, with each data set having an instantiation of a subclass. The `DataSet` has access to the subject list and the individual preprocessed imaging data. The main function of the Data set class is the  `get_data()` function, which provides the 
processed data `Y` in desired format `(N * P)` where `N` is the number of measurements (tasks, etc) `P` is the number of brain locations for a specific subject. 

### 2. Atlas map

To read out different data sets in a consistent anatomical location, we need to
decide an `Atlas` to be used. For each subejct, we need to have an `AtlasMap`, which provides the mapping function from the raw data space (per subject) to the common atlas space. We will have subclasses `AtlasMapMNI`, `AtlasMapFS32K`, and `AtlasMapSUIT3`. For each subject (and hence Dataset), there will be a seperate Instantiation (Object) of this class. If two dataset share the same raw data space for the same subject, they can rely on the same AtlasMap. 

### 3. Altas
Each group atlas also has some subject / study independent behaviors that will be implemented in the `Atlas` Class. Subclasses indicates the brain areas to be studied with a structural template (i.e MNI 152).

## Directory structure for derivatives 

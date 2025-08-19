Datasets
========

An increasing number of dataset are available within the Functional Fusion framework. We are focussing on datasets with a large variety of tasks within the same subjects. These datasets provide an especially powerful way to characterize brain organization. Below is a description of datasets that have been brought into the framework, some of which are openly available for download, while others are available upon request.

Contributing datasets
---------------------
You can bring your imaging dataset into the required format and use the framework for analysis (see :ref:`import`). If you want to make your datasets available, please contact us. We are happy to test your dataset, and list it here. Most of the processed datasets are stored on Zenodo, giving each preprocessed version a citable DOI. Zenodo also allows you to determine exactly with whom, when, and under which conditions, you want to share your dataset.

Using dataset
-------------
After download / organizing a dataset, you should place the Dataset folder (i.e. 'MTDB') into the Functional_Fusion base directory.
For automatic loading of the dataset, you need to add the following information into a file called ``dataset_description.tsv`` in the base directory. This is a tab-delimited text file with the following columns:
* name: Name of the dataset (e.g. 'MTDB')
* dir_name: Path to the dataset directory. Start with ``/`` for absolute path, or set relative to the Functional_Fusion base directory.
* class_name: Name of Dataset class (i.e. ``DataSetNative`` or ``DataSetMNIVol``)


Datasets that are available in the framework
------------------------------------------------

Multi-domain task battery (King et al., 2019)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
|    **Description**: A large battery of tasks that cover a wide range of cognitive domains, including perception, action, and memory. Includes hrs task data from 2 different tasks sets, collected over 4 days of scanning. Resting state data for 19 subjects is also available (20min).
|    **Reference**: King, M., Hernandez-Castillo, C. R., Poldrack, R. A., Ivry, R. B., & Diedrichsen, J. (2019). Functional boundaries in the human cerebellum revealed by a multi-domain task battery. Nature Neuroscience, 22(8), 1371–1378.
|    **No of subjects**: 24
|    **No of task conditions**: 47
|    **DOI/Link**: 10.5281/zenodo.16788784
|    **Availability**: Openly available
|    **Maintainer**: joern.diedrichsen@googlemail.com

Working memory and finger tapping (Shahshahani et al., 2024)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
|    **Description**: Working memory (digit span forward and backwards) and finger tapping with different speeds and forces.
|    **Reference**: Shahshahani L, King M, Nettekoven C, Ivry R, Diedrichsen J. Selective recruitment of the cerebellum evidenced by task-dependent gating of inputs. Elife. 2024;13: RP96386.
|    **No of subjects**: 16
|    **No of task conditions**: 17
|    **DOI/Link**: 10.5281/zenodo.16788634
|    **Availability**: Openly available
|    **Maintainer**: joern.diedrichsen@googlemail.com

HCP unrelated 100
^^^^^^^^^^^^^^^^^^^^^^^^^^
|    **Description**: Human connectome project, unrelated 100 dataset. Both task-based an resting state sessions for the same subjects.
|    **Reference**: Barch DM, Burgess GC, Harms MP, Petersen SE, Schlaggar BL, Corbetta M, et al. Function in the human connectome: task-fMRI and individual differences in behavior. Neuroimage. 2013;80: 169–189.
     **No of subjects**: 50
|    **No of task conditions**: 24
|    **DOI/Link**: 10.5281/zenodo.16903949
|    **Availability**: Openly available
|    **Maintainer**: jdiedric@uwo.ca

Individual Brain Charting (IBC, Pinho et al., 2021)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
|    **Description**: A deep-phenotyping dataset with a large number of task. Each group of tasks is measured in a separate imaging session.
|    **Reference**: Pinho AL, Amadon A, Ruest T, Fabre M, Dohmatob E, Denghien I, et al. Individual Brain Charting, a high-resolution fMRI dataset for cognitive mapping. Scientific data. 2018;5.
  |  **No of subjects**: 12
|    **No of task conditions**: 201
|    **DOI/Link**:
|    **Availability**: Openly available
|    **Maintainer**: joern.diedrichsen@googlemail.com


Nishomoto (Nakai and Nishimoto)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
|    **Description**: A dataset with 104 tasks, each presented once during each imaging run.
|    **Reference**: Nakai T, Nishimoto S. Quantitative models reveal the organization of diverse cognitive functions in the brain. Nat Commun. 2020;11.
  |  **No of subjects**: 6
|    **No of task conditions**: 104
|    **DOI/Link**:
|    **Availability**: Openly available
|    **Maintainer**: joern.diedrichsen@googlemail.com

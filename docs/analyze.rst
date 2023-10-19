Using Datasets
##############

To use Functional fusion framework you need access to a folder that holds the datasets and atlases.

.. code-block:: python

    import Functional_Fusion.dataset as ds
    import nibabel as nb
    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'


Loading Data
------------
Loading the data to get a ``n_subj x n_cond x n_voxels`` tensor:

.. code-block:: python

    X,info,dataset_obj = ds.get_dataset(base_dir,
            dataset='MDTB',
            atlas='fs32k',
            sess='all',
            type='CondAll')

You can specify subset of sessions, subjects, etc.


Using Datasets
==============

To use Functional fusion framework you need access to a folder that holds the datasets and atlases.

.. code-block:: python

    import Functional_Fusion.dataset as ds
    import Function_Fusion.utils as ut
    import nibabel as nb
    base_dir = ut.get_base_dir()


Loading Data
------------
Loading the data to get a ``n_subj x n_cond x n_voxels`` tensor:

.. code-block:: python

    X,info,dataset_obj = ds.get_dataset(base_dir,
            dataset='MDTB',
            atlas='fs32k',
            sess='all',
            type='CondRun')

You can specify subset of sessions, subjects, etc.

Aggregating data
----------------
If you want to average data across runs, you can use the get_dataset function with `type='CondAll'`, or alternatively aggregate the data the following way:

.. code-block:: python

    cinfo,C = ds.agg_data(info,['cond_num_uni'],['run','half','reg_num','names'])
    cdata = np.linalg.pinv(C) @ data

Group averaging data
--------------------
To produce the group-averaged dscalar files for a specfic atlas space and data type, just call:

.. code-block:: python

    dataset_obj.group_average_data(atlas='MNISymDentate1',ses_id='ses-s1',type='CondRun')

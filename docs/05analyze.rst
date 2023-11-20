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

Assessing reliability
---------------------

We can quickly calculate the reliability of the data across subjects and within subjects across sessions: 

.. code-block:: python

    for i,ses in enumerate(dataset.sessions):
        data,info = dataset.get_data('MNISymC3',ses,'CondHalf')
        rw = ds.reliability_within_subj(data,
                        part_vec=info.half,
                        cond_vec=info.reg_id)
        RW[:,i] = rw.mean(axis=1)
        RB[:,i] = ds.reliability_between_subj(data,
                        cond_vec=info.reg_id)
    pass



Decomposing variances
--------------------- 

A more general way of looking at reliability within and across subjects is to decompose the variance into components related to the group, the subject and the measurement error.

We have activity patterns or sets of activity patterns in a collection of :math:`NxP` matrices. You can have P=1 (activity profile), N=1 (activity pattern) or both >1 (sets or activity pattern). Let's call these :math:`\mathbf{y}_{i,j}`.  These can be either uncenterd (relative to rest) or centered to the mean of each voxel (condition difference). 

.. math::
    \mathbf{y}_{i,j} = \mathbf{g} + \mathbf{s}_i + \boldsymbol{\epsilon}_{i,j}


We want to estimate the norms (sum of squares) attached to these terms. 

.. math::
    \begin{array}{c}
    v_{g} = E(\mathbf{g}^T\mathbf{g})\\
    v_{s} = E(\mathbf{s}^T\mathbf{s})\\
    v_{\epsilon} = E(\mathbf{\epsilon}^T\mathbf{\epsilon})
    \end{array}

We assume that :math:`\mathbf{g}`, :math:`\mathbf{s}`,and :math:`\mathbf{\epsilon}` are mutually independent, i.e. :math:`E(\mathbf{g}^T\mathbf{s})=0`. 

First, let write out the expected values of the cross-subj, cross-run and within-run sums of squares. 



Across subjects: 

.. math::
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{k,l}) = v_{g}

Within subject, across runs: 

.. math::
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,k}) = v_{g} + v_{s}

Within observation:

.. math::
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,k}) =  v_{g} + v_{s} + v_{\epsilon}

To develop estimators for these quantities we replace the Expectation with the mean **over all possible pairings**.

ADD CODE EXAMPLE HERE


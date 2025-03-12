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

Group averarging data
---------------------
To produce the group-averaged dscalar files for a specfic atlas space and data type, just call:

.. code-block:: python

    dataset_obj.group_average_data(atlas='MNISymDentate1',ses_id='ses-s1',type='CondRun')


Assessing reliability within and across individuals
---------------------------------------------------


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



Decomposing variances (checking reliability)
--------------------------------------------

A more general way of looking at reliability within and across subjects is to decompose the variance into components related to the group, the subject and the measurement error.

We have activity patterns or sets of activity patterns in a collection of :math:`NxP` matrices. You can have P=1 (activity profile), N=1 (activity pattern) or N>1 and P>1 (activity matrix). Now we have measurements of these activity patterns from different subjects (i) and different runs or repetitions (j).
Let's call these :math:`\mathbf{y}_{i,j}`.  These can be either uncentered (relative to rest) or centered to the mean of each voxel (condition difference).

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
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,j}) =  v_{g} + v_{s} + v_{\epsilon}

To develop estimators for these quantities we replace the Expectation with the mean **over all possible pairings**.

ADD CODE EXAMPLE HERE


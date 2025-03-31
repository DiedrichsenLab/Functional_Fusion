Reliability
===========

In general the activity data from a subject is a :math:`NxP` matrix with *N* being the number of trials and *P* the number of voxels. The conditions are indicated by ``cond_vec`` and the partitions (independent measures) by ``part_vec``.  We have *K* conditions and *M* partitions.
The data can have ``P=1`` (activity profile), ``K=1`` (activity pattern) or ``N>1`` and ``P>1`` (activity matrix). In most cases, the ``data`` is a 3D-array with dimensions ``n_subjects x N x P``.


Reliability within individual
-----------------------------

As a measure of functional signal-to-noise ratio, we can calculate the reliability of activity patterns across runs. The measured activity pattern from subject *i* in run *j* (:math:`\mathbf{y}_{i,j}`) can be thought of consisting of the signal in subject i :math:`\mathbf{s}_{i}` and measurement noise :math:`\boldsymbol{\epsilon}_{i,j}`.

.. math::
    \mathbf{y}_{i,j} = \mathbf{s}_i + \boldsymbol{\epsilon}_{i,j}

For each subject, we can calculate the variances (or second moments) of these patterns across voxels and / or conditions:

.. math::
    \begin{array}{c}
    v_s = \mathbf{s}_i^T\mathbf{s}_i\\
    v_{\epsilon} = E(\mathbf{\epsilon}^T\mathbf{\epsilon})
    \end{array}

The correlation of the measured activity patterns across runs then is: 

.. math::
    r_{run} = \frac{v_{s}}{v_{s} + v_{\epsilon}}

This reliability measure is the reliability of one measurement (run). The (theoretical) reliability of the mean activity pattern (across *N* measurement runs) then would be: 

.. math::
    r_{whole} = \frac{v_{s}}{v_{s} + v_{\epsilon}/N}

With a bit of algebra, you can calculate the theoretical reliability of the mean activity pattern from the reliability across runs: 

.. math::
    r_{whole} = \frac{r_{run} N}{r_{run} (N-1) +1}

This is basically the same idea as Cronbach's alpha. 

We can quickly calculate the reliability of the data within each subjects across different runs:

.. code-block:: python

    data,info = dataset.get_data('MNISymC3',ses,'CondRun')
    rw = reliability.within_subj(data,
            part_vec=info.run,
            cond_vec=info.reg_id,
            separate='None',
            subtract_mean=True)

rw will be a vector of reliability values for each subject, if multiple subjects are present. 

Reliability across individuals
------------------------------ 
Similarly, we can calculate the reliability of the mean activity patterns across subjects. The underlying model here is that all subjects have a shared group pattern :math:`\mathbf{g}` and a subject-specific pattern :math:`\mathbf{s}_i`. The logic is the same as above, but now we are looking at the reliability of the group pattern across subjects.

.. math::
    \mathbf{y}_{i} = \mathbf{g} + \mathbf{s}_i 

.. code-block:: python

    # Get the data per subject: Note that the data is averaged across partitions in any case:  
    data,info = dataset.get_data('MNISymC3',ses,'CondAll')
    rw = reliability.within_subj(data,
            cond_vec=info.reg_id,
            separate='None',
            subtract_mean=True)

Decomposing variances (reliability within and across subject)
----------------------------------------------------------------------

A more general way of looking at reliability within and across subjects at the same time is to decompose the variance into components related to the group, the subject and the measurement error.

We have measurements of different subjects (i) and different runs or repetitions (j), :math:`\mathbf{y}_{i,j}`. We can think of these measurements as consisting of a group pattern :math:`\mathbf{g}`, a subject-specific pattern :math:`\mathbf{s}_i`, and measurement noise :math:`\boldsymbol{\epsilon}_{i,j}`.

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

.. code-block:: python

    # To get the group,subject, and run decomposition (fSNR) of the data:  
    data,info = dataset.get_data('MNISymC3',ses,'CondAll')
    var = rel.decompose_subj_group(data,
                cond_vec,
                part_vec,
                separate='subject_wise',
                subtract_mean=True)

Decomposing variances (subject-specific scaling)
------------------------------------------------
In situations where the data of subjects are scaled differently, the variance of group, individual, and measurement noise is also scaled accordingly.

.. math::
    \mathbf{y}_{i,j} = (\mathbf{g} + \mathbf{s}_i + \boldsymbol{\epsilon}_{i,j}) \times \text{scale}_{i}

Therefore, the variance components will be scaled too:

.. math::
    \begin{array}{c}
    v_{g,i} = v_{g} \times \text{scale}_{i}^2\\
    v_{s,i} = v_{s} \times \text{scale}_{i}^2\\
    v_{\epsilon,i} = v_{\epsilon} \times \text{scale}_{i}^2
    \end{array}

Considering the scaling factors, we can rewrite the expected values of the cross-subj, cross-run and within-run sums of squares:

Across subjects:

.. math::
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{k,l}) = \text{scale}_{i} \times \text{scale}_{k} \times v_{g}

Within subject, across runs:

.. math::
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,k}) = \text{scale}_{i}^2 (v_{g} + v_{s})

Within observation:

.. math::
    E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,j}) =  \text{scale}_{i}^2 (v_{g} + v_{s} + v_{\epsilon})

.. code-block:: python

    # To get the group,subject, and run decomposition (fSNR) of the data (subject-specific scaled):  
    var = decompose_variance_scaled(data)

Mean substraction
-----------------
All reliability functions have an optional input parameter ``subtract_mean``. The default setting is ``subtract_mean=True``. This means that the mean activity in each voxel in each partition (across conditions) is subtracted out before computing the variances or correlations. Thus reliability and noise estiamtes are based on **differences between conditions** but do not reflect the activation of a voxel relative to the implicit baseline. 

If you set ``subtract_mean=False``, the mean activity in each voxel in each partition is not subtracted out. This means that the reliability and noise estimates are based both on the **mean activity pattern**  across conditions, as well as **differences between conditions**. Usually, this leads to much higher reliabilities, as the mean activity pattern if often stronger than the differences between conditions. 

Separating the analysis by voxel and condition
----------------------------------------------
For all functions, you can specify the parameter ``separate``. The default setting is ``separate='None'``. This means that the reliability is calculated across all voxels and conditions. If you set ``separate='voxel_wise'``, the reliability is calculated for each voxel separately. If you set ``separate='cond_wise'``, the reliability is calculated for each condition separately. 

Leave-one-out reliability
-------------------------
We also provide a function that calculates the correlation of the pattern within each run with the average pattern for the other runs (``reliability.within_subj_loo``). The separate measures for each run are useful to spot a run that has bad signal-to-noise, or for which there was an error in the processing. 

Similarly, we also provide a function that calculates the correlation of the pattern for each subject with the average pattern for the other subjects (``reliability.between_subj_loo``). The subject-specific measure can be used for spotting outlier subjects or subjects for which an error occurred. The average of the loo-correlation can also serves as lower noise-ceiling for group models. 

*Note that the reliability measures across runs (or across subjects) are not strictly independent, so care needs to be taken when using these  measures in statistical tests.* 

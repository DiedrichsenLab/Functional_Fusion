# Decomposition of patterns into group, individual and noise

We have activity patterns or sets of activity patterns in a collection of $NxP$ matrices. You can have P=1 (activity profile), N=1 (activity pattern) or both >1 (sets or activity pattern). Let's call these $\mathbf{y}_{i,j}$.  These can be either uncenterd (relative to rest) or centered to the mean of each voxel (condition difference). 

$$
\mathbf{y}_{i,j} = \mathbf{g} + \mathbf{s}_i + \boldsymbol{\epsilon}_{i,j}
$$

We want to estimate the norms (sum of squares) attached to these terms. 
$$
v_{g} = E(\mathbf{g}^T\mathbf{g})\\
v_{s} = E(\mathbf{s}^T\mathbf{s})\\
v_{\epsilon} = E(\mathbf{\epsilon}^T\mathbf{\epsilon})
$$

We assume that $\mathbf{g}$, $\mathbf{s}$,and $\mathbf{\epsilon}$ are mutually independent, i.e. $E(\mathbf{g}^T\mathbf{s})=0$. 

First, let write out the expected values of the cross-subj, cross-run and within-run sums of squares. 



Across subjects: 
$$
E(\mathbf{y}_{i,j}^T\mathbf{y}_{k,l}) = v_{g}
$$
Within subject, across runs: 
$$
E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,k}) = v_{g} + v_{s}
$$
Within observation: 
$$
E(\mathbf{y}_{i,j}^T\mathbf{y}_{i,k}) =  v_{g} + v_{s} + v_{\epsilon}
$$

To develop estimators for these qualtities we replace the Expectation with the mean **overall all possible pairings**


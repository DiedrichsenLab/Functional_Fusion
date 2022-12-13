import numpy as np
from numpy.linalg import inv
import nibabel as nb


def sq_eucl_distances(coordA,coordB):
    D = coordA.reshape(3,-1,1)-coordB.reshape(3,1,-1)
    D = np.sum(D**2,axis=0)
    return D

def surf_from_cifti(ts_cifti,
                    struct_names=['CIFTI_STRUCTURE_CORTEX_LEFT',
                                  'CIFTI_STRUCTURE_CORTEX_RIGHT'],
                    mask_gii=None):
        """
        Gets the time series of cortical surface vertices (Left and Right)
        Args:
            ts_cifti (cifti obj) - cifti object of time series
            struct_names (list): the struct name of left and right cortex
            mask_gii: (list of Obj.) the mask gifti object for each hemisphere
                                     if given. Default is None,
                                     indicating no mask for return.
                      (list of str.) the file path of mask gifti for each
                                     hemisphere if given. Default is None,
                                     indicating no mask for return.
        Returns:
            cii (cifti object) - contains the time series for the cortex
        """
        # get brain axis models
        bmf = ts_cifti.header.get_axis(1)
        # print(dir(bmf))
        # get the data array with all the time points, all the structures
        ts_array = ts_cifti.get_fdata()
        ts_list = []
        for idx, (nam,slc,bm) in enumerate(bmf.iter_structures()):
            # just get the cortical surfaces
            if nam in struct_names:
                values = np.full((ts_array.shape[0],bmf.nvertices[nam]),np.nan)
                # get the values corresponding to the brain model
                values[:,bm.vertex] = ts_array[:, slc]
                if mask_gii is not None:
                    Xmask = mask_gii[struct_names.index(nam)]
                    if isinstance(Xmask, str):
                        Xmask = nb.load(Xmask).agg_data()
                    elif isinstance(Xmask, object):
                        Xmask = Xmask.agg_data()
                    else:
                        raise ValueError("The input mask_gii must be either a string of path"
                                         "or a nibabel gii image object!")
                    values = values[:, np.where(Xmask>0)[0]]
                ts_list.append(values)
            else:
                break
        return ts_list

def zstandarize_ts(X):
    X = X - X.mean(axis = 0, keepdims = True)
    X = X / np.sqrt(np.nansum(X**2, axis=0)/X.shape[0])
    return X

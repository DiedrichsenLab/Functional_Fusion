import numpy as np 
from numpy.linalg import inv

def affine_transform(x1, x2, x3, M):
    """
    Returns affine transform of x

    Args:
        x1 (np-array):
            X-coordinate of original
        x2 (np-array):
            Y-coordinate of original
        x3 (np-array):
            Z-coordinate of original
        M (2d-array):
            4x4 transformation matrix

    Returns:
        x1 (np-array):
            X-coordinate of transform
        x2 (np-array):
            Y-coordinate of transform
        x3 (np-array):
            Z-coordinate of transform

    """
    y1 = M[0,0]*x1 + M[0,1]*x2 + M[0,2]*x3 + M[0,3]
    y2 = M[1,0]*x1 + M[1,1]*x2 + M[1,2]*x3 + M[1,3]
    y3 = M[2,0]*x1 + M[2,1]*x2 + M[2,2]*x3 + M[2,3]
    return (y1,y2,y3)

def affine_transform_mat(x, M):
    """
    Returns affine transform of x

    Args:
        x (np-array):
            3xN array for original coordinates 
        M (2d-array):
            4x4 transformation matrix
    Returns:
        y (np-array):
            3xN array pof X-coordinate of transformed coordinaters 
    """
    y = M[0:3,0:3] @ x + M[0:3,3:]
    return (y)


def coords_to_linvidxs(coords,vol_def,mask=False):
    """
    Maps coordinates to linear voxel indices

    INPUT:
        coords (3xN matrix or Qx3xN array):
            (x,y,z) coordinates
        vol_def (nibabel object):
            Nibabel object with attributes .affine (4x4 voxel to coordinate transformation matrix from the images to be sampled (1-based)) and shape (1x3 volume dimension in voxels)
        mask (bool):
            If true, uses the mask image to restrict voxels (all outside = -1)
    OUTPUT:
        linvidxs (N-array or QxN matrix):
            Linear voxel indices
        good (bool) boolean array that tells you whether the index was in the mask
    """
    mat = inv(vol_def.affine)

    # Check that coordinate transformation matrix is 4x4
    if (mat.shape != (4,4)):
        raise(NameError('Error: Matrix should be 4x4'))

    rs = coords.shape
    if (rs[-2] != 3):
        raise(NameError('Coords need to be a (Kx) 3xN matrix'))

    # map to 3xP matrix (P coordinates)
    ijk = mat[0:3,0:3] @ coords + mat[0:3,3:]
    ijk = np.rint(ijk).astype(int)
    
    if ijk.ndim<=2:
        i = ijk[0]
        j = ijk[1]
        k = ijk[2]
    elif ijk.ndim==3:
        i = ijk[:,0]
        j = ijk[:,1]
        k = ijk[:,2]

    # Now set the indices out of range to 
    good = (i>=0) & (i<vol_def.shape[0]) & (j>=0) & (j<vol_def.shape[1]) &  (k[2]>=0) & (k[2]<vol_def.shape[2])
    
    linindx = np.ravel_multi_index((i,j,k),vol_def.shape,mode='clip')

    if mask:
        M=vol_def.get_fdata().ravel()
        good = good & (M[linindx]>0)
    
    return linindx, good

def sq_eucl_distances(coordA,coordB): 
    D = coordA.reshape(3,-1,1)-coordB.reshape(3,1,-1)
    D = np.sum(D**2,axis=0)
    return D
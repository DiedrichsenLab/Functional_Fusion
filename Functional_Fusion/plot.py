import matplotlib.pyplot as plt
import numpy as np
import nitools as nt
import nibabel as nb
import Functional_Fusion.atlas_map as am
import Functional_Fusion.util as ut
import nibabel as nb
import nilearn.plotting as nlp
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap

import nitools as nt
from numpy.linalg import inv

def ortho(data, voxel, fig=None, cursor=False, background=None, **kwargs):
    """Simple orthographic plot of a 3D array using matplotlib.

    Args:
        data (nd-array): 3D numpy array of image data 
        arg voxel (list): XYZ coordinates for each slice
        fig (plt.fig, optional):   Existing figure and axes for overlay plotting
        cursor (bool): Show a cursor at the `voxel` (default False)
        background (str): Color of background for slices (default None)
        **kwargs:  All other arguments are passed through to the `imshow` function.

    Returns:
        fig:  figure handle
        ax_x: Axis for x-slice (coronal)
        ax_y: Axis for y-slice (sagittal)
        ax_z: Axis for z-slice (axial)
    """

    voxel = [int(round(v)) for v in voxel]

    data = np.asanyarray(data, dtype=float)
    data[data <= 0] = np.nan

    x, y, z = voxel
    xslice = np.flipud(data[x, :, :].T)
    yslice = np.flipud(data[:, y, :].T)
    zslice = np.flipud(data[:, :, z].T)

    if fig is None:
        fig = plt.figure()
        xax = fig.add_subplot(1, 3, 1)
        yax = fig.add_subplot(1, 3, 2)
        zax = fig.add_subplot(1, 3, 3)
    else:
        fig, xax, yax, zax = fig

    # Plot black square the size of each slice
    if background == 'black':
        for ax, slc in zip((xax, yax, zax), (xslice, yslice, zslice)):
            ax.imshow(np.zeros_like(slc), cmap='gray')
    
    xax.imshow(xslice, **kwargs)
    yax.imshow(yslice, **kwargs)
    zax.imshow(zslice, **kwargs)

    if cursor:
        cargs = {'color': (0, 1, 0), 'linewidth': 1}
        xax.axvline(y, **cargs)
        xax.axhline(data.shape[2] - z, **cargs)
        yax.axvline(x, **cargs)
        yax.axhline(data.shape[2] - z, **cargs)
        zax.axvline(x, **cargs)
        zax.axhline(data.shape[1] - y, **cargs)

    for ax in (xax, yax, zax):
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout(pad=0)

    # Remove border from image
    for ax in (xax, yax, zax):
        for side in ('top', 'bottom', 'left', 'right'):
            ax.spines[side].set_visible(False)
    
    return (fig, xax, yax, zax)

def plot_dentate(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-31,-33,-35,-37,-39,-42],
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for detate nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-35,-37,-39,-42].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymDentate1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_dentate.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right dentate at the voxel coordinates
    c1 = np.array([[-25,-70,-43],[7,-70,-43]]).T # Lower left corner of image
    c2 = np.array([[-7,-43,-20],[25,-43,-20]]).T # Upper right corner of image
    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            fig = plt.figure(figsize=(2,10),facecolor='black')
        gridspec = fig.add_gridspec(6, 2,hspace=0.1,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=[z],
            bg_img=bg[i],
            black_bg=True,
            axes=axes[j,i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes


def plot_thalamus(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-31,-33,-7,-37,-38,7],   #-31, -33, -7, -37, -33, 7
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for detate nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-35,-37,-39,-42].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymThalamus1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_thalamus.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right dentate at the voxel coordinates
    c1 = np.array([[-34,-39,-8],[0,-39,-8]]).T # Lower left corner of image
    c2 = np.array([[1,5,19],[33,5,19]]).T # Upper right corner of image
    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            fig = plt.figure(figsize=(15,6),facecolor='black')
            #fig = plt.figure(figsize=(2,10),facecolor='black')
        #gridspec = fig.add_gridspec(6, 2,hspace=0.1,wspace=0.1)
        gridspec = fig.add_gridspec(1, 2,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=12,
            bg_img=bg[i],
            black_bg=True,
            #axes=axes[j,i],
            axes = axes[i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes

def plot_thalamus2(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-7,-5,-3,1, 8,12,15],  
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for detate nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-35,-37,-39,-42].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymThalamus1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_thalamus.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right dentate at the voxel coordinates
    c1 = np.array([[-34,-39,-10],[0,-39,-10]]).T # Lower left corner of image
    c2 = np.array([[0,5,18],[34,5,18]]).T # Upper right corner of image
    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            #fig = plt.figure(figsize=(15,6),facecolor='black')
            fig = plt.figure(figsize=(2,10),facecolor='black')
        gridspec = fig.add_gridspec(7, 2,hspace=0.1,wspace=0.1)
        #gridspec = fig.add_gridspec(1, 2,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=[z],
            bg_img=bg[i],
            black_bg=True,
            axes=axes[j,i],
            #axes = axes[i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes

def plot_olive(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-31,-29,-54,-31,-29,-47],
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for inferior olivary nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-60,-37,-39,-53].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymOlive1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_olive.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right olive at the voxel coordinates
    
    c1 = np.array([[-11,-45,-56],[0,-45,-56]]).T # Lower left corner of each image
    c2 = np.array([[0,-18,-29],[11,-18,-29]]).T # Upper right corner of each image

    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            fig = plt.figure(figsize=(15,6),facecolor='black')
            #fig = plt.figure(figsize=(2,10),facecolor='black')
        #gridspec = fig.add_gridspec(6, 2,hspace=0.1,wspace=0.1)
        gridspec = fig.add_gridspec(1, 2,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=12,
            bg_img=bg[i],
            black_bg=True,
            axes = axes[i],
            #axes=axes[j,i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes

def plot_olive2(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-47,-49, -50,-51,-52,-54],
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for inferior olivary nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-60,-37,-39,-53].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymOlive1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_olive.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right olive at the voxel coordinates
    
    c1 = np.array([[-11,-45,-56],[0,-45,-56]]).T # Lower left corner of each image
    c2 = np.array([[0,-18,-29],[11,-18,-29]]).T # Upper right corner of each image

    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            #fig = plt.figure(figsize=(15,6),facecolor='black')
            fig = plt.figure(figsize=(2,10),facecolor='black')
        gridspec = fig.add_gridspec(6, 2,hspace=0.1,wspace=0.1)
        #gridspec = fig.add_gridspec(1, 2,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=[z],
            bg_img=bg[i],
            black_bg=True,
            #axes = axes[i],
            axes=axes[j,i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes

def plot_rednucleus(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-13,-12,-10,-8,-7,-5],
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for inferior olivary nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-60,-37,-39,-53].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymRedNucleus1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_rednucleus.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right olive at the voxel coordinates
        
    c1 = np.array([[-20,-28,-16],[0,-28,-16]]).T # Lower left corner of each image
    c2 = np.array([[0,-6,-2],[20,-6,-2]]).T # Upper right corner of each image

    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            #fig = plt.figure(figsize=(15,6),facecolor='black')
            fig = plt.figure(figsize=(2,10),facecolor='black')
        gridspec = fig.add_gridspec(6, 2,hspace=0.1,wspace=0.1)
        #gridspec = fig.add_gridspec(1, 2,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=[z],
            bg_img=bg[i],
            black_bg=True,
            #axes = axes[i],
            axes=axes[j,i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes

def plot_pontine(data,
                 bg_img=None,
                 fig=None,
                 gridspec=None,
                 z_coords = [-47,-49, -50,-51,-52,-54],
                 cscale = [None,None],
                 cmap = 'cold_hot',
                 threshold = None):
    """Generate the plot for inferior olivary nucleus 
    For fine control of the visulization, see https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_img.html

    Args:
        data (ndarray): 
        bg_img (nifti1image): Background image. Defaults to None.
        fig (plt.figure): pre-specified matplotlib figure. 
        gridspec (Gridspec): A 6 x 2 Gridspec to plot the dentate data.
        z_coords (list): Z-coordinate slice to plot. Defaults to [-31,-33,-60,-37,-39,-53].
        cscale (list): [lower and upper] range for colorscale. None sets it to 2% percentile of data (asymmetric). 
        cmap (str, colormap, ndarray): Name, colormap, or Nx3 ndarray. Defaults to 'cold_hot'. 
        threshold (numeric): Single threshold: will plot data abs(y)> th
    
    Returns: 
        axes (array): 6 x2 array of subplots.
    """
    dn,_ = am.get_atlas('MNISymPontine1')
    if bg_img is None:
        adir = ut.default_atlas_dir
        bg_img = nb.load(adir + '/tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_res-1_olive.nii')
    
    # Project the functional data into the atlas space
    fcn_img = dn.data_to_nifti(data)

    # If cscale is not provided, calculate it from the data: 
    if cscale[0] is None:
        cscale[0] = np.percentile(data,2)
    if cscale[1] is None:
        cscale[1] = np.percentile(data,98)

    # Make a colormap from ndarray
    if isinstance(cmap,np.ndarray):
        cmap = ListedColormap(cmap)

    # Cut out the left and right olive at the voxel coordinates
    
    c1 = np.array([[-11,-45,-56],[0,-45,-56]]).T # Lower left corner of each image
    c2 = np.array([[0,-18,-29],[11,-18,-29]]).T # Upper right corner of each image

    v1 = nt.affine_transform_mat(c1,inv(bg_img.affine)).astype(int)
    v2 = nt.affine_transform_mat(c2,inv(bg_img.affine)).astype(int)

    bg = [] # Slice background data
    fc = [] # Sliced functional data
    for i in range(2):
        bg.append(bg_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])
        fc.append(fcn_img.slicer[v1[0,i]:v2[0,i]+1,v1[1,i]:v2[1,i]+1,v1[2,i]:v2[2,i]+1])

    # Initialize the figure and axes if not provided.
    if gridspec is None:
        if fig is None:
            #fig = plt.figure(figsize=(15,6),facecolor='black')
            fig = plt.figure(figsize=(2,10),facecolor='black')
        gridspec = fig.add_gridspec(6, 2,hspace=0.1,wspace=0.1)
        #gridspec = fig.add_gridspec(1, 2,wspace=0.1)
    
    # axes
    axes = gridspec.subplots()

    # Now use the nibabel plotting functions to plot the images
    for i in range(2):
        for j,z in enumerate(z_coords):
            nlp.plot_img(fc[i],
            display_mode="z",
            cut_coords=[z],
            bg_img=bg[i],
            black_bg=True,
            #axes = axes[i],
            axes=axes[j,i],
            threshold=threshold,
            vmin=cscale[0],
            vmax=cscale[1],
            cmap=cmap,
            annotate=False)

    return axes


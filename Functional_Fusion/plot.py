import matplotlib.pyplot as plt
import numpy as np


def ortho(data, voxel, fig=None, cursor=False, background=None, **kwargs):
    """Simple orthographic plot of a 3D array using matplotlib.

    :arg data:   3D numpy array
    :arg voxel:  XYZ coordinates for each slice
    :arg fig:    Existing figure and axes for overlay plotting
    :arg cursor: Show a cursor at the `voxel`
    All other arguments are passed through to the `imshow` function.

    :returns:   The figure and orthogaxes (which can be passed back in as the
                `fig` argument to plot overlays).
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

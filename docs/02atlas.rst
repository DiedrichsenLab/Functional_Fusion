Atlases, Regions, and Parcellations
===================================

The Atlas class defines the atlas space, the voxel resolution, and the set of voxels or vertices that will form the basis of the analysis.
Subclasses are ``AtlasVolumetric``, and ``AtlasSurface``. Each Atlas contains  :math:`P` locations (vertices or voxels) that are being sampled for the analysis.

Atlas Spaces
------------
Each Atlas is defined in a specific space. The atlas space defines how data is transformed between Native space and the atlas, or between different atlases. ``Functional_Fusion`` follows the convention of atlas space naming used in `Templateflow.org <https://www.templateflow.org/>`_.

The following spaces are defined so far:
    * ``MNI152NLin6Asym``: MNI152NLin6Asym space, the standard MNI space template in FSL
    * ``MNI152NLin2009cSym``: Symmetric version of the MNI152NLin2009c space
    * ``MNI152NLin2009cAsym``: Asymmetric version of the MNI152NLin2009c space
    * ``SUIT``: The SUIT space, a cerebellar template space
    * ``MNI152NLin2009cSymC``: Symmetric Cerebellum-only atlas (aligned to)
    * ``fs32k``: Symmetric surface-based space for the left and right hemispheres
    * ``fsaverage6``: Freesurfer fsaverage6 space - needs to be still implemented

Templates in different resolutions and the mask images that define the Atlases (see below) can be found in the directory ``Functional_Fusion/Atlases/tpl-<space>``.

Standard Atlases
----------------
Atlases are defined by a specific space and a mask image that defines the resolution and the voxels / vertices that are included in the atlas.
Standard atlases are indicated by short strings that indicate the Atlas. The following Atlases are defined so far:
    * ``SUIT3``:  Cerebellum in the SUIT space (3mm voxels)
    * ``SUIT2``:  Cerebellum in the SUIT space (2mm voxels)
    * ``SUIT1``:  Cerebellum in the SUIT space (1mm voxels)
    * ``MNISymC3``: Cerebellum in the MNI152NLin2009cSym space (3mm voxels)
    * ``MNISymC2``: Cerebellum in the MNI152NLin2009cSym space (2mm voxels)
    * ``fs32k``: The left and right hemispheres, using the symmetric medial wall mask
    * ``fs32kAsym``: The left and right hemispheres, using the asymmetric medial wall mask from the Human Connectome Project (HCP)
    * ``MNIAsymBg2``: Basal Ganglia in the MNI152NLin6Asym space (2mm voxels)
    * ``MNIAsymBg1``: Basal Ganglia in the MNI152NLin6Asym space (1mm voxels)

All these atlases are defined in ``Functional_Fusion/Atlases/atlas_description.json``.
In this file, each atlas is defined by a dictionary with the following fields:
    * ``type``: Class-name for the atlas object ("AtlasVolumetric", "AtlasSurface", ...)
    * ``dir``: Directory where the atlas files are stored: ``Functional_Fusion/Atlases/<dir>``
    * ``mask``: Name(s) of mask `NIFTI`/`GIFTI` files that define the atlas space.
    * ``space``: Space in which coordinates are defined (e.g. "MNI152NLin6AsymC")
    * ``structure``: `CIFTI`-structure name - needs to be a valid, predefined `CIFTI` structure (e.g. "cerebellum")

Optional fields:
    * ``res``: resolution (not used yet) - voxel locations are defined by mask image
    * ``normspace``: Norm space for suit.vol2surf -

Getting atlases
---------------
You can get an atlas by calling:

.. code-block:: python

    import Functional_Fusion.atlas_map as am

    atlas, ainf = am.get_atlas(atlas_str)

``atlas_str`` is a string that defines the atlas (see above). The function returns a tuple, where ``atlas`` is an Atlas object, and ``ainf`` is a dictionary with information about the atlas.

Defining new Atlases or Regions
-------------------------------
If an Atlas or Region is not defined in the repository, you can obtain a new Atlas object by calling the constructor of the Atlas class.

.. code-block:: python

    atlas = am.AtlasVolumetric('name', 'mask.nii', structure='thalamus', space='MNI152NLin6Asym')

``mask.nii`` is a `NIFTI` file that defines the voxels of the atlas, ``structure`` is the `CIFTI`-structure name, and ``space`` is the space in which the atlas is defined.


If an Atlas or Region is not defined in the repository, you can obtain a new Atlas object by calling the constructor of the Atlas class.

.. code-block:: python

    atlas = am.AtlasSurface('name', ['mask_L.gii', 'mask_L.gii'],structure=['cortex_left', 'cortex_right'], space='fs32k')

A surface-based atlas contains both hemisphere and therefore needs two `GIFTI` mask files and two structure names. You can also define an atlas on a single hemisphere.

Another way of defining an atlas is getting a ``subatlas`` of an existing atlas. This is often the case when you want to define an ROI for a specific region within an Atlas. Here are some examples for volumetric atlases:

.. code-block:: python

    # Volumetric subatlas examples
    atlas_vol = am.AtlasVolumetric(...)

    # Define a new atlas that is a spherical ROI in an existing volumetric atlas
        # center: 3-vector of center of the sphere in mm
        # raduis: radius of the sphere in mm
    region1 = atlas_vol.get_subatlas_sphere('region1', center, radius)

    # Define a new atlas that is defined by an image with an existing volumetric atlas
    region2 = atlas_vol.get_subatlas_image('region2', 'mask.nii')

And here is an example for a surface-based atlas:

.. code-block:: python

    # Surface-based subatlas example
    atlas_surf = am.AtlasSurface(...)

    # Define a new atlas that is defined by a set of nodes of an existing surface-based atlas (left hemisphere)
    atlas_left = atlas_surf.get_hemisphere(0)
    region3 = atlas_left.get_subatlas_image('region3', 'mask_L.gii')

    # Define a new atlas that is defined by a set of nodes of an existing surface-based atlas (right hemisphere)
    atlas_right = atlas_surf.get_hemisphere(1)
    region4 = atlas_left.get_subatlas_image('region3', 'mask_R.gii')

Reading data from NIFTI, GIFTI, or CIFTI files
----------------------------------------------

You can use the function ``read_data`` to get the data from a file in the atlas space.
Here an example for a volumetric atlas:

.. code-block:: python

    # Get the atlas
    atlas, ainf = am.get_atlas('SUIT3')
    # Read data from NIFTI file using linear interpolation (1)
    X = atlas.read_data('my_nift_file', interpolation=1)

And here one for a surface-based atlas:

.. code-block:: python

    # Get the atlas
    atlas,ainf = am.get_atlas('fs32k')
    # Read data from two GIFTI files for left and right hemisphere
    # (In case of a single-hemisphere atlas, only one GIFTI file is needed)
    files = ['myfile_hemi-L.func.gii','myfile_hemi-R.func.gii']
    X = atlas.read_data(files)
    # Read data from a single CIFTI file that combines left and right hemisphere
    X = atlas.read_data('my_single_cifti_file.dscalar.nii')


Writing atlas data to NIFTI, GIFTI, or CIFTI files
--------------------------------------------------
The functions ``data_to_nifti`` and ``data_to_cifti`` can be used to write data to a file in the atlas space.
A typical use case is to read data from a `NIFTI` file, does some computations, and writes it to a `CIFTI` file could look like this:

.. code-block:: python

    # Write data to a CIFTI file:
    cifti = atlas.data_to_cifti(Y)
    nb.save(cifti,'my_cifti_file.dscalar.nii')

    # Write data to a NIFTI file:
    nifti = atlas.data_to_nifti(Y)
    nb.save(nifti,'my_nifti_file.nii')

Using parcellations
-------------------
If you have a parcellation of your atlas, you can use a dseg.nii or label.gii ROI files to read into your Atlas, and then summarize your extracted data within those Parcels.
Here an example how to define parcels for a cortical atlas, and to average some data within those parcels:

.. code-block:: python

    # Get the atlas
    atlas, ainf = am.get_atlas('fs32k')
    # Get the label (1-K) for each vertex. 0 means not assigned
    gii_files = [atlas_dir + '/tpl-fs32k/Icosahedron1442.L.label.gii',
                 atlas_dir + '/tpl-fs32k/Icosahedron1442.R.label.gii']
    label_vec,labels = atlas.get_parcel(gii_files)
    # Average the data (ignoring Nans) in each ROI
    Yn = ds.agg_parcels(Y,label_vec,fcn=np.nanmean)

Saving parcellation results as pscalar CIFTI files
--------------------------------------------------
`CIFTI` files are very handy, in that they cannot only store volume and surface data, but also the data for parcels defined in the volume or on the surface. The connectome workbench displays these files correctly, without having to project them back into the full space.

And here is a full example for an ROI-analysis for a volumetric (cerebellar) atlas:

.. code-block:: python

    # if not already done, import nibabel
    import nibabel as nb

    # get the atlas
    atlas, ainf = am.get_atlas('SUIT2', atlas_dir)

    # Load the ROI file and define labels
    roi_files = atlas_dir + '/tpl-SUIT/atl-Anatom_space-SUIT_dseg.nii'
    label_vec, labels = atlas.get_parcel(roi_files)

    # Average some  data within each ROI
    Yn = ds.agg_parcels(Y, label_vec, fcn=np.nanmean)

    # create parcel axis for the cerebellum (will be used as column axis in pscalar file)
    p_axis = atlas.get_parcel_axis()

    # generate row axis with
    row_axis = nb.cifti2.ScalarAxis(row_labels)

    # Make the CIFTI file and save
    header = nb.Cifti2Header.from_axes((row_axis, p_axis))
    cifti_img = nb.Cifti2Image(Y, header=header)
    nb.save(cift_img,'myROIresult.pscalar.nii')

Transforming data between atlas spaces
-------------------------------------
Data extracted in one Atlas space can be directly transformed into another atlas space. For this, the two Atlasses need to cover the same brain structure. Currently, direct deformation is only possible between two volumetric atlases. The deformation depends on the `xfm` file found in the template directory of the target space. If the file does not exist, raise an issue on Github.

.. code-block:: python

    atlas_src,_ = am.get_atlas('SUIT3')
    atlas_trg,_ = am.get_atlas('MNISymC2')
    data = atlas_src.read_data('file_in_suit_space.dscalar.nii')
    data_new = am.deform_data(data,atlas_src, atlas_trg, interpolation=1)
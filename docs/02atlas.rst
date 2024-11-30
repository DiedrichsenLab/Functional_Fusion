Atlases and Parcellations
#########################

Atlas Spaces
------------

The Atlas class defines the atlas space, the voxel resolution, and the set of voxels or vertices that will form the basis of the analysis.
Subclasses are ``AtlasVolumetric``, and ``AtlasSurface``. Each Atlas contains  ``P`` locations (vertices or voxels) that are being sampled for the analysis. Atlases are indicated by short strings that indicate the Atlas. The following Atlases are defined so far:

* ``SUIT3``:  Cerebellum in SUIT space (3mm voxels)
* ``SUIT2``:  Cerebellum in SUIT space (2mm voxels)
* ``SUIT1``:  Cerebellum in SUIT space (1mm voxels)
* ``MNISymC3``: Cerebellum in MNI152NLin2009cSym space (3mm voxels)
* ``MNISymC2``: Cerebellum in MNI152NLin2009cSym space (2mm voxels)
* ``fs32k``: Left and Right hemisphere, using identical Medial wall mask
* ``fs32kAsym``: Left and Right hemisphere, using asymmetric medial wall mask  from HCP.
* ``fs32kAsym``: Left and Right hemisphere, using asymmetric medial wall mask  from HCP.
* ``fs32kAsym``: Left and Right hemisphere, using asymmetric medial wall mask  from HCP.
* ``MNIAsymBg1``: Basal Ganglia in MNI152NLin6Asym space (1mm voxels)
* ``MNIAsymBg2``: Basal Ganglia in MNI152NLin6Asym space (2mm voxels)

These atlases are defined in ``Functional_Fusion/Atlases/atlas_description.json``.
To add a new atlas, you need to add a entry into that json file. The Atlas needs the following information:

* ``type``: Class-name for the atlas object ("AtlasVolumetric", "AtlasSurface", ...)
* ``dir``: Directory where the atlas files are stored: Functional_Fusion/Atlases/<dir>)
* ``mask``: Name(s) of mask niftis/gifti files that define the atlas space.
* ``space``: Space in which coordinates are defined (e.g. "MNI152NLin6AsymC")
* ``structure``: cifti-structure name - needs to be a valid, predefined cifti structure (e.g. "cerebellum")

Optional fields:

* ``res``: resolution (not used yet) -voxel locations are defined by mask image
* ``normspace``: Norm space for suit.vol2surf -

Getting atlases
---------------
You can get an atlas by calling:

.. code-block:: python

    atlas,ainf = am.get_atlas(atlas_str)

``atlas_str`` is a string that defines the atlas (see above). The function returns a tuple, where ``atlas`` is an Atlas object, and ``ainf`` is a dictionary with information about the atlas.


Reading data from Nifti, Gifti, or Cifti files
----------------------------------------------

You can use the function ``read_data`` to get the data from a file in the atlas space.
Here an example for a volumetric atlas:

.. code-block:: python

    # Get the atlas
    atlas,ainf = am.get_atlas('SUIT3')
    # Read data from Nifti file using linear interpolation (1)
    X = atlas.read_data('my_nift_file',interpolation=1)

And here one for a surface atlas:

.. code-block:: python

    # Get the atlas
    atlas,ainf = am.get_atlas('fs32k')
    # Read data from two gifti files for left and right hemisphere
    files = ['myfile_hemi-L.func.gii','myfile_hemi-R.func.gii']
    X = atlas.read_data(files)
    # Read data from a single cifti file that combines left and right hemisphere
    X = atlas.read_data('my_single_cifti_file.dscalar.nii'')


Writing atlas data to Nifti, Gifti, or Cifti files
--------------------------------------------------
The functions ``data_to_nifti`` and ``data_to_cifti`` can be used to write data to a file in the atlas space.
A typical use case is to read data from a Nifti file, does some computations, and writes it to a Cifti file could look like this:

.. code-block:: python

    cifti = atlas.data_to_cifti(Y)
    nb.save(cifti,'my_cifti_file.dscalar.nii')
    # Write data to Nifti file:
    nifti = atlas.data_to_nifti(Y)
    nb.save(nifti,'my_nifti_file.nii')

Using parcellations
-------------------
You can also ise dseg.nii or label.gii ROI files to define a parcelation to do your computation on
an ROI-level
Here an example how to define ROIs for the cortical atlas, and to average some data within those ROIs:

.. code-block:: python

    # Get the atlas
    atlas,ainf = am.get_atlas('fs32k')
    # Get the label (1-K) for each vertex. 0 means not assigned
    gii_files = [atlas_dir + '/tpl-fs32k/Icosahedron1442.L.label.gii',
                atlas_dir + '/tpl-fs32k/Icosahedron1442.R.label.gii']
    label_vec,labels = atlas.get_parcel(gii_files)
    # Average the data (ignoring Nans) in each ROI
    Yn = ds.agg_parcels(Y,label_vec,fcn=np.nanmean)

Saving parcellation results as pscalar cifti files
--------------------------------------------------
Cifti files are very handy, in that they cannot only store volume and surface data, but also the data for ROIs defined in the volume or on the surface. The connectome workbench displays these files correctly, without having to project them back into the full space.

And here a full example for an ROI-analysis for a volumetric (cerebellar) atlas:

.. code-block:: python

    # get the atlas
    atlas,ainf = am.get_atlas('SUIT2',atlas_dir)
    # Load the ROI file and define labels
    roi_files = atlas_dir + '/tpl-SUIT/atl-Anatom_space-SUIT_dseg.nii'
    label_vec,labels = atlas.get_parcel(roi_files)
    # Average some  data within each ROI
    Yn = ds.agg_parcels(Y,label_vec,fcn=np.nanmean)
    # create parcel axis for the cerebellum (will be used as column axis in pscalar file)
    p_axis = atlas.get_parcel_axis()
    # generate row axis with
    row_axis = nb.cifti2.ScalarAxis(row_labels)
    # Make the cifti file and save
    header = nb.Cifti2Header.from_axes((row_axis, p_axis))
    cifti_img = nb.Cifti2Image(Y, header=header)
    nb.save(cift_img,'myROIresult.pscalar.nii')

Transforming data between atlas space
-------------------------------------
Data extracted in one Atlas space can be directly transformed into another atlas space. For this, the two Atlasses need to cover the same brain structure. Currently, direct deformation is only possible between two volumetric atlases. The deformation depends on the `xfm` file found in the template directory of the target space. If the file does not exist, raise an issue on Github.

.. code-block:: python

    atlas_src,_ = am.get_atlas('SUIT3')
    atlas_trg,_ = am.get_atlas('MNISymC2')
    data = atlas_src.read_data('file_in_suit_space.dscalar.nii')
    data_new = am.deform_data(data,atlas_src,atlas_trg,interpolation=1)
Data Analysis
#############

To use Functional fusion framework you need access to a folder that holds the datasets and atlases.

.. code-block:: python

    import Functional_Fusion.dataset as ds
    import Functional_Fusion.atlas_map as am
    import nibabel as nb
    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
    atlas_dir = base_dir + '/Atlases'


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


About Atlas Spaces
------------------

The Atlas class defines the atlas space, the voxel resolution, and the set of voxels or vertices that will form the basis of the analysis.
Subclasses are ``AtlasVolumetric``, and ``AtlasSurface``. Each Atlas contains  ``P`` locations (vertices or voxels) that are being sampled for the analysis. Atlases are indicated by short strings that indicate the Atlas. The following Atlases are defined so far:

* ``SUIT3``:  Cerebellum in SUIT space (3mm voxels)
* ``SUIT2``:  Cerebellum in SUIT space (2mm voxels)
* ``SUIT1``:  Cerebellum in SUIT space (1mm voxels)
* ``MNISymC3``: Cerebellum in MNI152NLin2009cSym space (3mm voxels)
* ``MNISymC2``: Cerebellum in MNI152NLin2009cSym space (2mm voxels)
* ``fs32k``: Left and Right hemisphere, using identical Medial wall mask
* ``fs32k_Asym``: Left and Right hemisphere, using asymmetric medial wall mask from HCP.

Getting atlases
---------------
You can get an atlas by calling:

.. code-block:: python

    atlas,ainf = am.get_atlas(atlas_str,atlas_dir)

``atlas`` is an Atlas object, and ``ainf`` is a dictionary with information about the atlas.


Reading and writing from atlas space to Nifti, Gifti, or Cifti files
--------------------------------------------------------------------

You can use the function ``read_data`` to get the data from a file in the atlas space.
The functions ``data_to_nifti`` and ``data_to_cifti`` can be used to write data to a file in the atlas space.
A typical use case is to read data from a Nifti file, does some computations, and writes it to a Cifti file could look like this:

.. code-block:: python

    # Get the atlas
    atlas,ainf = am.get_atlas('SUIT3',atlas_dir)
    # Read data from Nifti file using linear interpolation (1)
    X = atlas.read_data('my_nift_file',interpolation=1)
    # Do some computations
    Y = mycomputation(X)
    # Write data to Cifti file:
    cifti = atlas.data_to_cifti(Y)
    nb.save(cifti,'my_cifti_file.nii')
    # Write data to Nifti file:
    nifti = atlas.data_to_nifti(Y)
    nb.save(nifti,'my_nifti_file.nii')


Using parcellations
-------------------
You can also ise dseg.nii or label.gii ROI files to define a parcelation to do your computation on
an ROI-level
Here an example how to define ROIs for the cortical atlas, and to average some data within those ROIs:

.. code-block:: python

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


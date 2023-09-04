Data Analysis
#############

To use Functional fusion framework you need access to a folder that holds the datasets and atlases. 

.. code-block:: python

    import Functional_Fusion.dataset as ds
    import Functional_Fusion.atlas_map as am
    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
    atlas_dir = base_dir + '/Atlases'
    

Loading Data
------------
Loading the data to get a ``n_subj x n_cond x n_voxels`` tensor :

.. code-block:: python

    X,info,dataset_obj = ds.get_dataset(base_dir,
            dataset='MDTB',
            atlas='fs32k',
            sess='all',
            type='CondAll')





About Atlas Spaces
------------------

The Atlas class defines the group-atlas that the data is analyzed in. Subclasses are ``AtlasVolumetric``, and ``AtlasSurface``. The altases is in a specific space, but also has a mask that defines the ``P`` locations (vertices or voxels) that are being sampled for the analysis. Atlases are indicated by short strings that indicate the Atlas. d altlasses so far:

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

    atlas_map.get_atlas(string).

Using parcellations: 
--------------------

.. code-block:: python

    gii_files = [atlas_dir + '/tpl-fs32k/Icosahedron1442.L.label.gii',
                atlas_dir + '/tpl-fs32k/Icosahedron1442.R.label.gii']
    label_vec,labels = atlas.get_parcel(gii_files)
    Yn = ds.agg_parcels(Y,label_vec)



Transforming calculations into Nifti/Cifti/Gifti files
------------------------------------------------------
NA 

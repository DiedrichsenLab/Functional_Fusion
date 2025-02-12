Atlasmaps
=========

An ``Atlasmap`` defines the mapping between an ``Atlas`` (or region) in group space and the native space of each subject (or any other space the original data is stored in). For volume-based atlases, this mapping is usually through a (series of) deformation field(s), and a specific interpolation mode (or smoothing). For surface-based atlases, it is defined by a white and pial surface of the individual subject and a mapping rule (which cortical depth to use).
For fully-integrated datasets, the definition and application of the ``Atlasmap`` taken care of automatically in ``dataset.extract_all``. However, for custom ROI-analysis, it is sometimes useful to use an atlasmap directly to extract data for a specific ROI and subject.

Example of volume-based ROI analysis
-------------------------------------

.. code-block:: python

    # Import atlas map
    import Functional_Fusion.atlas_map as am

    # Define the the region, get only left hemisphere
    atlas,_ = am.get_atlas('MNI152NLin6AsymC')
    # Set the Gifti file for the region (func.gii or label.gii)
    subatlas = atlas_left.get_subatlas_image('Path_to_roi_img.nii')


The subatlas will now have the $P$ locations in voxel space. You can use the subatlas.data_to_nifti() function to save data in that group space. For mapping data into the group space, we need to define an AtlasMap.

.. code-block:: python

    #  Define atlas map
    deform = '/sub-01_to_atlasspace_xfm.nii'                 # Deformation file
    mask = glm_dir + '/sub-01/mask.nii'                      # Mask in functional space
    amap = am.AtlasMapDeform(subatlas.voxels,deform,mask) # Atlas map
    amap.build(interpolation=1)  # Using Trilinear interpolation (0 for nearest neighbor, 2 for smoothing)


You can the proceed with data extract as shown below.


Example of surface-based ROI analysis
-------------------------------------

.. code-block:: python

    # Import atlas map
    import Functional_Fusion.atlas_map as am
    import nitools as nt

    # Define the the region, get only left hemisphere
    atlas,_ = am.get_atlas('fs32k')
    atlas_left = atlas.get_hemisphere(0)
    # Equivalently you could have used
    atlas_left,_ = am.get_atlas('fs32k_L')

    # Set the Gifti file for the region (func.gii or label.gii)
    subatlas = atlas_left.get_subatlas_image('Path_to_roi_img.gii')

The subatlas will now have the $P$ locations in vertex space. You can use the subatlas.data_to_cifti() function to save data in that group space. For mapping data into the group space, we need to define an AtlasMap.

.. code-block:: python

    #  Define atlas map
    white = surf_dir + '/sub-01/sub-01.L.white.32k.surf.gii' # Individual white surface
    pial = surf_dir + '/sub-01/sub-01.L.pial.32k.surf.gii'   # Invividual pial surface
    mask = glm_dir + '/sub-01/mask.nii'                      # Mask in functional space
    amap = am.AtlasMapSurf(subatlas.vertex[0],white,pial,mask) # Atlas map
    amap.build()

Data Extraction using atlas maps
--------------------------------

Once the Atlas map is built (surface or volume based), you can use it to extract data from the native space of the subject.

* The function ``extract_data_native()`` will extract the data from all the voxel in native space of the subject that map to group space.
* The function ``extract_data_group()`` will extract the data in group space.
* The function ``map_native_to_group()`` will map the data from native to group space.

.. code-block:: python

    dnames = ['beta_0001.nii','beta_0002.nii','beta_0003.nii'] # Data files can be 3d- or 4d-niftis
    n_data = amap.extract_data_native(dnames)
    g_data = amap.map_native_to_group(n_data) 
    g_data = amap.extract_data_group(dnames) # Results in the same as the above two lines

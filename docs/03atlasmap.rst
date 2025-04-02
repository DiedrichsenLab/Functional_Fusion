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
    # Get a Nifti-file of the group atlas. If it is a 0/1 ROI image:
    subatlas = atlas.get_subatlas_image('Path_to_roi_img.nii')


If you have a discrete segmentation volumetric atlas with different ROIs in it, you can pick out specific values from the file: 

.. code-block:: python

    # Get a the areas 18 from the dseg file: 
    subatlas = atlas_left.get_subatlas_image('myatlas_dseg.nii',label_value=18)
    # Or if you multiple possible values 
    subatlas = atlas_left.get_subatlas_image('myatlas_dseg.nii',label_value=[2,18])



The subatlas will now have the ``P`` locations in voxel space. You can use the ``subatlas.data_to_nifti()`` function to save data in that group space. For mapping data into the group space, we need to define an ``AtlasMapDeform``.

.. code-block:: python

    #  Define atlas map
    deform = '/sub-01_to_atlasspace_xfm.nii'                 # Deformation file
    mask = glm_dir + '/sub-01/mask.nii'                      # Mask in functional space
    amap = am.AtlasMapDeform(atlas.world,deform,mask) # Atlas map
    amap.build(interpolation=1)  # Using Trilinear interpolation (0 for nearest neighbor, 2 for smoothing)
    # save the ROI mask in native space
    amap.save_as_image('/sub-01/ROI_mask.nii') 


You can the proceed with data extract as shown below.


Example of surface-based ROI analysis
-------------------------------------

The first step to define a surface-based ROI is to get the atlas for the corresponding hemisphere for the group surface atlas. 

.. code-block:: python

    # Import atlas map
    import Functional_Fusion.atlas_map as am
    import nitools as nt

    # Define the the region, get only left hemisphere
    atlas,_ = am.get_atlas('fs32k')
    atlas_left = atlas.get_hemisphere(0)
    # Equivalently you could have used
    atlas_left,_ = am.get_atlas('fs32k_L')

A surface-based ROI is usually defined in a gifti- or cifti-file that indicates whether the surface node is part of the ROI or not (0/1). Sometime we have discrete parcellation files (*.label.gii or _dseg.nii) that indicates multiple ROIS with integer numbers. As for the volume-based ROI you can also specify a label_value to pick out a specific (set of) ROIs from a discrete segementation atlas. 

.. code-block:: python

    # Set the Gifti file for the region (func.gii or label.gii)
    # This one uses any value >0 as part of the ROI
    subatlas = atlas_left.get_subatlas_image('Path_to_roi_img.gii')
    # Here an example of using one specific value
    subatlas = atlas_left.get_subatlas_image('Path_to_roi_img.gii', label_value=18)


The subatlas will now have the ``P`` locations in vertex group space. You can use the ``subatlas.data_to_cifti()`` function to save data in that group space. 

For mapping data between group space and individual space, we need to define an ``AtlasMapSurf``. This is done over the individual pial and whilte surface. 

.. code-block:: python

    #  Define atlas map
    white = surf_dir + '/sub-01/sub-01.L.white.32k.surf.gii' # Individual white surface
    pial = surf_dir + '/sub-01/sub-01.L.pial.32k.surf.gii'   # Invividual pial surface
    mask = glm_dir + '/sub-01/mask.nii'                      # Mask in functional space for that subject
    amap = am.AtlasMapSurf(subatlas.vertex[0],white,pial,mask) # Atlas map
    # Compute the voxels in native space 
    amap.build()
    # save the ROI mask in native space for checking only 
    amap.save_as_image('/sub-01/ROI_mask.nii') 

Data Extraction using atlas maps
--------------------------------

Once the Atlas map is built (surface or volume-based), you can use it to extract data from the native space of the subject.

* The function ``extract_data_native()`` will extract the data from all the voxel in native space of the subject that map to group space.
* The function ``extract_data_group()`` will extract the data in group space.
* The function ``map_native_to_group()`` will map the data from native to group space.
* The function ``save_as_image()`` saves the ROI as a 1/0 mask in native space.

.. code-block:: python

    dnames = ['beta_0001.nii','beta_0002.nii','beta_0003.nii'] # Data files that you want to map can be 3d- or 4d-niftis
        
    # This extract all the relevant voxels in native space (use for RSA)
    n_data = amap.extract_data_native(dnames)

    # this statement maps the data to group space 
    g_data = amap.extract_data_group(dnames)

    # Actually, the mapping to group space consists of the following 2 lines of code: 
    n_data = amap.extract_data_native(dnames)
    # This maps native data to group space 
    g_data = amap.map_native_to_group(n_data) 

The advantage of ussing map_native_to_group is that you can do some computation on data in native space and then map and save it in group space for subsequent analysis. 
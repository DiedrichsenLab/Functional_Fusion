Data Analysis
#############

Data
====

Loading Data
------------



Atlases
=======
About Atlas Spaces
------------------

The Atlas class defines the group-atlas that the data is analyzed in. Subclasses are `AtlasVolumetric`, and `AtlasSurface`. The altases is in a specific space, but also has a mask that defines the `P` locations (vertices or voxels) that are being sampled for the analysis. Atlases are indicated by short strings that indicate the Atlas. d altlasses so far:

* `SUIT3`:  Cerebellum in SUIT space (3mm voxels)
* `SUIT2`:  Cerebellum in SUIT space (2mm voxels)
* `SUIT1`:  Cerebellum in SUIT space (1mm voxels)
* `MNISymC3`: Cerebellum in MNI152NLin2009cSym space (3mm voxels)
* `MNISymC2`: Cerebellum in MNI152NLin2009cSym space (2mm voxels)
* `fs32k`: Left and Right hemisphere, using identical Medial wall mask 
* `fs32k_Asym`: Left and Right hemisphere, using asymmetric medial wall mask from HCP. 

Getting atlases
---------------
You can get an atlas by calling `atlas_map.get_atlas(string)`. Define

Transforming calculations into Nifti/Cifti/Gifti files
------------------------------------------------------


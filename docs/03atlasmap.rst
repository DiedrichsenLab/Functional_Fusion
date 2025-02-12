Atlasmaps
=========

An ``Atlasmap`` defines the mapping between an ``Atlas`` (or region) in group space and the native space of each subject (or any other space the original data is stored in). For volume-based atlases, this mapping is usually through a (series of) deformation field(s), and a specific interpolation mode (or smoothing). For surface-based atlases, it is defined by a white and pial surface of the individual subject and a mapping rule (which cortical depth to use).
For fully-integrated datasets, the definition and application of the ``Atlasmap`` taken care of automatically in ``dataset.extract_all``. However, for custom ROI-analysis, it is sometimes useful to use an atlasmap directly to extract data for a specific ROI and subject.

Here a full example of the definition and data extract for a volume-based ROI:







Here a full example of the definition and data extract for a surface-based ROI:

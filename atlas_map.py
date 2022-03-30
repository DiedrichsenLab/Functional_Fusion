#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The functions of atlas mapping

Created on 3/30/2022 at 3:00 PM
Author: dzhi
"""
import numpy as np
import nilearn as nl
import os

PATH = 'mapping'


def atlas_map(atlas, resolution):
    """The function to return the atlas map as required
    Args:
        atlas: which atlas is being used
        resolution: the resolution that wants to map

    Returns: the atlas mapper
    """
    mapper = os.path.join(PATH, "file_%s_%d.txt" % (atlas, resolution))
    return mapper

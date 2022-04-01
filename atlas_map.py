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


class atlas_map:
    def __int__(self, name, sub):
        self.name = name
        self.sub = sub

    def build(self, xyz, nii_file):
        """
        Args:
            xyz:
            nii_file:

        Returns: xyz in subject space

        """
        pass

    def save(self, file_name):
        pass

    def load(self, file_name):
        pass


# def atlas_map(subject, atlas, resolution):
#     """The function to return the atlas map as required
#     Args:
#         subject: give the subject
#         atlas: which atlas is being used
#         resolution: the resolution that wants to map
#
#     Returns: the atlas mapper
#     """
#     mapper = os.path.join(PATH, "file_%s_%d.txt" % (atlas, resolution))
#     return mapper

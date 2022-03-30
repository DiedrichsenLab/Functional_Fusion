#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data fusion project dataset module

The class for converting and mapping raw data from multi-dataset
to a standard data structure that can be used in Diedrichsen lab

Created on 3/30/2022 at 12:21 PM
Author: dzhi
"""
import numpy as np
import torch as pt
import nilearn as nl


class dataset:
    def __init__(self, atlas, raw_data, target='betas'):
        # Constructor of dataset class, and the class
        # knows which atlas will be used for this individual's
        # raw data when initializing.
        self.atlas = atlas
        self.data = raw_data
        self._data_processing(target=target)

    def _minimal_preprocessing(self):
        # TODO: the minimal preprocessing pipeline
        #  to process the raw input data
        pass

    def glm(self):
        pass

    def connectivity(self):
        pass

    def _data_processing(self, target='betas', return_data=False):
        """The minimum preprocessing pipeline to process the
           raw individual data
        Args:
            target: 'betas' - to get beta estimates
                    'conn' - to get the connectivity fingerprint

        Returns:
        """
        self._minimal_preprocessing()
        if target == 'betas':
            self.glm()
        elif target == 'conn':
            self.connectivity()
        else:
            raise ValueError('the value of target is invalid.')

        if return_data:
            return self.data

    def get_data(self, atlas_map):
        """the main function to output the processed data
        Args:
            atlas_map: atlas mapper to find the voxels
                       associated to the vertices

        Returns: Y - the processed individual data
        """
        pass

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
import pandas as pd


class DataSet:
    def __init__(self, base_dir):
        """DataSet class:
        Implements the interface for each of the data set
        Note that the actual preprocessing and glm estimate
        do not have to be performed with functionality provided by
        this class. The class is just a instrument to present the user with
        a uniform interface of how to get subject info

        Args:
            basedir (str): _description_
        """
        self.base_dir  = base_dir
        self.surface_dir = base_dir + '/{0}/surface'
        self.anatomical_dir = base_dir + '/{0}/anatomical'
        self.contrast_dir = base_dir + '/{0}/contrast'
        self.suit_dir = base_dir + '/{0}/suit'

    def get_participants(self):
        """ returns a data frame with all participants
        available in the study. The fields in the data frame correspond to the
        standard columns in participant.tsv.
        https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html
        Returns:
            Pinfo (pandas data frame): participant information in standard bids format
        """
        Pinfo = pd.read_csv(self.base_dir + '/participant.tsv')
        return Pinfo

    def get_data(self, participant_id, atlas_map):
        """the main function to output the processed data
        Args:
            participant_id: standard participant_id
            atlas_map: AtlasMAP to find the voxels

        Returns:
            Y (np.ndarray):
                A N x P numpy array of aggregated data
            T (pd.DataFrame):
                A data frame with information about the N numbers provide
        """
        pass

# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
import Functional_Fusion.dataset as ds
import nibabel as nb
from matrix import indicator
import sys
import os

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'


if __name__ == "__main__":
    # make_mdtb_suit()
    # test_atlas_sym()
    data,info,ds = ds.get_dataset(base_dir,'Demand',atlas='MNISymC3')
    pass


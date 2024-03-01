# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
from Functional_Fusion.matrix import indicator
import nibabel as nb

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'


def test_extract_mdtb(): 
    dataset = ds.DataSetMDTB(base_dir + '/MDTB')
    dataset.extract_all(ses_id='ses-s1',
                        type='CondAll',
                        atlas='SUIT3',
                        smooth=2,
                        subj=[0])
    pass


if __name__ == "__main__":
    test_extract_mdtb()
 
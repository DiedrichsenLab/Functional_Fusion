# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
from Functional_Fusion.matrix import indicator
import nibabel as nb

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'


def test_extract_mdtb():
    dataset = ds.DataSetMDTB(base_dir + '/MDTB')
    dataset.extract_all(ses_id='ses-s2',
                        type='CondAll',
                        atlas='MNISymC3',
                        smooth=2,
                        subj=[21,22,23])
    pass

def test_extract_dmcc():
    dataset = ds.DataSetDmcc(base_dir + '/DMCC')
    dataset.extract_all(ses_id='ses-axcpt_bas',
                        type='CondAll',
                        atlas='MNIAsymBg2',
                        smooth=1,
                        subj=[0])
    pass

def test_extract_somatotopic():
    dataset = ds.DataSetSomatotopic(base_dir + '/Somatotopic')
    dataset.extract_all(ses_id='ses-motor',
                        type='CondAll',
                        atlas='MNIAsymBg2',
                        smooth=1,
                        subj=[0])
    pass


if __name__ == "__main__":
    test_extract_mdtb()
    # test_extract_dmcc()
    # test_extract_somatotopic()

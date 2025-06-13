"""
Script for testing the extraction of data from the MDTB dataset with a universal condense_data in the dataset class.
"""

from pathlib import Path
from Functional_Fusion.dataset import DataSetMDTB
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as util

base_dir = util.get_base_dir()


def test_extract(dataset,sess,space,type):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    mydataset.extract_all(ses_id=sess, type=type, atlas=space)



if __name__ == "__main__":
    # social
    # test_extract('social', 'ses-social', 'fs32k', 'CondRun')
    # test_extract('social', 'ses-social', 'MNISymC3', 'CondRun')
    # test_extract('social', 'ses-social', 'fs32k', 'CondHalf')
    # test_extract('social', 'ses-social', 'MNISymC3', 'CondHalf')
    # test_extract('social', 'ses-social', 'fs32k', 'CondAll')
    # test_extract('social', 'ses-social', 'MNISymC3', 'CondAll')

    # language
    test_extract('language', 'ses-localizer', 'fs32k', 'CondRun')
    test_extract('language', 'ses-localizer', 'MNISymC3', 'CondRun')
    test_extract('language', 'ses-localizer', 'fs32k', 'CondHalf')
    test_extract('language', 'ses-localizer', 'MNISymC3', 'CondHalf')
    test_extract('language', 'ses-localizer', 'fs32k', 'CondAll')
    test_extract('language', 'ses-localizer', 'MNISymC3', 'CondAll')
    test_extract('language', 'ses-localizerfm', 'fs32k', 'CondRun')
    test_extract('language', 'ses-localizerfm', 'MNISymC3', 'CondRun')
    test_extract('language', 'ses-localizerfm', 'fs32k', 'CondHalf')
    test_extract('language', 'ses-localizerfm', 'MNISymC3', 'CondHalf')
    test_extract('language', 'ses-localizerfm', 'fs32k', 'CondAll')
    test_extract('language', 'ses-localizerfm', 'MNISymC3', 'CondAll')

    # language but TaskAll and taskHalf and TaskRun
    test_extract('language', 'ses-localizer', 'fs32k', 'TaskAll')
    test_extract('language', 'ses-localizer', 'MNISymC3', 'TaskAll')
    test_extract('language', 'ses-localizerfm', 'fs32k', 'TaskAll')
    test_extract('language', 'ses-localizerfm', 'MNISymC3', 'TaskAll')
    test_extract('language', 'ses-localizer', 'fs32k', 'TaskHalf')
    test_extract('language', 'ses-localizer', 'MNISymC3', 'TaskHalf')
    test_extract('language', 'ses-localizerfm', 'fs32k', 'TaskHalf')
    test_extract('language', 'ses-localizerfm', 'MNISymC3', 'TaskHalf')
    test_extract('language', 'ses-localizer', 'fs32k', 'TaskRun')
    test_extract('language', 'ses-localizer', 'MNISymC3', 'TaskRun')
    test_extract('language', 'ses-localizerfm', 'fs32k', 'TaskRun')
    test_extract('language', 'ses-localizerfm', 'MNISymC3', 'TaskRun')





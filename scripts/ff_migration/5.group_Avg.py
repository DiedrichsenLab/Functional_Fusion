"""
Script for making group avgs
"""

from pathlib import Path
from Functional_Fusion.dataset import DataSetMDTB
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as util

base_dir = util.get_base_dir()

def make_avg(dataset,sess,space,type):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    mydataset.group_average_data( ses_id=sess,
                               type=type,
                               atlas=space,
                               subj=None)
    

if __name__ == "__main__": 

     # IBC with all its session
    # sessions = ['ses-archi',
    #                      'ses-clips4',
    #                      'ses-enumeration',
    #                      'ses-hcp1', 'ses-hcp2',
    #                      'ses-lyon1', 'ses-lyon2',
    #                      'ses-mathlang',
    #                      'ses-mtt1', 'ses-mtt2',
    #                      'ses-preference',
    #                      'ses-rsvplanguage',
    #                      'ses-spatialnavigation',
    #                      'ses-tom']
    # data_types = ['CondHalf','CondRun','CondAll']
    # spaces = ['fs32k','MNISymC3']

    # for session in sessions:
    #     for data_type in data_types:
    #         for space in spaces:
    #             make_avg('IBC', session, space, data_type)


    # social
    make_avg('social', 'ses-social', 'fs32k', 'CondRun')
    make_avg('social', 'ses-social', 'MNISymC3', 'CondRun')
    make_avg('social', 'ses-social', 'fs32k', 'CondHalf')
    make_avg('social', 'ses-social', 'MNISymC3', 'CondHalf')
    make_avg('social', 'ses-social', 'fs32k', 'CondAll')
    make_avg('social', 'ses-social', 'MNISymC3', 'CondAll')

    # language
    make_avg('language', 'ses-localizer', 'fs32k', 'CondRun')
    make_avg('language', 'ses-localizer', 'MNISymC3', 'CondRun')
    make_avg('language', 'ses-localizer', 'fs32k', 'CondHalf')
    make_avg('language', 'ses-localizer', 'MNISymC3', 'CondHalf')
    make_avg('language', 'ses-localizer', 'fs32k', 'CondAll')
    make_avg('language', 'ses-localizer', 'MNISymC3', 'CondAll')
    make_avg('language', 'ses-localizerfm', 'fs32k', 'CondRun')
    make_avg('language', 'ses-localizerfm', 'MNISymC3', 'CondRun')
    make_avg('language', 'ses-localizerfm', 'fs32k', 'CondHalf')
    make_avg('language', 'ses-localizerfm', 'MNISymC3', 'CondHalf')
    make_avg('language', 'ses-localizerfm', 'fs32k', 'CondAll')
    make_avg('language', 'ses-localizerfm', 'MNISymC3', 'CondAll')



    pass

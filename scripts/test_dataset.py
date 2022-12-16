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

def reliability_ibc():
    dataset = ds.DataSetIBC(base_dir + '/IBC')
    # Specify the fields you want to have / check 
    T = dataset.get_participants()
    dataset.sessions = dataset.sessions[9:12]
    num_sess = len(dataset.sessions)
    RW = np.empty((T.shape[0],num_sess))
    RB = np.empty((T.shape[0],num_sess))
    Missing = np.empty((T.shape[0],num_sess))

    for i,ses in enumerate(dataset.sessions):
        data,info = dataset.get_data('MNISymC3',ses,'CondHalf')
        m = np.isnan(data).sum(axis=1)
        Missing[:,i] = (m>0).sum(axis=1)
        rw = ds.reliability_within_subj(data,part_vec=info.half,cond_vec=info.reg_id)
        RW[:,i] = rw.mean(axis=1)
        RB[:,i] = ds.reliability_between_subj(data,cond_vec=info.reg_id)
    pass

if __name__ == "__main__":
    # make_mdtb_suit()
    reliability_ibc()
    # data,info,ds = ds.get_dataset(base_dir,'Demand',atlas='MNISymC3')
    pass


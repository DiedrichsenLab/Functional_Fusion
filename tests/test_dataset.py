# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import numpy as np
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as util
from Functional_Fusion.matrix import indicator
import nibabel as nb

base_dir = util.get_base_dir()

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

def test_get_data():
    dataset,info,myds= ds.get_dataset(base_dir,'MDTB',atlas='MNISymC3',subj=[0,1,2])
    dataset,info,myds= ds.get_dataset(base_dir,'Demand',atlas='MNISymC3',subj=[0,1,2])
    dataset,info,myds= ds.get_dataset(base_dir,'Pontine',atlas='MNISymC3',subj=[0,1,2])
    dataset,info,myds= ds.get_dataset(base_dir,'HCPur100',atlas='MNISymC3',subj=[0,1,2])
    pass

def test_extract(dataset,sess,space,type):
    mydataset = ds.get_dataset_class(base_dir,dataset)
    mydataset.extract_all(ses_id=sess, type=type, atlas=space)

 


def test_decompose(): 
    N = 5
    R = 3 
    C = 10
    P = 80 

    data = np.random.randn(N,R,C,P)
    variances = ds.decompose_pattern_into_group_indiv_noise(data,'global')
    variances1 = ds.decompose_pattern_into_group_indiv_noise(data,'voxel_wise')
    variances2 = ds.decompose_pattern_into_group_indiv_noise(data,'condition_wise')
    
    pass


if __name__ == "__main__":
    # make_mdtb_suit()
    # test_decompose()
    # test_get_data()
    # test_extract('Demand','ses-01','fs32k','CondHalf')
    # test_extract('Demand','ses-01','MNISymC3','CondAll')
    # test_extract('MDTB','ses-s1','fs32k','CondRun')
    # test_extract('MDTB','ses-s2','fs32k','CondRun')
    # test_extract('MDTB','ses-s2','MNISymC3','CondRun')
    # test_extract('MDTB','ses-s1','fs32k','CondHalf')
    # test_extract('WMFS','ses-01','fs32k','CondHalf')
    # test_extract('WMFS','ses-02','fs32k','CondHalf')
    # test_extract('WMFS','ses-01','fs32k','CondRun')
    # test_extract('WMFS','ses-02','fs32k','CondRun')
    test_extract('WMFS','ses-01','fs32k','CondAll')
    test_extract('WMFS','ses-02','fs32k','CondAll')
    test_extract('WMFS','ses-01','MNISymC3','CondHalf')
    test_extract('WMFS','ses-02','MNISymC3','CondHalf')
    # test_extract('WMFS','ses-01','MNISymC3','CondRun')
    # test_extract('WMFS','ses-02','MNISymC3','CondRun')
    test_extract('WMFS','ses-01','MNISymC3','CondAll')
    test_extract('WMFS','ses-02','MNISymC3','CondAll')

    # data,info,ds = ds.get_dataset(base_dir,'Demand',atlas='MNISymC3')
    pass


import numpy as np
import os 
import shutil
import pandas as pd
import nibabel as nb
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as ut
# compares versions of the new ald old functional fusion 
base_dir = '/Volumes/diedrichsen_data$/data'

def compare_data(dataset,sess, subj,space,type):
    name = f'{subj}_space-{space}_{sess}_{type}.dscalar.nii'
    if (dataset =='Language'): 
        oname = f'{subj}_space-{space}_{sess}_cond_{type}.dscalar.nii'
    else:
        oname = name
    old_path  = os.path.join(base_dir, 'FunctionalFusion', dataset, 'derivatives',subj,'data',oname)
    new_path = os.path.join(base_dir, 'FunctionalFusion_new', dataset, 'derivatives','ffextract',subj,name) 
    A= nb.load(old_path).get_fdata()
    B = nb.load(new_path).get_fdata()

    baseline_old = np.nanstd(np.nanmean(A,axis=0))
    baseline_new = np.nanstd(np.nanmean(B,axis=0))
    missing_old = np.sum(np.isnan(np.sum(A,axis=0)))
    missing_new = np.sum(np.isnan(np.sum(B,axis=0)))

    A = A.flatten()
    B = B.flatten()

    C = A+B
    R= np.corrcoef(A[~np.isnan(C)],B[~np.isnan(C)])
    return R[0,1],baseline_old,baseline_new,missing_old,missing_new

def compare_all(
    datasets = ['MDTB','Demand','Nishimoto','Somatotopic','IBC','Language','WMFS','Social'],
    spaces = ['fs32k','MNISymC3'],
    type = 'CondHalf'):
    D = pd.DataFrame()
    for i,dataset in enumerate(datasets):
        for space in spaces:
            myds = ds.get_dataset_class(ut.get_base_dir(),dataset)
            T=myds.get_participants()
            for s in range(len(T)):
                mysubj = T.participant_id.iloc[s]                   
                for mysess in myds.sessions:
                    if mysess != 'ses-rest':
                        R,bo,bn,mo,mn= compare_data(dataset,mysess,mysubj,space,type)
                        print(f'{dataset} {mysess} {mysubj} {space} {type}: {R:.3f}')
                        d = {'dataset': dataset, 'session': mysess, 'subject': mysubj,
                            'space': space,
                            'type': type,
                            'session' : mysess,
                            'participant_id': mysubj,
                            'correlation': R,
                            'baseline_old': bo,
                            'baseline_new': bn,
                            'missing_old': mo,
                            'missing_new': mn}
                        D = pd.concat([D, pd.DataFrame([d])], ignore_index=True)
    return D 

if __name__ == '__main__':
    D= compare_all(['Demand'],type='CondAll')
    pass 
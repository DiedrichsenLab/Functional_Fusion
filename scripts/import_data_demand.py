# Imports all the data from the pe / sigmasquares files from individual 
# runs. Concatenates them into dscalar cifti-files per run
from pathlib import Path
import numpy as np
import nibabel as nb
import nitools as nt
from Functional_Fusion.dataset import DataSetDemand 
import pandas as pd
import os

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/FunctionalFusion/Demand'
orig_dir = base_dir + '/Cerebellum/Demand'
phase_enc = ['AP','PA','AP','PA']
def import_data():
    dataset = DataSetDemand(data_dir)
    T = dataset.get_participants()
    for j,sid in enumerate(T.participant_id):
        dirw = dataset.estimates_dir.format(sid) + '/ses-01'
        if not os.path.exists(dirw):
           os.makedirs(dirw)
        info = pd.read_csv(data_dir+f'/ses-01_reginfo.tsv',sep='\t')
        # get the correct pe files and concatenate them into a dscalar
        data =[];
        old_sid = T.orig_id[j]
        for i,row in info.iterrows():
            print(i)
            rdir = f'tfMRI_visual_r{row.run}_{phase_enc[row.run-1]}_hp200_s2_level1_MSMAll_DeDrift_Resample_hp0_clean.feat'
            filename = orig_dir + f'/{old_sid}/{rdir}/GrayordinatesStats/pe{row.pe_id}.dtseries.nii'
            cifti=nb.load(filename)
            data.append(cifti.get_fdata())
        bm = cifti.header.get_axis(1)
        row_axis=nb.cifti2.ScalarAxis(info.cond_name)
        header = nb.Cifti2Header.from_axes((row_axis,bm))
        cifti_img = nb.Cifti2Image(dataobj=np.concatenate(data),header=header)
        outname = dirw + f'/{sid}_ses-01_beta.dscalar.nii'
        nb.save(cifti_img,outname)

        outname = dirw + f'/{sid}_ses-01_reginfo.tsv'
        info['sn']=[sid]*info.shape[0]
        info.to_csv(outname,sep='\t',index=False)

        # Now add the sigma-square images into a seperate file
        data = []
        names = []
        for r in [1,2,3,4]:
            rdir = f'tfMRI_visual_r{row.run}_{phase_enc[row.run-1]}_hp200_s2_level1_MSMAll_DeDrift_Resample_hp0_clean.feat'
            filename = orig_dir + f'/{old_sid}/{rdir}/GrayordinatesStats/sigmasquareds.dtseries.nii'
            cifti=nb.load(filename)
            data.append(cifti.get_fdata())
        bm = cifti.header.get_axis(1)
        row_axis=nb.cifti2.ScalarAxis(['mean_resms'])
        header = nb.Cifti2Header.from_axes((row_axis,bm))
        data = np.concatenate(data).mean(axis=0,keepdims=True)
        cifti_img = nb.Cifti2Image(dataobj=data,header=header)
        outname = dirw + f'/{sid}_ses-01_resms.dscalar.nii'
        nb.save(cifti_img,outname)


if __name__ == '__main__':
    import_data()
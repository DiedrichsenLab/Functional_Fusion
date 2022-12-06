# Imports all the data from the pe / sigmasquares files from individual 
# runs. Concatenates them into dscalar cifti-files per run
from pathlib import Path
import numpy as np
import nibabel as nb
import nitools as nt
from Functional_Fusion.dataset import DataSetDemand 
import pandas as pd

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/FunctionalFusion/Demand'
orig_dir = base_dir + '/Cerebellum/Demand/raw'

def import_data():
    dataset = DataSetDemand(data_dir)
    T = dataset.get_participants()
    for j,sid in enumerate(T.participant_id):
        dirw = dataset.estimates_dir.format(sid) + '/ses-01'
        info = pd.read_csv(dirw+f'/{sid}_ses-01_reginfo.tsv',sep='\t')
        # get the correct pe files and concatenate them into a dscalar
        data =[];
        old_sid = f'S{j+1:02d}'
        for i,row in info.iterrows():
            print(i)
            filename = orig_dir + f'/{old_sid}/level1_run{row.run}/GrayordinatesStats/pe{row.pe_id}.dtseries.nii'
            cifti=nb.load(filename)
            data.append(cifti.get_fdata())
        bm = cifti.header.get_axis(1)
        row_axis=nb.cifti2.ScalarAxis(info.cond_name)
        header = nb.Cifti2Header.from_axes((row_axis,bm))
        cifti_img = nb.Cifti2Image(dataobj=np.concatenate(data),header=header)
        outname = dirw + f'/{sid}_ses-01_pe.dscalar.nii'
        nb.save(cifti_img,outname)

        # Now add the sigma-square images into a seperate file
        data = []
        names = []
        for r in [1,2,3,4]:
            filename = orig_dir + f'/{old_sid}/level1_run{row.run}/GrayordinatesStats/sigmasquareds.dtseries.nii'
            cifti=nb.load(filename)
            data.append(cifti.get_fdata())
            names.append(f'run{r}')
        bm = cifti.header.get_axis(1)
        row_axis=nb.cifti2.ScalarAxis(names)
        header = nb.Cifti2Header.from_axes((row_axis,bm))
        cifti_img = nb.Cifti2Image(dataobj=np.concatenate(data),header=header)
        outname = dirw + f'/{sid}_ses-01_sigmasquared.dscalar.nii'
        nb.save(cifti_img,outname)


if __name__ == '__main__':
    import_data()
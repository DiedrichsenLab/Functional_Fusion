"""
Script for importing the Pontine dataset to general format.

Created Sep 2022
Author: caro nettekoven
"""

import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import scipy.io as sio
from import_data import *
from Functional_Fusion.dataset import DataSetSomatotopic
import Functional_Fusion.atlas_map as am
import shutil
import nibabel as nb

base_dir = '/Volumes/diedrichsen_data$/data'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data'


src_base_dir = Path(base_dir + '/Cerebellum//Somatotopic')
dest_base_dir = Path(base_dir + '/FunctionalFusion/Somatotopic')


def import_data():
    dataset = DataSetSomatotopic(str(dest_base_dir))
    info = pd.read_csv(dest_base_dir / 'ses-all_reginfo.tsv', sep='\t')
    T = dataset.get_participants()
    for _, id in T.iterrows():
        print(f'Importing {id.participant_id}')
        # for session in np.arange(4):
        for session in [3]:
            resms = []
            for run in np.arange(info.run.max()):
                # for pe in info.pe_id.unique():

                #     src = src_base_dir / \
                #         f'raw/Functional/{id.orig_id}/{id.orig_id}_sess{session+1:02d}_MOTOR{run+1}/pe{pe}.nii.gz'
                    

                #     reg=info[info.pe_id == pe].reg_id.iloc[0]
                #     dest = dest_base_dir / \
                #         f'derivatives/{id.participant_id}/estimates/ses-{session+1:02d}/{id.participant_id}_ses-{session+1:02d}_run-{run+1:02d}_reg-{reg:02d}_beta.nii'
                    
                #     # Make folder
                #     dest.parent.mkdir(parents=True, exist_ok=True)

                #     # Copy func file to destination folder and rename
                #     if  ~dest.exists():
                #         try:
                #             shutil.copyfile(src,
                #                     dest)
                #         except:
                #             print('skipping ' + str(src))


                
                src = src_base_dir / \
                    f'raw/Functional/{id.orig_id}/{id.orig_id}_sess{session+1:02d}_MOTOR{run+1}/sigmasquareds.nii.gz'
            #     resms_img = nb.load(src)
            #     resms.append(resms_img.get_fdata())

            # resms = np.mean(resms,axis=0)
            # nifti_img = nb.Nifti1Image(dataobj=resms, affine=resms_img.affine)
            outname = dest_base_dir / \
                f'derivatives/{id.participant_id}/estimates/ses-{session+1:02d}/{id.participant_id}_ses-{session+1:02d}_resms.nii'
            # nb.save(nifti_img, outname)
            
            # Import the reginfo.tsv file from the general file
            info.insert(loc=0, column='sn', value=[id.participant_id]*info.shape[0])
            info.to_csv(dest_base_dir /
                        f'derivatives/{id.participant_id}/estimates/ses-{session+1:02d}/{id.participant_id}_ses-{session+1:02d}_reginfo.tsv', sep='\t', index=False)

            pass


        



if __name__ == '__main__':
    
    # --- Importing Estimates ---
    import_data()



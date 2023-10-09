# Script for importing the IBC data set from super_cerebellum to general format.
import os
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetIBC
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import subprocess

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    raise NameError('Could not find base_dir')

data_dir = base_dir + '/IBC'
atlas_dir = base_dir + '/Atlases'


def show_ibc_group(ses_id='ses-hcp1', type='CondHalf', atlas='MNISymC3',
                   cond=0, info_column='names', savefig=False):
    if (atlas == 'MNISymC3'):
        mask = atlas_dir + '/tpl-MNI152NLIn2000cSymC/tpl-MNISymC_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum', mask_img=mask)
    
    ibc_dataset = DataSetIBC(data_dir)
    C = nb.load(ibc_dataset.data_dir.split('/{0}')[0] +
                f'/group/group_{ses_id}_space-{atlas}_{type}.dscalar.nii')
    D = pd.read_csv(ibc_dataset.data_dir.split('/{0}')[0] +
                    f'/group/group_{ses_id}_info-{type}.tsv', sep='\t')
    X = C.get_fdata()

    if cond == 'all':
        conditions = D[info_column]
        # -- each in seperate figures --
        dest_dir = ibc_dataset.data_dir.split('/{0}')[0] + f'/group/figures/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        for i, c in enumerate(conditions):
            Nifti = suit_atlas.data_to_nifti(X[i, :])
            surf_data = suit.flatmap.vol_to_surf(Nifti, atlas[:-1])
            fig = suit.flatmap.plot(
                surf_data, render='matplotlib', new_figure=True)
            fig.set_title(c)
            # save figure
            if savefig:
                plt.savefig(dest_dir + f'group_{ses_id}_{c}.png')
            plt.clf()
            pass

    else:
        Nifti = suit_atlas.data_to_nifti(X[cond, :])
        surf_data = suit.flatmap.vol_to_surf(Nifti)
        fig = suit.flatmap.plot(surf_data, render='plotly')
        fig.show()
        print(f'Showing {D.cond_name[cond]}')
        pass

def extract_all(atlas='MNISym3'):
    ibc_dataset = DataSetIBC(data_dir)
    info = ibc_dataset.get_participants()
    for ses in ibc_dataset.sessions[9:]:
        print(f'extracting {ses}')
        if atlas == 'fs32k':
            ibc_dataset.extract_all_fs32k(ses,type='CondHalf')
        else:
            ibc_dataset.extract_all_suit(ses,type='CondHalf',atlas=atlas)

def smooth_ibc_fs32k(type='CondHalf', smooth=1):
    myatlas, _ = am.get_atlas('fs32k', atlas_dir)
    ds = DataSetIBC(data_dir)
    T = ds.get_participants()

    # get the surfaces for smoothing
    surf_L = ds.atlas_dir + f'/tpl-fs32k/fs_LR.32k.L.midthickness.surf.gii'
    surf_R = ds.atlas_dir + f'/tpl-fs32k/fs_LR.32k.R.midthickness.surf.gii'

    for ses_id in ds.sessions:
        for s in T.participant_id:
            print(f'- Smoothing data for {s} fs32k {ses_id} in {smooth} mm ...')
            # Load the unsmoothed data and fill nan with zeros
            C = nb.load(ds.data_dir.format(s)
                        + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii')
            mask = np.isnan(C.get_fdata())
            C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
            nb.save(C, 'tmp.dscalar.nii')

            dest_dir = ds.data_dir.format(s)
            cifti_out = dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_' \
                                   f'desc-sm{smooth}.dscalar.nii'

            # Write in smoothed surface data (filled with 0)
            smooth_cmd = f"wb_command -cifti-smoothing tmp.dscalar.nii " \
                         f"{smooth} {smooth} COLUMN {cifti_out} " \
                         f"-left-surface {surf_L} -right-surface {surf_R} " \
                         f"-fix-zeros-surface"
            subprocess.run(smooth_cmd, shell=True)
            os.remove("tmp.dscalar.nii")

            # Replace 0s back to NaN (we don't want the 0s impact model learning)
            C = nb.load(cifti_out)
            data = C.get_fdata()
            data[mask] = np.nan
            C = nb.Cifti2Image(dataobj=data, header=C.header)
            nb.save(C, cifti_out)

def copy_currentAsOld():
    """This function copies the regressor info file as "_old.tsv"
    and creat a new one with correct format in the current name.

    Returns:
        None
    Notes:
        I'm not sure what "ses-self" does, so I just correct all
        separate sessions reginfo file. Ana, if you want to correct
        the file in ses-self, you can easily modify this function
        to do so.           --dzhi
    """
    ibc_dataset = DataSetIBC(base_dir + '/IBC')
    sess = ibc_dataset.sessions
    T = ibc_dataset.get_participants()
    for sub in T.participant_id:
        for s in sess:
            print(f'Correcting {sub} {s} ...')
            dirw = ibc_dataset.estimates_dir.format(sub) + f'/{s}'
            df = pd.read_csv(dirw + f'/{sub}_{s}_reginfo.tsv', delimiter='\t')

            # 1. Save an old copy
            df.to_csv(dirw + f'/{sub}_{s}_reginfo_old.tsv' ,sep='\t', index=False)

            # 2. convert all values to int (As there are float
            # values in old file, we don't want that)
            df = df.astype({'reg_id': 'int', 'reg_num': 'int',
                            'half': 'int', 'cond_num_uni': 'int'})

            # 4. Switch the column name `reg_id` and `reg_num`
            df.rename(columns={'reg_id': 'reg_num', 'reg_num': 'reg_id'}, inplace=True)

            # Save the new file
            df.to_csv(dirw + f'/{sub}_{s}_reginfo.tsv' ,sep='\t', index=False)
            del df

    pass

def correct_condHalf():
    """This function corrects the Cond 'half' to be either
    1 or 2, instead of multiple numbers in IBC

    Returns:
        None
    Notes:
        The reason of not integrating this function in the
        copy_currentAsOld() is that we may want to have another
        ways of dealing with this multiple `half`. Right now,
        it just assign all odd number as 1, evens as 2 which not
        sure is the ideal solution.         --dzhi
    """
    ibc_dataset = DataSetIBC(base_dir + '/IBC')
    sess = ibc_dataset.sessions
    T = ibc_dataset.get_participants()
    for sub in T.participant_id:
        for s in sess:
            print(f'Correcting {sub} {s} ...')
            dirw = ibc_dataset.estimates_dir.format(sub) + f'/{s}'
            df = pd.read_csv(dirw + f'/{sub}_{s}_reginfo.tsv', delimiter='\t')

            # change all odd value in `half` to be 1, all even to be 2
            df.loc[df['half'] % 2 == 1, 'half'] = 1
            df.loc[df['half'] % 2 == 0, 'half'] = 2

            # Save the new file
            df.to_csv(dirw + f'/{sub}_{s}_reginfo.tsv', sep='\t', index=False)
            del df

    pass

if __name__ == "__main__":
    # copy_currentAsOld()
    # correct_condHalf()
    # extract_all('SUIT3')
    # extract_all('MNISymC3')
    # extract_all('MNISymC2')
    # extract_all('MNISymC2')

    ################# Smooth IBC fs32k data #################
    for s in [1,2,3,4,5,6,7]:
        smooth_ibc_fs32k(type='CondHalf', smooth=s)

    # dataset = DataSetIBC(data_dir)
    # for session in dataset.sessions:
    #     dataset.group_average_data(atlas='MNISymC2', ses_id=session)
    #
    # dataset.plot_cerebellum(savefig=True, atlas='MNISymC2', colorbar=True)

 
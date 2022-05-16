# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import shutil
from pathlib import Path
import mat73
import numpy as np
import atlas_map as am
from dataset import DataSetMDTB
import nibabel as nb
import SUITPy as suit


base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'

data_dir = base_dir + '/MDTB'
atlas_dir = base_dir + '/Atlases'

def get_mdtb_suit(ses_id='ses-s1',type='CondSes'):
    # Make the atlas object
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    # initialize the data set object
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    for s in T.participant_id:
        print(f'Atlasmap {s}')
        deform = mdtb_dataset.suit_dir.format(s) + f'/{s}_space-SUIT_xfm.nii'
        mask = mdtb_dataset.suit_dir.format(s) + f'/{s}_desc-cereb_mask.nii'
        atlas_map = am.AtlasMapDeform(mdtb_dataset, suit_atlas, s,deform, mask)
        atlas_map.build(smooth=2.0)
        print(f'Extract {s}')
        data,info,names = mdtb_dataset.get_data(s,[atlas_map],
                                                ses_id=ses_id,
                                                type=type)
        C=am.data_to_cifti(data,[atlas_map],names)
        dest_dir = mdtb_dataset.data_dir.format(s)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        nb.save(C, dest_dir + f'/{s}_space-SUIT3_{ses_id}_{type}.dscalar.nii')
        info.to_csv(dest_dir + f'/{s}_{ses_id}_info-{type}.tsv',sep='\t')

def show_mdtb_suit(subj,sess,cond): 
    mask = atlas_dir + '/tpl-SUIT/tpl-SUIT_res-3_gmcmask.nii'
    suit_atlas = am.AtlasVolumetric('cerebellum',mask_img=mask)
    mdtb_dataset = DataSetMDTB(data_dir)
    T = mdtb_dataset.get_participants()
    s = T.participant_id[subj]
    ses = f'ses-s{sess}'
    C = nb.load(mdtb_dataset.data_dir.format(s) + f'/{s}_space-SUIT3_{ses}_CondSes.dscalar.nii')
    D = pd.read_csv(mdtb_dataset.data_dir.format(s) + f'/{s}_{ses}_info-CondSes.tsv',sep='\t')
    X = C.get_fdata()
    Nifti = suit_atlas.data_to_nifti(X)
    surf_data = suit.flatmap.vol_to_surf(Nifti)
    fig = suit.flatmap.plot(surf_data[:,cond],render='plotly')
    fig.show()
    print(f'Showing {D.cond_name[cond]}')
    pass 

def get_mdtb_fs32k():
    # Make the atlas object
    atlas =[]
    bm_name = ['cortex_left','cortex_right']
    for i,hem in enumerate(['L','R']):
        mask = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-{hem}_mask.label.gii'
        atlas.append(am.AtlasSurface(bm_name[i],mask_gii=mask))
    # initialize the data set object
    mdtb_dataset = DataSetMDTB(data_dir)

    # create and calculate the atlas map for each participant
    T = mdtb_dataset.get_participants()
    for s in T.participant_id:
        atlas_maps = []
        data = []
        for i,hem in enumerate(['L','R']):
            adir = mdtb_dataset.anatomical_dir.format(s)
            pial = adir + f'/{s}_space-32k_hemi-{hem}_pial.surf.gii'
            white = adir + f'/{s}_space-32k_hemi-{hem}_white.surf.gii'
            mask = adir + f'/{s}_desc-brain_mask.nii'
            atlas_maps.append(am.AtlasMapSurf(mdtb_dataset, atlas[i],
                            s,white,pial, mask))
            atlas_maps[i].build()
            # data = mdtb_dataset.get_data(s,[A])
            # data_files=mdtb_dataset.get_data_fnames(s,'ses-s1')
            data.append(np.random.normal(0,1,(100,atlas_maps[i].P))) # am.get_data(data_files,atlas_maps)
        im = am.data_to_cifti(data,atlas_maps)
        nb.save(im,atlas_dir + '/tpl-fs32k/tpl_gs32k_func.dscalar.nii')
        pass


if __name__ == "__main__":
    get_mdtb_suit(ses_id='ses-s1',type='CondRun')
    pass


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
# Script for importing the Pontine data set to general format.
import pandas as pd
from pathlib import Path
import numpy as np
import sys, os, time
import Functional_Fusion.atlas_map as am
import Functional_Fusion.dataset as ds
import Functional_Fusion.util as fut
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import ants
import nitools as nt

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'


def ants_transform_to_deformation_field(xfm,src_img,ref_img):
    """Converts the ANTs transform to a deformation field"""
    ref_ants = ants.image_read(ref_img)
    src_img = nb.load(src_img)

    # Make a image in src space with the voxel coordinates 
    Xv,Yv,Zv = np.meshgrid(range(src_img.shape[0]),
                           range(src_img.shape[1]),
                           range(src_img.shape[2]),indexing='ij')
    X,Y,Z = nt.affine_transform(Xv,Yv,Zv,src_img.affine)
    data = np.stack([X,Y,Z],axis=3)
    coord_img  = nb.Nifti1Image(data,src_img.affine)
    coord_ants = ants.from_nibabel(coord_img)

    # Apply the transform to this coordinate image
    y_ants = ants.apply_transforms(ref_ants, coord_ants, xfm, 
                          interpolator='linear', 
                          imagetype=3, 
                          defaultvalue=np.nan)
    # Insert a dummy 4th dimension
    y_img = y_ants.to_nibabel()
    data = y_img.get_fdata()
    data = np.expand_dims(data,axis=3)
    y_img  = nb.Nifti1Image(data,y_img.affine)
    return y_img

def normalize_bold_all(dataset,template,mask=None,tot='SyNCC',kwargs={},subj=None): 
    if subj is None:
        subj = dataset.get_participants().participant_id
    
    # Temporary directories for cehcking (delete later)
    def_dir = f'{dataset.base_dir}/deformed/'
    xfm_dir = f'{dataset.base_dir}/transforms/'
    if not os.path.isdir(def_dir):
        os.mkdir(def_dir)
    if not os.path.isdir(xfm_dir):
        os.mkdir(xfm_dir)

    trg = f'{base_dir}/Atlases/{template}'
    trg_img = ants.image_read(trg)
    for src_sub in subj[1:]:
        print(f'Processing {src_sub}')
        src = f'{dataset.anatomical_dir.format(src_sub)}/{src_sub}_meanbold.nii'
        prefix = f'{xfm_dir}/xfm_{src_sub}_'
        src_img = ants.image_read(src)
        mytx = ants.registration(fixed=trg_img , 
            moving=src_img,
            type_of_transform=tot,
            outprefix=prefix,
            write_composite_transform=True,
            **kwargs)
        fname_w = f'{def_dir}/w{src_sub}.nii'
        wsrc_img = mytx['warpedmovout']
        ants.image_write(wsrc_img,fname_w)
        
        # Translater into deformation field
        xfm = f'{xfm_dir}/xfm_{src_sub}_composite.h5'
        def_field = ants_transform_to_deformation_field([xfm],src,trg)
        fname=f'{dataset.anatomical_dir.format(src_sub)}/{src_sub}_space-MNI2009cSym_xfm.nii'
        nb.save(def_field,fname)
        pass


if __name__ == "__main__":
    mdtb = ds.get_dataset_class(base_dir,'MDTB')
    normalize_bold_all(mdtb,'tpl-MNI152NLin2009cSym/tpl-MNI152NLin2009cSym_bold.nii')


    pass

# Script for distorting a new SUIT surface for the MNISymC template
import numpy as np 
import nibabel as nb
import SUITPy as suit
import nitools as nt

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
atlas_dir = base_dir + '/Atlases'
suit_dir = '/Users/jdiedrichsen/Python/SUITPy/SUITPy/surfaces'

def deform_surface(surface,xfm_file):
    gifti=nb.load(surface)
    xfm_img = nb.load(xfm_file)
    coords = gifti.agg_data('NIFTI_INTENT_POINTSET')
    new_coords = nt.sample_image(xfm_img,coords[:,0],coords[:,1],coords[:,2],1)
    new_coords=new_coords.squeeze()
    new_coords[np.isnan(new_coords)]=0
    gifti.darrays[0].data=new_coords.astype(np.float32)
    return gifti

def deform_suit_surfaces():
    adir = atlas_dir +'/tpl-MNI152NLin2000cSymC'
    def_name = adir + '/tpl-MNI152NLin2009cSymC_space-SUIT_xfm.nii'
    in_files = ['PIAL_SUIT.surf.gii','WHITE_SUIT.surf.gii']
    out_files = ['PIAL_MNISymC.surf.gii','WHITE_MNISymC.surf.gii']
    for i,fn in enumerate(in_files):
        gii=deform_surface(suit_dir + '/' + fn,def_name)
        nb.save(gii,suit_dir + '/' + out_files[i])
    pass 


if __name__ == "__main__":
    # reslice_SUIT()
    deform_suit_surfaces()

    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass
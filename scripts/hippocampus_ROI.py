# Builds a single hippocampus functional ROI dseg image (ROIs 1-5, left and right
# combined) from the functional 5-cluster from the PLatchi atlas , using the
# MNIAsymHippocampus_L / _R atlases.

import numpy as np
import nibabel as nb
from scipy.spatial import cKDTree
import Functional_Fusion.atlas_map as am

data_dir = "/Users/jdiedrichsen/Diedrichsenlab Dropbox/Joern Diedrichsen/projects/Hippocampus_ROI"
n_clusters = 5

def hem_label_vector(atlas, hem):
    label_vector = np.zeros((atlas.P,), dtype=np.int32)
    for i in range(1, n_clusters + 1):
        fname = f"{data_dir}/{hem}h_func_{n_clusters}solution_Cluster_F{i}.nii.gz"
        data = atlas.read_data(fname, interpolation=0)
        label_vector[data > 0] = i

    # Assign unlabeled atlas voxels the ROI of the spatially closest labeled voxel
    assigned = label_vector > 0
    coords = atlas.world.T
    tree = cKDTree(coords[assigned])
    _, nearest = tree.query(coords[~assigned])
    label_vector[~assigned] = label_vector[assigned][nearest]
    return label_vector

def make_dseg():
    atlas_l, _ = am.get_atlas("MNIAsymHippocampus_L")
    atlas_r, _ = am.get_atlas("MNIAsymHippocampus_R")
    label_l = hem_label_vector(atlas_l, "L")
    label_r = hem_label_vector(atlas_r, "R")

    img_l = atlas_l.data_to_nifti(label_l)
    img_r = atlas_r.data_to_nifti(label_r)
    X = img_l.get_fdata().astype(np.int32) + img_r.get_fdata().astype(np.int32)
    dseg_img = nb.Nifti1Image(X, img_l.affine)

    out_name = f"{data_dir}/atl-Platchi{n_clusters}_space-MNI152NLin6Asym_dseg.nii"
    nb.save(dseg_img, out_name)
    print(f"Saved {out_name}")

if __name__ == "__main__":
    make_dseg()

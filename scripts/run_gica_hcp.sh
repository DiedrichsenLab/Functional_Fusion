#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  run_gica_hcp.sh
#
# Description:  Script to run group ICA on HCP resting-state data
#               
#
# Author:       Caroline Nettekoven, 2022
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# Dependencies
if [ -d /Volumes/diedrichsen_data$/data ]; then
workdir=/Volumes/diedrichsen_data$/data
elif [ -d /srv/diedrichsen/data ]; then 
workdir=/srv/diedrichsen/data
else
echo "Workdir not found. Mount or connect to server and try again."
fi
# --------------------------------------------------------------------------------------------------------
output_dir=${workdir}/FunctionalFusion/HCP/group_ica
script_dir=~/code/functional_fusion
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --- Make input list for niftis --- 
ls ${workdir}/HCP_UR100_rfMRI/*/MNINonLinear/Results/rfMRI_REST?_??/rfMRI_REST?_??_hp2000_clean.nii.gz > ${script_dir}/input_list_hcp_rsfmri.txt


# --- Run group ICA --- 
tr=0.72

dim=50

fsl_sub \
melodic \
-i ${script_dir}/input_list_hcp_rsfmri.txt \
-o ${output_dir}/dim_${dim} \
--tr=${tr} \
--nobet \
-a concat \
-m ${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz \
--report \
--Oall \
-d ${dim}


fsl_sub \
melodic \
-i ${script_dir}/input_list_hcp_rsfmri.txt \
-o ${output_dir}/dim_auto \
--tr=${tr} \
--nobet \
-a concat \
-m ${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask.nii.gz \
--report \
--Oall
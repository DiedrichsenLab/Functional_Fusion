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
tr=0.72
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# --- Make input list for niftis ---
ls ${workdir}/HCP_UR100_rfMRI/*/MNINonLinear/Results/rfMRI_REST?_??/rfMRI_REST?_??_hp2000_clean.nii.gz > ${script_dir}/input_list_hcp_rsfmri.txt


# --- Run group ICA ---
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

# --- Extract all components from melodic_IC ---
fslsplit \
${output_dir}/dim_auto/melodic_IC.nii.gz ${output_dir}/dim_auto/component_ -t

# --- Merge signal components from melodic_IC ---
cat ${output_dir}/dim_auto/classified_components.txt | grep "Signal" | awk '{print $1}' > ${output_dir}/dim_auto/signal_components.txt

for i in `cat ${output_dir}/dim_auto/signal_components.txt`; do
component=${i%%,}
component=$(( component - 1 ))
echo  ${output_dir}/dim_auto/components/component_`zeropad ${component} 4`.nii.gz >> ${output_dir}/dim_auto/signal_components_imgs.txt
done

fslmerge -t ${output_dir}/dim_auto/signal_components.nii.gz `cat ${output_dir}/dim_auto/signal_components_imgs.txt`

# --- Resample into 1 mm volume space ---
# flirt \
# -in ${output_dir}/dim_auto/signal_components.nii.gz \
# -ref ${workdir}/FunctionalFusion/Atlases/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_T1w.nii \
# -out ${output_dir}/dim_auto/signal_components_1mm.nii.gz \
# -nosearch \
# -applyisoxfm 1 \
# -usesqform


# --- Resample all signal componnents into surface space ---

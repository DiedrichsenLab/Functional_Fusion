#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  rest_mdtb.sh
#
# Description:
#               Script to preprocess MDTB Resting state data
#
# Author:       Caroline Nettekoven, 2022
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -d /srv/diedrichsen/data ]; then
base_dir=/srv/diedrichsen/data
elif [ -d /Volumes/diedrichsen_data$/data ]; then
base_dir=/Volumes/diedrichsen_data$/data 
else
echo "Server directory not found."
fi

data_dir=${base_dir}/Cerebellum/super_cerebellum/resting_state/imaging_data/
design_dir=~/code/Python/Functional_Fusion/preprocessing/design_files/
runs=("01" "02")
# ----
# ------------------------------------------------------------------------------
# Step 1: Add TR to header of nifti files (needed for single-subject ICA with fsl)
for i in ${data_dir}/s*; do
    for run in $runs; do
        subject=s${i#*imaging_data//s}
        if [ ! -f ${data_dir}/${subject}/rrun_${run}.nii.gz ] && [ -f ${data_dir}/${subject}/rrun_${run}.nii ]; then
            echo "Adding TR to header of ${subject} ${run}"
            fslhd -x ${data_dir}/${subject}/rrun_${run}.nii > myhdr.txt
            sed "s/dt =.*/dt = \'1\'/" myhdr.txt > myhdr2.txt
            fslcreatehd myhdr2.txt ${data_dir}/${subject}/rrun_${run}
            rm myhdr.txt myhdr2.txt
            else
            echo "${subject} ${run} incomplete"
            echo " "
        fi
    done
done    
# ------------------------------------------------------------------------------

for i in ${data_dir}/s*; do
    for run in $runs; do
        subject=${i#*imaging_data//s}
        if [ -f ${data_dir}/${subject}/rrun_${run}.nii.gz ] && [ ! -d ${data_dir}/${subject}/run${run}.ica ]; then
            
            cp ${design_dir}/ssica_template.fsf ${design_dir}/rest_${subject}_run-${run}.fsf
            # Fill in subject-and-run-specific parameters
            sed -i  "s/XX/${subject}/g" "${design_dir}/rest_${subject}_run-${run}.fsf"
            sed -i  "s/YY/${run}/g" "${design_dir}/rest_${subject}_run-${run}.fsf"

            if [ ! -d ${data_dir}/${subject}/run${run}.ica/ ] ; then
                echo "Running SS melodic for subjec ${subject} ${run}"
                feat "${design_dir}/rest_${subject}_run-${run}.fsf"   
            else
                echo "SS melodic already run for ${subject} ${run}"
            fi

            if [ -d ${data_dir}/${subject}/run${run}.ica/ ]; then 
            firefox ${data_dir}/${subject}/run${run}.ica/report.html &
            fi

            else
            echo "${subject} ${run} incomplete"
            echo " "
        fi
    done
done
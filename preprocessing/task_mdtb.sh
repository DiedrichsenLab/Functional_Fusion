#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  task_mdtb.sh
#
# Description:
#               Script to preprocess MDTB Task like rest data
#
# Author:       Caroline Nettekoven, 2024
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

data_dir=${base_dir}/Cerebellum/super_cerebellum/sc1/imaging_data/
design_dir=~/code/Python/Functional_Fusion/preprocessing/design_files/
runs=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32")
# ----
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

for i in ${data_dir}/s*; do
    for run in $runs; do
        subject=${i#*imaging_data//s}
        if [ -f ${data_dir}/${subject}/rrun_${run}.nii.gz ] && [ ! -d ${data_dir}/${subject}/run${run}.ica ]; then
            
            cp ${design_dir}/ssica_template_task.fsf ${design_dir}/task_${subject}_run-${run}.fsf
            # Fill in subject-and-run-specific parameters
            sed -i  "s/XX/${subject}/g" "${design_dir}/task_${subject}_run-${run}.fsf"
            sed -i  "s/YY/${run}/g" "${design_dir}/task_${subject}_run-${run}.fsf"

            if [ ! -d ${data_dir}/${subject}/run${run}.ica/ ] ; then
                echo "Running SS melodic for subjec ${subject} ${run}"
                feat "${design_dir}/task_${subject}_run-${run}.fsf"   
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

#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  bids_somatotopic.sh
#
# Description:
#               Script to sort the somatotopic data into BIDS-compatible format
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
# ----
cd ${SCRIPTDIR}

src_data_dir=${base_dir}/Cerebellum/Somatotopic/raw/originals/MNI_defaced_8subjects
dest_data_dir=${base_dir}/Cerebellum/Somatotopic/raw/


# ------ Register MNI space ROIs (1mm resolution MNI) to Degeneration_40 template ------
for image in ${src_data_dir}/S*.nii.gz; do
    sub=${image##*/S}
    sub=${sub%%_*}
    bids_name=`zeropad ${sub} 2`
    echo ${image} ${bids_name}


    if [ ! -d ${dest_data_dir}/sub-${bids_name}/anat ]; then
    mkdir -p ${dest_data_dir}/sub-${bids_name}/anat
    fi

    # mv ${image} ${dest_data_dir}/sub-${sub}/anat/.
    mv ${image} ${dest_data_dir}/sub-${bids_name}/anat/sub-${bids_name}_T1w.nii.gz

done

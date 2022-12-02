#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  find_ac_coordinates.sh
#
# Description:
#               Script to find AC coordinates for anatomical data of a whole dataset
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

data_dir=${base_dir}/Cerebellum/Somatotopic/raw/


# ------ Register MNI space ROIs (1mm resolution MNI) to Degeneration_40 template ------
for subject in ${data_dir}/sub-*; do
sub=${subject##*sub-}
# strip leading zeros (reverse zero-padding)
sub=$((10#${sub}))
anat_image=S${sub}_anat_mni_underlay_defaced.nii.gz

fsleyes ${subject}/anat/${anat_image}

done

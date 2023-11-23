#!/bin/bash
# ------------------------------------------------------------------------------
# Script name:  classify_components.sh
#
# Description:
#               Script to loop through MDTB single-subject ICA components and classify them
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
subjects=(`ls -d ${data_dir}/s??`)
runs=("01" "02")
# ----

# ------ Classify single-subject ICA components ------ 
for i in $subjects; do
    subject=${i##*s}
    # test if subject is in ${design_dir}/classified_subjects.csv
    if grep -q s${subject} ${design_dir}/classified_subjects.csv; then
        run=`grep s${subject} ${design_dir}/classified_subjects.csv | awk '{print $2}'`
        if [ -d $i/run${run}.ica/filtered_func_data.ica ] && [ ! -f $i/run${run}.ica/filtered_func_data.ica/classifications.txt ]; then
                echo CLassifying $i ${run}
                
                # Open motion parameter plots either from report_prestats.html or generate them from rp_run_${run}.txt
                if [ -f  $i/rp_run_${run}.txt ] && [ ! -f $i/rot_${run}.png ]; then
                    fsl_tsplot -i $i/rp_run_${run}.txt -t 'SPM estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 800 -h 300 -o $i/rot_${run}.png
                    fsl_tsplot -i $i/rp_run_${run}.txt -t 'SPM estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 800 -h 300 -o $i/trans_${run}.png
                    fsl_tsplot -i $i/rp_run_${run}.txt -t 'SPM estimated rotations and translations (mm)' -u 1 --start=1 --finish=6 -a "x(rot),y(rot),z(rot),x(trans),y(trans),z(trans)" -w 800 -h 300 -o $i/rot_trans_${run}.png
                fi
                if [ -f $i/rot_trans_${run}.png ]; then
                    open $i/*_${run}.png
                else 
                    open $i/run${run}.ica/report_prestats.html
                fi
                
                fsleyes --scene melodic -ad $i/run${run}.ica/filtered_func_data.ica

            else
                echo "Already classified $i ${run}"
        fi
    else
        echo "Skipping $i"
    fi
done


# Ploting the movement parameters for scan processed with SPM:
# The SPM-output text file starting with rp_ and containing the run number (01 or 02) contains the movement parameters of the scan (previously estiated during SPM motion correction).
# The file contains six columns. These represent position (not movement) of the subject in terms of X, Y, and Z displacement and X, Y, and Z axis rotation with respect to the first scan entered.
# The first three columns are in millimeters and the last three are in degrees.

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

# # Rename all .nii.gz files to end in _hdr.nii.gz
# for i in /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/imaging_data/s*/*.nii.gz; do
#     mv "$i" "${i%.nii.gz}_hdr.nii.gz"
# done




 for i in /Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s*; do
    echo $i;
    fsleyes --voxelLoc 91 73 147 $i/run01.feat/reg/highres $i/run01.feat/reg/example_func ;
done





i=/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s10
echo $i;
fsleyes --voxelLoc 91 73 147 $i/run01.feat/reg/highres $i/run01.feat/reg/example_func 



# ------ CEREBELLUM ----------------


fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s06/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s06/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s06-run01_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s06-run01_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &


fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s07/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s07/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s07-run01_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s07-run01_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &


fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s08/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s08/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s08-run02_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s08-run02_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &


fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s10/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s10/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s10-run02_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/c_s10-run02_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &





# ------ OCCIPITAL ----------------

fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s06/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s06/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s06-run01_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s06-run01_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &

fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s07/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s07/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s07-run01_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s07-run01_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &

fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s08/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s08/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s08-run02_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s08-run02_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &

fsleyes \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s10/run01.feat/reg/example_func.nii.gz \
/Volumes/diedrichsen_data$/data/FunctionalFusion/../Cerebellum/super_cerebellum/resting_state/imaging_data/s10/run01.feat/reg/highres.nii.gz \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s10-run02_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
/Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_s10-run02_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &


# Pull the images up automatically
for i in /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/imaging_data/s1*; do
    echo $i;
    subject=${i##*imaging_data/}
    subject=${subject%/}
    fsleyes --voxelLoc 91 73 147 $i/run01.feat/reg/highres $i/run01.feat/reg/example_func  \
    /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_${subject}-run01_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
    /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_${subject}-run01_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &
done


# Take screenshots automatically
for i in /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/imaging_data/s*; do
    echo $i;
    subject=${i##*imaging_data/}
    subject=${subject%/}
    fsleyes render -of ~/Documents/Projects/LocalizationCerebellum/figures/corrmaps/${subject}_prefix.png --voxelLoc 91 73 147 $i/run01.feat/reg/highres $i/run01.feat/reg/example_func  \
    /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_${subject}-run01_prefix.nii.gz -dr 0.1 0.3 -cm red-yellow \
    
    fsleyes render -of ~/Documents/Projects/LocalizationCerebellum/figures/corrmaps/${subject}_postfix.png --voxelLoc 91 73 147 $i/run01.feat/reg/highres $i/run01.feat/reg/example_func  \
    /Volumes/diedrichsen_data$/data/Cerebellum/super_cerebellum/resting_state/fix_ica/corrmaps/o_${subject}-run01_postfix.nii.gz -dr 0.1 0.3 -cm red-yellow &
done
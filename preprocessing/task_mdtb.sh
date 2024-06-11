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






# ------- Check if melodic_IC is empty -------
for i in `ls -d ${data_dir}/s*/*.feat/`; do
         # File to check
        file=${i}/filtered_func_data.ica/melodic_IC.nii.gz

        # Check if the file is 0 bytes
        if [ ! -s "$file" ]; then
            echo "\n\n============= $file is empty. rerunning melodic...============="
            
            # Change directory to .feat folder
            cd ${i}

            # Run melodic
            melodic \
            -i filtered_func_data \
            -o filtered_func_data+.ica \
            -v --nobet --bgthreshold=1 --tr=1.0 -d 0 --mmthresh=0.5 
        else 
            echo "$file SUCCESS"
        fi

    done
done



# ------- Run Feat if no feat directory exists -------
cd ${design_dir}
# Only loop through subjects 02-05

# for n in `seq 2 5`; do --> RUNNING last on Caro GPU
# for n in `seq 6 10`; do --> RUNNING on Caro GPU
# for n in `seq 11 15`; do  --> RUNNING on Marco Basic
# for n in `seq 16 20`; do --> RUNNING on Marco Heavy
# for n in `seq 21 25`; do --> RUNNING on Marco GPU
# for n in `seq 26 30`; do --> RUNNING on Caro Heavy
# for n in `seq 31 33`; do --> RUNNING on Caro Basic
    subject=`zeropad $n 2`
    for run in `seq 1 32`; do
        run=`zeropad $run 2`
        if [ -f ${data_dir}/s${subject}/rrun_${run}_hdr.nii.gz ] && [ ! -d ${data_dir}/s${subject}/run${run}.feat ] && [ ! -d ${data_dir}/s${subject}/run${run}.correct ]; then
            echo "Running Feat for subject ${subject} run ${run}"
            feat ssica_task_${subject}_run-${run}.fsf
            # Remove superfluous files after feat
            rm ${data_dir}/s${subject}/run${run}.feat/prefiltered_func_*
        elif [ -f ${data_dir}/s${subject}/rrun_${run}_hdr.nii.gz ] ; then
            echo "Feat already run for subject ${subject} run ${run}"
        
        fi
    done
done


# for i in `ls ${data_dir}/s*/*hdr.nii.gz`; do
#          # Feat folder to check
#         run=${i#*rrun_}
#         run=${run%_hdr.nii.gz}
#         subject=${i#*imaging_data//s}
#         subject=${subject%/*}

#         # Check if the file is 0 bytes
#         if [ ! -s "$file" ]; then
#             echo "\n\n============= $file is empty. rerunning melodic...============="

#             # Run feat
#             # feat \
#             ls ssica_task_${subject}_run-${run}.fsf
#         else 
#             echo "$file SUCCESS"
#         fi

# done


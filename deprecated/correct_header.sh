
# Step 1: Add TR to header of nifti files (needed for single-subject ICA with fsl)
# (Deprecated since this is done by the Python function now)
for i in ${data_dir}/s*; do
    for run in $runs; do
        subject=s${i#*imaging_data//s}
        if [ ! -f ${data_dir}/${subject}/rrun_${run}.nii.gz ] && [ -f ${data_dir}/${subject}/rrun_${run}.nii ]; then
            echo "Adding TR to header of ${subject} ${run}"
            fslhd -x ${data_dir}/${subject}/rrun_${run}.nii > myhdr.txt
            sed "s/dt =.*/dt = \'1\'/" myhdr.txt > myhdr2.txt
            fslcreatehd myhdr2.txt ${data_dir}/${subject}/rrun_${run}_hdr
            rm myhdr.txt myhdr2.txt
            else
            echo "${subject} ${run} incomplete"
            echo " "
        fi
    done
done    
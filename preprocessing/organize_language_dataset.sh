raw_dir=/Volumes/diedrichsen_data$/data/FunctionalFusion/Language/raw

for i in ${raw_dir}/*.nii; do
    sub=${i##*sub-} 
    sub=${sub%%_*}
    echo ${sub}
    if [ ! -d ${raw_dir}/sub-${sub} ]; then
    mkdir ${raw_dir}/sub-${sub} 
    fi
    mv ${raw_dir}/sub-${sub}*.nii ${raw_dir}/sub-${sub}/.
done

cd ${raw_dir}
echo "participant_id" > ${raw_dir}/../participants.tsv
ls -1 -d sub-* >> ${raw_dir}/../participants.tsv

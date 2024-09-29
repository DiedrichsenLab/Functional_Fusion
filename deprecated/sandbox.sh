# base_dir="/Volumes/diedrichsen_data$/data/"
base_dir="/cifs/diedrichsen/data/"
for i in ${base_dir}/FunctionalFusion/MDTB/derivatives/sub*/data/*FixBeta*; do 
    echo $i
    # Rename the "FixBeta" part of the file to "FixNewCondRun"
    newname=$(echo $i | sed 's/FixBeta/FixNewCondRun/')
    echo $newname
    mv $i $newname

done

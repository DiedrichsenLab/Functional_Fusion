# base_dir="/Volumes/diedrichsen_data$/data/"
base_dir="/cifs/diedrichsen/data/"
for i in ${base_dir}/FunctionalFusion/MDTB/derivatives/sub*/data/*FixBeta*; do 
    echo $i
    # Rename the "FixBeta" part of the file to "FixNewCondRun"
    newname=$(echo $i | sed 's/FixBeta/FixNewCondRun/')
    echo $newname
    mv $i $newname

done


base_dir="/cifs/diedrichsen/data/"
for i in ${base_dir}/FunctionalFusion/MDTB/derivatives/sub*/data/*ses-rest_Tseries*; do 
    echo $i
    # Rename the "FixBeta" part of the file to "FixNewCondRun"
    newname=$(echo $i | sed 's/Tseries/FixTseries/')
    echo $newname
    mv $i $newname

done


base_dir="/cifs/diedrichsen/data/"
for i in ${base_dir}/FunctionalFusion/MDTB/derivatives/sub*/data/*ses-rest_*Net69Run*; do 
    echo $i
    # Rename the "FixBeta" part of the file to "FixNewCondRun"
    newname=$(echo $i | sed 's/Net69Run/Net69FixRun/')
    echo $newname
    mv $i $newname

done


# for all files named '/Volumes/diedrichsen_data$/data//FunctionalFusion/MDTB/derivatives/sub-*/estimates/ses-rest/sub-*_ses-restrrun_*.nii', rename as /Volumes/diedrichsen_data$/data//FunctionalFusion/MDTB/derivatives/sub-*/estimates/ses-rest/sub-*_ses-rest_run-*.nii
base_dir="/cifs/diedrichsen/data/"
for i in ${base_dir}/FunctionalFusion/MDTB/derivatives/sub*/estimates/ses-rest/sub-*_ses-restrrun_*.nii; do 
    echo $i
    # Rename the "FixBeta" part of the file to "FixNewCondRun"
    newname=$(echo $i | sed 's/restrrun_/rest_run-/')
    echo $newname
    mv $i $newname

done
for i in ../imaging_data/s*/*/hand_labels_noise.txt; do
    echo $i; n=`cat $i | wc | awk '{print $1}'`;
    echo $((n-2));
    fslinfo ${i%%hand*}/filtered_func_data.ica/melodic_IC.nii.gz | grep "dim4" | head -1 | awk '{print $2}' ;
done

for i in *.csv; do 
cat $i | wc | awk '{print $1}';
n=`cat ${i%%.csv}.txt | wc | awk '{print $1}'` ; echo $((n-2)); 
done


for i in *.csv; do 
cat $i | wc | awk '{print $1}';
tail -2 ${i%%.csv}.txt | head -1 | awk '{print $1}';
done
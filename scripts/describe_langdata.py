import pandas as pd

# run in bash: ls /Volumes/diedrichsen_data$/data/FunctionalFusion/Language/raw/*.nii > filenames.txt
# add header to filenames.txt
lang_info = '/Volumes/diedrichsen_data$/data/FunctionalFusion/Language/raw/filenames.txt'
info = pd.read_csv(lang_info,sep='_', header=0, names=['subject', 'task', 'contrast', 'imagetype'])

info.subject = info.subject.astype('category')
info.task = info.task.astype('category')
info.contrast = info.contrast.astype('category')
info.imagetype = info.imagetype.astype('category')

info.subject.unique() # n=41
info.task.unique() # 5 tasks : ['MDloc', 'ProdE1', 'ProdE2', 'ProdE3', 'langloc']
info.contrast.unique() # 8 contrasts: ['H-E', 'NProd', 'S-N', 'SComp', 'SProd', 'VisEvSem', 'WComp', 'WProd']
info.imagetype.unique() # 2 image types: ['con.nii', 't.nii']]

info.subject.value_counts() # 28, 16, 8, 2
info.task.value_counts() 
# ProdE1     336
# ProdE3     168
# langloc     82
# MDloc       78
# ProdE2      44
info.contrast.value_counts() 
# SProd       106
# WProd       106
# NProd        84
# SComp        84
# VisEvSem     84
# WComp        84
# S-N          82
# H-E          78
info.imagetype.value_counts() # 354 each

# --- For visualisation only: replace nans with zero to display in fsleyes ---
# -- pull images for subject 756 (has 16 images) into visualisation folder with nans set to 0.
# mkdir /Volumes/diedrichsen_data$/data/FunctionalFusion/Language/raw/visualisation
# for i in /Volumes/diedrichsen_data$/data/FunctionalFusion/Language/raw/756_*; do fslmaths $i -nan ${i%%756*}/visualisation/${i##*raw/}; done


# --- For bids-standard: prepend sub- before subject id number ---
# for i in /Volumes/diedrichsen_data$/data/FunctionalFusion/Language/raw/*; do
# echo mv $i ${i%%raw/*}raw/sub-${i##*raw/};
# done

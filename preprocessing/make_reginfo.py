# Import /Users/callithrix/Desktop/reginfo_lang.tsv
import pandas as pd
import os
import numpy as np

# Repeat the content of the reginfo file for each run, replace 'run' with integer 1-8
# time_id and timepoint should be continuous
# Concate the content of the reginfo file for each run
# Save the concatenated content to a new file



# Load the reginfo file
reginfo = pd.read_csv('/Users/callithrix/Desktop/reginfo_lang.tsv', sep='\t')

# Loop through runs 1-8
for run in np.arange(1, 9):
    # Repeat the content of the reginfo file for each run, replace 'run' with integer 1-8
    reginfo_run = reginfo.copy()
    reginfo_run['run'] = run
    # Concate the content of the reginfo file for each run
    if run == 1:
        reginfo_concat = reginfo_run
    else:
        reginfo_concat = pd.concat([reginfo_concat, reginfo_run])

time_id = np.arange(1, reginfo_concat.shape[0] + 1)
timepoint = [f'T{i:04d}' for i in time_id]
reginfo_concat['time_id'] = time_id
reginfo_concat['timepoint'] = timepoint

# Save the concatenated content to a new file
reginfo_concat.to_csv('/Users/callithrix/Desktop/reginfo_lang_concat.tsv', sep='\t', index=False)

# Copy the file to the folder of one subject that will have the 'general info file'
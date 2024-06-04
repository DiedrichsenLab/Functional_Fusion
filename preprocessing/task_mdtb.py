import ProbabilisticParcellation.util as ut
import nibabel as nib
from pathlib import Path
import os.path as op
import subprocess
import random
from itertools import product
import pandas as pd
from datetime import datetime
import os
import shutil
import numpy as np

data_dir = f'{ut.base_dir}/../Cerebellum/super_cerebellum/'
fusion_dir = Path(f'{ut.base_dir}/MDTB/')
runs = np.arange(1, 33)
# get zeropadded numbers from 1 to 32
runs = [f'{run:02d}' for run in runs]
sessions = ["sc1", "sc2"]




if __name__ == "__main__":
    T = pd.read_csv(f'{fusion_dir}/participants.tsv', delimiter='\t')
    # --- Copy the raw runs into estimates ---
    for subject in T.iterrows():
        subject = subject[1].participant_id
        for session in sessions:
            task_dir = Path(f'{data_dir}/{session}/imaging_data/s{subject[-2:]}')
            for run in runs:
                task_file = f"{str(task_dir)}/rrun_{run}.nii"
                if op.exists(task_file):
                    subprocess.run(
                        ['cp', task_file, f"{fusion_dir}/derivatives/{subject}/{subject}_ses-{session}_run-{run}.nii"])
                print(f'Copied {run} for {subject} in {session}')
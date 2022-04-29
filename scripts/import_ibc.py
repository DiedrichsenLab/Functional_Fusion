#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to transfer IBC data from Drago to CBS

Author: Ana Luisa Pinho

Created: April 2022
Last update: April 2022
"""

import os
import glob
import subprocess
import pandas as pd

from pathlib import Path


# session_names = ['archi', 'hcp1', 'hcp2', 'rsvp']
# subjects_numbers = [1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]
session_names = ['archi']
subjects_numbers = [1]

drago = 'agrilopi@drago:/storage/store2/data/ibc/'
drago_derivatives = drago + 'derivatives/'
home = str(Path.home())
cbs = os.path.join(home, 'diedrichsen_data/data/FunctionalFusion/ibc')
cbs_derivatives = os.path.join(cbs, 'derivatives')

sess_map = 'ibc_sessions_map.tsv'
sess_struct = 'ibc_sessions_structure.tsv'
task_conditions = 'ibc_conditions.tsv'
dfm = pd.read_csv(open(sess_map), sep='\t', index_col=0)
dfs = pd.read_csv(open(sess_struct), sep='\t')
dfc = pd.read_csv(open(task_conditions), sep='\t')

subjects_list = ['sub-%02d' % s for s in subjects_numbers]


def import_estimates(sub, sname, df1, df2, df3):
    session = df1[df1[sub].values == sname].index.values[0]
    cbs_path = cbs_derivatives + '/' + sub + '/estimates/' + session
    if not os.path.exists(cbs_path):
        os.makedirs(cbs_path)
    else:
        for ng in glob.glob(cbs_path + '/*.nii.gz'):
            os.remove(ng)
    runs = df2[df2.session == sname].srun.values
    tasks = df2[df2.session == sname].task.values
    phasedir = df2[df2.session == sname].phase.values
    session_folder = drago_derivatives + sub + '/' + session
    for rn, tk, ph in zip(runs, tasks, phasedir):
        zfolder = session_folder + '/res_stats_' + tk + '_' + ph + \
            '/z_score_maps/'
        conditions = df3[df3.task == tk].condition.values
        regressors = df3[df3.task == tk].regressor.values
        for cond, reg in zip(conditions, regressors):
            drago_file = zfolder + cond + '.nii.gz'
            print(drago_file)
            p = subprocess.Popen(["scp", '-o BatchMode=yes', drago_file,
                                  cbs_path])
            p.wait()
            f = cbs_path + '/' + cond + '.nii.gz'
            print(f)
            ff = cbs_path + '/' + sub + '_' + session + '_run-' + \
                '%02d' % rn + '_reg-' + '%02d' % reg + '_zmap.nii.gz'
            print(ff)
            os.rename(f, ff)


if __name__ == "__main__":
    for subject in subjects_list:
        for session_name in session_names:
            import_estimates(subject, session_name, dfm, dfs, dfc)

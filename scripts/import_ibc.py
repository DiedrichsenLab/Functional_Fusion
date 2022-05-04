#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to transfer IBC data from Drago to CBS

Author: Ana Luisa Pinho

Created: April 2022
Last update: May 2022
"""

import os
import glob
import re
import subprocess
import numpy as np
import pandas as pd

from pathlib import Path


subjects_numbers = [1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]
# subjects_numbers = [1, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]

# session_names = ['archi', 'hcp1', 'hcp2', 'rsvp-language']
# session_names = ['archi', 'hcp1', 'hcp2', 'rsvp-language', 'mtt1', 'mtt2',
#                  'preference', 'tom', 'enumeration', 'self', 'clips4',
#                  'lyon1', 'lyon2']


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


def transfer_estimates(sub, sname, df1, df2, df3):
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
        if tk == 'RSVPLanguage':
            zfolder = session_folder + '/res_stats_' + tk + \
                '_%02d_' % (rn - 1) + ph + '/z_score_maps/'
        elif tk in ['MTTWE', 'MTTNS']:
            zfolder = session_folder + '/res_stats_' + tk + \
                '%d_' % rn + ph + '/z_score_maps/'
        elif sub == 'sub-11' and sname == 'preference' and rn == 6:
            tk = 'PreferenceFaces'
            zfolder = session_folder + '/res_stats_' + tk + '_' + ph + \
                '_run-01/z_score_maps/'
        elif sub == 'sub-11' and sname == 'preference' and rn == 7:
            zfolder = session_folder + '/res_stats_' + tk + '_' + ph + \
                '_run-02/z_score_maps/'
        else:
            zfolder = session_folder + '/res_stats_' + tk + '_' + ph + \
                '/z_score_maps/'
        if tk in ['VSTM' + '%d' % s for s in np.arange(1, 3)]:
            conditions = df3[df3.task == 'VSTM'].condition.values
            regressors = df3[df3.task == 'VSTM'].regressor.values
        elif tk in ['Self' + '%d' % s for s in np.arange(1, 5)]:
            conditions = df3[df3.task == 'Self'].condition.values
            regressors = df3[df3.task == 'Self'].regressor.values
        elif tk in ['WedgeClock', 'WedgeAnti']:
            conditions = df3[df3.task == 'Wedge'].condition.values
            regressors = df3[df3.task == 'Wedge'].regressor.values
        elif tk in ['ExpRing', 'ContRing']:
            conditions = df3[df3.task == 'Ring'].condition.values
            regressors = df3[df3.task == 'Ring'].regressor.values
        else:
            conditions = df3[df3.task == tk].condition.values
            regressors = df3[df3.task == tk].regressor.values
        for cond, reg in zip(conditions, regressors):
            drago_file = zfolder + cond + '.nii.gz'
            with subprocess.Popen(["scp", '-o BatchMode=yes', drago_file,
                                  cbs_path]) as p:
                p.wait()
            f = cbs_path + '/' + cond + '.nii.gz'
            print(f)
            ff = cbs_path + '/' + sub + '_' + session + '_run-' + \
                '%02d' % rn + '_reg-' + '%02d' % reg + '_zmap.nii.gz'
            print(ff)
            os.rename(f, ff)


def transfer_anat(pt):
    cbs_anatpath = cbs_derivatives + '/' + pt + '/anat/'
    if not os.path.exists(cbs_anatpath):
        os.makedirs(cbs_anatpath)
    else:
        for ng in glob.glob(cbs_anatpath + '/*.nii.gz'):
            os.remove(ng)
    drago_anatsess = drago_derivatives + pt + '/ses-00/anat/'
    drago_anatfile = drago_anatsess + pt + '_ses-00_T1w.nii.gz'
    w_drago_anatfile = drago_anatsess + 'w' + pt + '_ses-00_T1w.nii.gz'
    for afile in [drago_anatfile, w_drago_anatfile]:
        with subprocess.Popen(["scp", '-o BatchMode=yes', afile,
                               cbs_anatpath]) as a:
            a.wait()
        tag = re.match(
            drago_anatsess + '(.*)_ses-00_T1w.nii.gz', afile).groups()[0]
        t1 = cbs_anatpath + tag + '_ses-00_T1w.nii.gz'
        print(t1)
        new_t1 = cbs_anatpath + tag + '_T1w.nii.gz'
        print(new_t1)
        os.rename(t1, new_t1)


if __name__ == "__main__":
    for subject in subjects_list:
        # Import derivatives
        # for session_name in session_names:
        #     transfer_estimates(subject, session_name, dfm, dfs, dfc)
        # Import T1w images
        transfer_anat(subject)

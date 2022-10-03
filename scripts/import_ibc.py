#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to transfer IBC data from Drago to CBS

Author: Ana Luisa Pinho

Created: April 2022
Last update: October 2022
"""

import os
import glob
import re
import subprocess
import numpy as np
import pandas as pd

from pathlib import Path

import nibabel as nb
from nilearn.image import load_img, mean_img


# ############################# FUNCTIONS ##############################


def transfer_t1w(sub, source_derivatives, target_derivatives):
    target_anatpath = os.path.join(target_derivatives, sub, 'anat')
    if not os.path.exists(target_anatpath):
        os.makedirs(target_anatpath)
    else:
        for ng in glob.glob(target_anatpath + '/sub-*_T1w.nii*'):
            os.remove(ng)
    source_anatsess = os.path.join(source_derivatives, sub, 'ses-00/anat/')
    source_anatfile = os.path.join(source_anatsess, sub + '_ses-00_T1w.nii.gz')
    print(source_anatfile)

    with subprocess.Popen(["scp", '-o BatchMode=yes', source_anatfile,
                           target_anatpath]) as a:
        a.wait()
    t1 = os.path.join(target_anatpath, sub + '_ses-00_T1w.nii.gz')
    new_t1 = os.path.join(
        target_anatpath,
        sub + '_space-native_desc-resampled_T1w.nii.gz')
    print(t1)
    print(new_t1)
    os.rename(t1, new_t1)


def source_funcdata(sub, sname, original_sourcepath, destination_sourcepath,
                    df1, df2, data_type='bold.nii.gz'):
    """
    Data_type accepts 'bold.nii.gz' or 'events.tsv' strings as arguments.
    Default is 'bold.nii.gz'.
    """
    session = df1[df1[sub].values == sname].index.values[0]
    func_folder = os.path.join(original_sourcepath, sub, session, 'func')
    target_dir = os.path.join(destination_sourcepath, sub, 'func', session)
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    else:
        for ng in glob.glob(target_dir + '/*' + data_type):
            os.remove(ng)
    runs = df2[df2.session == sname].srun.values
    tasks = df2[df2.session == sname].task.values
    phasedir = df2[df2.session == sname].phase.values
    for i, (rn, tk, ph) in enumerate(zip(runs, tasks, phasedir)):
        if tk == 'RSVPLanguage':
            epi_fname = sub + '_' + session + '_task-' + tk + '_dir-' + ph + \
                '_run-%02d' % (rn - 1) + '_' + data_type
        elif tk == 'RestingState' and data_type == 'events.tsv':
            continue
        elif tk in ['MTTWE', 'MTTNS']:
            epi_fname = sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_run-%02d' % (rn - 2) + '_' + data_type
        elif sub == 'sub-11' and sname == 'preference' and rn == 6:
            tk = 'PreferenceFaces'
            epi_fname = sub + '_' + session + '_task-' + tk + '_dir-' + ph + \
                '_run-01_' + data_type
        elif sub == 'sub-11' and sname == 'preference' and rn == 7:
            epi_fname = sub + '_' + session + '_task-' + tk + '_dir-' + ph + \
                '_run-02_' + data_type
        elif tk in ['VSTM' + '%d' % s for s in np.arange(1, 3)]:
            if tk == 'VSTM1':
                epi_fname = sub + '_' + session + '_task-VSTM_dir-' + ph + \
                    '_run-01_' + data_type
            else:
                assert tk == 'VSTM2'
                epi_fname = sub + '_' + session + '_task-VSTM_dir-' + ph + \
                    '_run-02_' + data_type
        elif tk in ['Self' + '%d' % s for s in np.arange(1, 5)]:
            epi_fname = sub + '_' + session + '_task-Self_dir-' + \
                ph + '_run-%02d' % rn + '_' + data_type
        elif tk in ['MathLanguage1', 'SpatialNavigation']:
            epi_fname = sub + '_' + session + '_task-' + tk + '_dir-' + ph + \
                '_run-%02d' % rn + '_' + data_type
        elif tk == 'MathLanguage2':
            epi_fname = sub + '_' + session + '_task-' + tk + '_dir-' + ph + \
                '_run-%02d' % (rn - 4) + '_' + data_type
        else:
            epi_fname = sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_' + data_type
        epi_path = os.path.join(func_folder, epi_fname)
        with subprocess.Popen(["scp", '-o BatchMode=yes', epi_path,
                               target_dir]) as epi:
            epi.wait()
        target_epi = os.path.join(target_dir, epi_fname)
        new_epi_path = os.path.join(
            target_dir,
            sub + '_' + session + '_run-' + '%02d' % rn + '_' + data_type)
        print(target_epi)
        print(new_epi_path)
        os.rename(target_epi, new_epi_path)


def wepi(sub, sname, source_derivatives, target_derivatives, df1, df2,
         first_run_only = False):
    session = df1[df1[sub].values == sname].index.values[0]
    func_folder = os.path.join(source_derivatives, sub, session, 'func')
    target_dir = os.path.join(target_derivatives, sub, 'func', session)
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    else:
        for ng in glob.glob(target_dir + '/*bold.nii.gz'):
            os.remove(ng)
    runs = df2[df2.session == sname].srun.values
    tasks = df2[df2.session == sname].task.values
    phasedir = df2[df2.session == sname].phase.values
    for i, (rn, tk, ph) in enumerate(zip(runs, tasks, phasedir)):
        if first_run_only is True and i > 0:
            break
        if tk == 'RSVPLanguage':
            wepi_fname = 'wrdc' + sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_run-%02d' % (rn - 1) + '_bold.nii.gz'
        elif tk in ['MTTWE', 'MTTNS']:
            wepi_fname = 'wrdc' + sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_run-%02d' % (rn - 2) + '_bold.nii.gz'
        elif sub == 'sub-11' and sname == 'preference' and rn == 6:
            tk = 'PreferenceFaces'
            wepi_fname = 'wrdc' + sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_run-01_bold.nii.gz'
        elif sub == 'sub-11' and sname == 'preference' and rn == 7:
            wepi_fname = 'wrdc' + sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_run-02_bold.nii.gz'
        elif tk in ['VSTM' + '%d' % s for s in np.arange(1, 3)]:
            if tk == 'VSTM1':
                wepi_fname = 'wrdc' + sub + '_' + session + \
                    '_task-VSTM_dir-' + ph + '_run-01_bold.nii.gz'
            else:
                assert tk == 'VSTM2'
                wepi_fname = 'wrdc' + sub + '_' + session + \
                    '_task-VSTM_dir-' + ph + '_run-02_bold.nii.gz'
        elif tk in ['Self' + '%d' % s for s in np.arange(1, 5)]:
            wepi_fname = 'wrdc' + sub + '_' + session + '_task-Self_dir-' + \
                ph + '_run-%02d' % rn + '_bold.nii.gz'
        else:
            wepi_fname = 'wrdc' + sub + '_' + session + '_task-' + tk + \
                '_dir-' + ph + '_bold.nii.gz'
        wepi = os.path.join(func_folder, wepi_fname)
        with subprocess.Popen(["scp", '-o BatchMode=yes', wepi,
                               target_dir]) as epi:
            epi.wait()
        target_epi = os.path.join(target_dir, wepi_fname)
        new_target_epi = os.path.join(
            target_dir,
            sub + '_' + session + '_run-' + '%02d' % rn + '_bold.nii.gz')
        print(target_epi)
        print(new_target_epi)
        os.rename(target_epi, new_target_epi)


def compute_wmeanepi(sub, sname, target_derivatives, df1):
    session = df1[df1[sub].values == sname].index.values[0]
    sdir = os.path.join(target_derivatives, sub, 'func', session)
    wepis_paths = glob.glob(sdir + '/*_mean.nii.gz')
    for wepi_path in wepis_paths:
        wepi_fname = re.match(
            sdir + '/(.*)_bold.nii.gz', wepi_path).groups()[0]
        if wepi_fname == 'sub-14_ses-01_run-03':
            img = nb.load(os.path.join(
                sdir, 'sub-14_ses-01_run-03_bold.nii.gz'))
            X = img.get_fdata()
            Y = np.nanmean(X, axis=3)
            meanimg = nb.Nifti1Image(Y, img.affine)
            meanimg_path = os.path.join(
                sdir, 'sub-01_ses-03_run-01_mean.nii.gz')
            nb.save(meanimg, meanimg_path)
            print(meanimg_path)
        else:
            wepi = load_img(wepi_path)
            wmeanepi = mean_img(wepi)
            wmeanepi_fullpath = os.path.join(sdir, wepi_fname + '_mean.nii.gz')
            wmeanepi.to_filename(os.path.join(wmeanepi_fullpath))
            print(wmeanepi_fullpath)


def transfer_estimates(sub, sname, source_derivatives, target_derivatives,
                       df1, df2, df3):
    session = df1[df1[sub].values == sname].index.values[0]
    target_path = os.path.join(target_derivatives, sub, 'estimates', session)
    if not os.path.exists(target_path):
        os.makedirs(target_path)
    else:
        for ng in glob.glob(target_path + '/*.nii.gz'):
            os.remove(ng)
    runs = df2[df2.session == sname].srun.values
    tasks = df2[df2.session == sname].task.values
    phasedir = df2[df2.session == sname].phase.values
    session_folder = os.path.join(source_derivatives, sub, session)
    for rn, tk, ph in zip(runs, tasks, phasedir):
        if tk == 'RSVPLanguage':
            zfolder = os.path.join(
                session_folder,
                'res_stats_' + tk + '_%02d_' % (rn - 1) + ph,
                'z_score_maps')
        elif tk in ['MTTWE', 'MTTNS']:
            zfolder = os.path.join(
                session_folder,
                'res_stats_' + tk + '%d_' % (rn - 2) + ph,
                'z_score_maps')
        elif sub == 'sub-11' and sname == 'preference' and rn == 6:
            tk = 'PreferenceFaces'
            zfolder = os.path.join(
                session_folder,
                'res_stats_' + tk + '_' + ph + '_run-01',
                'z_score_maps')
        elif sub == 'sub-11' and sname == 'preference' and rn == 7:
            zfolder = os.path.join(
                session_folder,
                'res_stats_' + tk + '_' + ph + '_run-02',
                'z_score_maps')
        else:
            zfolder = os.path.join(
                session_folder,
                'res_stats_' + tk + '_' + ph,
                'z_score_maps')
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
            source_file = zfolder + cond + '.nii.gz'
            with subprocess.Popen(["scp", '-o BatchMode=yes', source_file,
                                  target_path]) as p:
                p.wait()
            f = os.path.join(target_path, cond + '.nii.gz')
            ff = os.path.join(
                target_path,
                sub + '_' + session + '_run-' + '%02d' % rn + '_reg-' + \
                '%02d' % reg + '_zmap.nii.gz')
            print(f)
            print(ff)
            os.rename(f, ff)


def generate_sessinfo(sub, sname, target_derivatives, df1, df2, df3):
    session = df1[df1[sub].values == sname].index.values[0]
    ifolder = os.path.join(target_derivatives, sub, 'estimates', session)
    subject_no = sub[4:]
    sess_no = session[4:]
    if not os.path.exists(ifolder):
        os.makedirs(ifolder)
    else:
        if glob.glob(ifolder + '/*.tsv'):
            for ng in glob.glob(ifolder + '/*.tsv'):
                os.remove(ng)
    run_numbers = df2[df2.session == sname].srun.values
    task_names = df2[df2.session == sname].task.values
    n_repetitions = df2[df2.session == sname].nrep.values
    sessinfo = np.empty((0, 7))
    for rnum, tname, n_rep in zip(run_numbers, task_names, n_repetitions):
        if sub == 'sub-11' and rnum == 6 and \
           tname == 'PreferencePaintings':
            tname = 'PreferenceFaces'
        condition_names = df3[df3.task == tname].condition.tolist()
        reg_numbers = df3[df3.task == tname].regressor.tolist()
        sub_rep = np.repeat(subject_no, len(condition_names)).tolist()
        sess_rep = np.repeat(sess_no, len(condition_names)).tolist()
        rnum_rep = np.repeat(rnum, len(condition_names)).tolist()
        tname_rep = np.repeat(tname, len(condition_names)).tolist()
        nrep_rep = np.repeat(n_rep, len(condition_names)).tolist()
        rstack = np.vstack((sub_rep, sess_rep, rnum_rep, tname_rep,
                            condition_names, reg_numbers, nrep_rep)).T
        sessinfo = np.vstack((sessinfo, rstack))
    dff = pd.DataFrame(sessinfo, columns = ['sn', 'sess', 'run', 'task_name',
                                            'cond_name', 'reg_num', 'n_rep'])
    dff_fname = sub + '_' + session + '_reginfo.tsv'
    dff_path = os.path.join(ifolder, dff_fname)
    dff.to_csv(dff_path, sep='\t', index=False)


def rename_sessions(sub, sname, folder_path, df1):
    session = df1[df1[sub].values == sname].index.values[0]
    old_folder = os.path.join(folder_path, sub, 'func', session)
    new_folder = os.path.join(folder_path, sub, 'func', sname)
    print(old_folder)
    print(new_folder)
    os.rename(old_folder, new_folder)


# ############################### INPUTS ###############################

# subjects_numbers = [1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]
subjects_numbers = [1, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]

# session_names = ['archi', 'hcp1', 'hcp2', 'rsvp-language']
session_names = ['mtt1', 'mtt2', 'preference', 'tom', 'enumeration', 'self',
                 'clips4', 'lyon1', 'lyon2', 'mathlang',
                 'spatial-navigation']


# ############################# PARAMETERS #############################

drago = 'agrilopi@drago:/storage/store2/data/ibc'
drago_sourcedata = os.path.join(drago, 'sourcedata')
drago_derivatives = os.path.join(drago, 'derivatives')

home = str(Path.home())
cbs = os.path.join(home, 'diedrichsen_data/data/ibc')
cbs_sourcedata = os.path.join(cbs, 'raw')
cbs_derivatives = os.path.join(cbs, 'derivatives')

sess_map = 'ibc_sessions_map.tsv'
sess_struct = 'ibc_sessions_structure.tsv'
task_conditions = 'ibc_conditions.tsv'
dfm = pd.read_csv(open(sess_map), sep='\t', index_col=0)
dfs = pd.read_csv(open(sess_struct), sep='\t')
dfc = pd.read_csv(open(task_conditions), sep='\t')

subjects_list = ['sub-%02d' % s for s in subjects_numbers]

# ############################# RUN ####################################

if __name__ == "__main__":
    for subject in subjects_list:
        ## Import T1w in the native space ##
        # transfer_t1w(subject, drago_derivatives, cbs_sourcedata)

        for session_name in session_names:
            ## Import raw EPI ##
            # source_funcdata(subject, session_name, drago_sourcedata,
            #                 cbs_sourcedata, dfm, dfs)

            ## Import normalized-EPI ##
            # wepi(subject, session_name, drago_derivatives, cbs_derivatives,
            #      dfm, dfs, first_run_only = True)

            ## Compute mean normalized-EPI ##
            # compute_wmeanepi(subject, session_name, cbs_derivatives, dfm)

            ## Import paradigm descriptors ##
            # source_funcdata(subject, session_name, drago_sourcedata,
            #                 cbs_sourcedata, dfm, dfs, data_type='events.tsv')

            ## Import derivatives ##
            # transfer_estimates(subject, session_name, drago_derivatives,
            #                    cbs_derivatives, dfm, dfs, dfc)

            ## Generate tsv files with session info ##
            # generate_sessinfo(subject, session_name, cbs_derivatives, dfm,
            #                   dfs, dfc)

            rename_sessions(subject, session_name, cbs_sourcedata, dfm)

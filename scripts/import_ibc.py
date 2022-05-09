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

from nilearn.image import load_img, mean_img


# ############################# FUNCTIONS ##############################

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
            ff = cbs_path + '/' + sub + '_' + session + '_run-' + \
                '%02d' % rn + '_reg-' + '%02d' % reg + '_zmap.nii.gz'
            print(f)
            print(ff)
            os.rename(f, ff)


def generate_sessinfo(individual, sesstag, derivatives, df1, df2, df3):
    session_number = df1[df1[individual].values == sesstag].index.values[0]
    ifolder = derivatives + '/' + individual + '/estimates/' + \
        session_number
    if not os.path.exists(ifolder):
        os.makedirs(ifolder)
    else:
        if not glob.glob(ifolder + '/*.tsv'):
            for ng in glob.glob(ifolder + '/*.tsv'):
                os.remove(ng)
    run_numbers = df2[df2.session == sesstag].srun.values
    task_names = df2[df2.session == sesstag].task.values
    sessinfo = np.empty((0, 4))
    for rnum, tname in zip(run_numbers, task_names):
        if individual == 'sub-11' and rnum == 6 and \
           tname == 'PreferencePaintings':
            tname = 'PreferenceFaces'
        condition_names = df3[df3.task == tname].condition.tolist()
        reg_numbers = df3[df3.task == tname].regressor.tolist()
        rnum_rep = np.repeat(rnum, len(condition_names)).tolist()
        tname_rep = np.repeat(tname, len(condition_names)).tolist()
        rstack = np.vstack((rnum_rep, tname_rep, condition_names,
                            reg_numbers)).T
        sessinfo = np.vstack((sessinfo, rstack))
    dff = pd.DataFrame(sessinfo, columns = ['run','task_name','cond_name',
                                            'reg_num'])
    dff_fname = individual + '_' + session_number + '_reginfo.tsv'
    dff_path = os.path.join(ifolder, dff_fname)
    dff.to_csv(dff_path, sep='\t', index=False)


def epi(sbj, sess_id, df1, df2, first_run_only = False):
    sess = df1[df1[sbj].values == sess_id].index.values[0]
    cbs_dir = cbs_derivatives + '/' + sbj + '/func/' + sess
    if not os.path.exists(cbs_dir):
        os.makedirs(cbs_dir)
    else:
        for ng in glob.glob(cbs_dir + '/*bold.nii.gz'):
            os.remove(ng)
    func_folder = drago_derivatives + sbj + '/' + sess + '/func/'
    runs = df2[df2.session == sess_id].srun.values
    tasks = df2[df2.session == sess_id].task.values
    phasedir = df2[df2.session == sess_id].phase.values
    for i, (rn, tk, ph) in enumerate(zip(runs, tasks, phasedir)):
        if first_run_only == True and i > 0:
            break
        if tk == 'RSVPLanguage':
            wepi_fname = 'wrdc' + sbj + '_' + sess + '_task-' + tk + \
                '_dir-' + ph + '_run-' + '%02d' % (rn - 1) + '_bold.nii.gz'
        elif tk in ['MTTWE', 'MTTNS']:
            wepi_fname = 'wrdc' + sbj + '_' + sess + '_task-' + tk + \
                '_dir-' + ph + '_run-' + '%02d' % rn + '_bold.nii.gz'
        elif tk in ['VSTM' + '%d' % s for s in np.arange(1, 3)]:
            wepi_fname = 'wrdc' + sbj + '_' + sess + '_task-VSTM_dir-' + \
                ph + '_run-' + '%02d' % rn + '_bold.nii.gz'
        elif tk in ['Self' + '%d' % s for s in np.arange(1, 5)]:
            wepi_fname = 'wrdc' + sbj + '_' + sess + '_task-Self_dir-' + \
                ph + '_run-' + '%02d' % rn + '_bold.nii.gz'
        else:
            wepi_fname = 'wrdc' + sbj + '_' + sess + '_task-' + tk + \
                '_dir-' + ph + '_bold.nii.gz'
        wepi = func_folder + wepi_fname
        with subprocess.Popen(["scp", '-o BatchMode=yes', wepi,
                               cbs_dir]) as epi:
            epi.wait()
        cbs_epi = cbs_dir + '/' + wepi_fname
        new_cbs_epi = cbs_dir + '/' + sbj + '_' + sess + '_run-' + \
            '%02d' % rn + '_bold.nii.gz'
        print(cbs_epi)
        print(new_cbs_epi)
        os.rename(cbs_epi, new_cbs_epi)


def compute_wmeanepi(subj, sess_id, derivatives, df1):
    sess = df1[df1[subj].values == sess_id].index.values[0]
    sdir = derivatives + '/' + subj + '/func/' + sess
    wepis_paths = glob.glob(sdir + '/*_bold.nii.gz')
    for wepi_path in wepis_paths:
        wepi_fname = re.match(
            sdir + '/(.*)_bold.nii.gz', wepi_path).groups()[0]
        if wepi_fname == 'sub-14_ses-01_run-03':
            continue
        wepi = load_img(wepi_path)
        wmeanepi = mean_img(wepi)
        wmeanepi_fullpath = os.path.join(sdir, wepi_fname + '_mean.nii.gz')
        wmeanepi.to_filename(os.path.join(wmeanepi_fullpath))
        print(wmeanepi_fullpath)


def transfer_anat(pt):
    cbs_anatpath = cbs_derivatives + '/' + pt + '/anat/'
    if not os.path.exists(cbs_anatpath):
        os.makedirs(cbs_anatpath)
    else:
        for ng in glob.glob(cbs_anatpath + '*_T1w.nii.gz'):
            os.remove(ng)
    drago_anatsess = drago_derivatives + pt + '/ses-00/anat/'
    drago_anatfile = drago_anatsess + pt + '_ses-00_T1w.nii.gz'
    w_drago_anatfile = drago_anatsess + 'w' + pt + '_ses-00_T1w.nii.gz'
    for afile in [drago_anatfile, w_drago_anatfile]:
        with subprocess.Popen(["scp", '-o BatchMode=yes', afile,
                               cbs_anatpath]) as a:
            a.wait()
        if afile == drago_anatfile:
            t1 = cbs_anatpath + pt + '_ses-00_T1w.nii.gz'
            new_t1 = cbs_anatpath + pt + '_space-native_T1w.nii.gz'
        else:
            assert afile == w_drago_anatfile
            t1 = cbs_anatpath + 'w' + pt + '_ses-00_T1w.nii.gz'
            new_t1 = cbs_anatpath + pt + '_T1w.nii.gz'
        print(t1)
        print(new_t1)
        os.rename(t1, new_t1)


def transfer_cmasks(pt, derivatives):
    cbs_anatpath = derivatives + '/' + pt + '/anat/'
    if not os.path.exists(cbs_anatpath):
        os.makedirs(cbs_anatpath)
    else:
        for ng in glob.glob(cbs_anatpath + 'mwc*.nii.gz'):
            os.remove(ng)
    drago_anatsess = drago_derivatives + pt + '/ses-00/anat/'
    drago_c1 = drago_anatsess + 'mwc1' + pt + '_ses-00_T1w.nii.gz'
    drago_c2 = drago_anatsess + 'mwc2' + pt + '_ses-00_T1w.nii.gz'
    for cfile in [drago_c1, drago_c2]:
        with subprocess.Popen(["scp", '-o BatchMode=yes', cfile,
                               cbs_anatpath]) as c:
            c.wait()
        if cfile == drago_c1:
            cmask = cbs_anatpath + 'mwc1' + pt + '_ses-00_T1w.nii.gz'
            new_cmask = cbs_anatpath + pt + '_mask-c1_T1w.nii.gz'
        else:
            assert cfile == drago_c2
            cmask = cbs_anatpath + 'mwc2' + pt + '_ses-00_T1w.nii.gz'
            new_cmask = cbs_anatpath + pt + '_mask-c2_T1w.nii.gz'
        print(cmask)
        print(new_cmask)
        os.rename(cmask, new_cmask)


def transfer_meshes(participant):
    cbs_meshfolder = cbs_derivatives + '/' + participant + '/anat/'
    if not os.path.exists(cbs_meshfolder):
        os.makedirs(cbs_meshfolder)
    else:
        for ng in glob.glob(cbs_meshfolder + '*.surf'):
            os.remove(ng)
    drago_meshfolder = drago_derivatives + participant + '/ses-00/anat/' + \
        participant + '/surf/'
    hemispheres = ['lh', 'rh']
    meshes = ['orig', 'pial', 'sulc', 'white']
    for hemi in hemispheres:
        for mesh in meshes:
            drago_meshfile = drago_meshfolder + hemi + '.' + mesh
            with subprocess.Popen(["scp", '-o BatchMode=yes', drago_meshfile,
                                   cbs_meshfolder]) as m:
                m.wait()
            cbs_meshfile = cbs_meshfolder + hemi + '.' + mesh
            if hemi == 'lh':
                new_cbs_meshfile = cbs_meshfolder + participant + \
                    '_hemi-L_' + mesh + '.surf'
            else:
                assert hemi == 'rh'
                new_cbs_meshfile = cbs_meshfolder + participant + \
                    '_hemi-R_' + mesh + '.surf'
            print(cbs_meshfile)
            print(new_cbs_meshfile)
            os.rename(cbs_meshfile, new_cbs_meshfile)


# ############################### INPUTS ###############################

# subjects_numbers = [1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]
# subjects_numbers = [1, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]
subjects_numbers = [4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]
# subjects_numbers = [1]

# session_names = ['archi', 'hcp1', 'hcp2', 'rsvp-language']
session_names = ['mtt1', 'mtt2', 'preference', 'tom', 'enumeration', 'self',
                 'clips4', 'lyon1', 'lyon2']
# session_names = ['self', 'clips4', 'lyon1', 'lyon2']
# session_names = ['enumeration']


# ############################# PARAMETERS #############################

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


# ############################# RUN ####################################

if __name__ == "__main__":
    for subject in subjects_list:
        # ## Import T1w images ##
        # transfer_anat(subject)

        # ## Import cmasks ##
        # transfer_cmasks(subject, cbs_derivatives)

        # ## Import Freesurfer meshes ##
        # transfer_meshes(subject)

        for session_name in session_names:
            ## Import mean EPI ##
            epi(subject, session_name, dfm, dfs, first_run_only = True)

            # ## Import derivatives ##
            # transfer_estimates(subject, session_name, dfm, dfs, dfc)

            # ## Generate tsv files with session info ##
            # generate_sessinfo(subject, session_name, cbs_derivatives, dfm,
            #                   dfs, dfc)

            # ## Compute mean EPI ##
            compute_wmeanepi(subject, session_name, cbs_derivatives, dfm)

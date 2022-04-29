import os
import glob
import re
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

sessions_map = 'ibc_sessions.csv'
df = pd.read_csv(open(sessions_map), sep=',', index_col=0)

subjects_list = ['sub-%02d' % s for s in subjects_numbers]


def import_func_derivatives(sub, sname):
    session = df[df[sub].values == sname].index.values[0]
    drago_files = drago_derivatives + sub + '/' + session + \
        '/*ffx/stat_maps/*.nii.gz'
    target_path = cbs_derivatives + '/' + sub + '/estimates/' + session
    if not os.path.exists(target_path):
        os.makedirs(target_path)
    else:
        for ng in glob.glob(target_path + '/*.nii.gz'):
            os.remove(ng)
    p = subprocess.Popen(["scp", '-o BatchMode=yes', drago_files, target_path])
    p.wait()
    original_files = target_path + '/*.nii.gz'
    for f in glob.glob(original_files):
        original_fname = re.match(target_path + '/(.*).nii.gz', f).groups()[0]
        final_fname = sub + '_' + session + '_reg-' + original_fname + \
            '_zmaps.nii.gz'
        ff = target_path + '/' + final_fname
        os.rename(f, ff)


if __name__ == "__main__":
    for subject in subjects_list:
        for session_name in session_names:
            import_func_derivatives(subject, session_name)

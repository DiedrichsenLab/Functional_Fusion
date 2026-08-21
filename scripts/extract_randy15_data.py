# Script for preparing all data from Randy 15 subject dataset
import pandas as pd
import shutil, time, subprocess
from pathlib import Path
import nitools as nt
import numpy as np
import Functional_Fusion.atlas_map as am
from Functional_Fusion.dataset import DataSetRANDY15
import Functional_Fusion.dataset as ds
import nibabel as nb
import SUITPy as suit
import os, sys, re
import scipy.io as spio
import matplotlib.pyplot as plt
import Functional_Fusion.connectivity as conn
import Functional_Fusion.util as ut

ERIS_DIR = '/home/dzhi/eris_mount'

if not Path(ERIS_DIR).exists():
    ERIS_DIR = '/data/tge'
if not Path(ERIS_DIR).exists():
    raise (NameError('Could not find ERIStwo datashare'))

base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = ERIS_DIR + '/Tian/UKBB_full/imaging'
if not Path(base_dir).exists():
    base_dir = '/srv/diedrichsen/data/FunctionalFusion'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion'
if not Path(base_dir).exists():
    print('diedrichsen data server not mounted')

randy15_dir = ERIS_DIR + '/dzhi/projects/RANDY15'
atlas_dir = base_dir + '/Atlases'
hem_name = ['cortex_left', 'cortex_right']
data_dir = ERIS_DIR + '/Tian/RANDY15'

def extract_randy15_timeseries(ses_id='ses-rest1', type='Tseries', atlas='MNISymC3'):
    hcp_dataset = DataSetRANDY15(data_dir)
    hcp_dataset.extract_all(ses_id, type, atlas)


def make_info(type='Tseries', ses_id='ses-rest1'):
    """Adding an extra column 'run_id' to the tsv file

    Args:
        type: file type
        ses_id: session id

    Returns:
        Write in the modified tsv file
    """
    hcp_dataset = DataSetRANDY15(randy15_dir)

    T = pd.read_csv(hcp_dataset.base_dir + '/participants.tsv', sep='\t')
    for p, participant_id in enumerate(T.participant_id):
        # Make info
        dest_dir = hcp_dataset.base_dir + \
            f'/derivatives/{participant_id}/data/'
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        info = pd.read_csv(dest_dir + f'{participant_id}_{ses_id}_info-{type}.tsv',
                           sep='\t')

        info['run_id'] = info['run'].copy()
        if ses_id == 'ses-rest2':
            info['run_id'] = info['run_id'] + 2

        info.to_csv(
            dest_dir + f'{participant_id}_{ses_id}_info-{type}.tsv', sep='\t', index=False)



def get_hcp_fs32k_rsfc(type='Net69Run', space='MNISymC3', ses_id='ses-rest1',
                       subj_list=None, smooth=None, kernel=None, thres=None, keeptop=False):
    # Load dataset
    if subj_list is None:
        dset = DataSetRANDY15(randy15_dir)
    else:
        dset = DataSetRANDY15(randy15_dir, subj_id_file=subj_list)
    T = dset.get_participants()
    
    # Load the cortical networks
    target, type = re.findall(r'[A-Z]+[a-z0-9]*', type)      
    res = target[3:]
    
    if target[:3] == 'Net':
        net = nb.load(randy15_dir + '/derivatives/group' +
                      f'/{target}_space-fs32k.dscalar.nii')
    elif target[:3] == 'Ico':
        net = [atlas_dir + f'/tpl-fs32k/Icosahedron{res}.L.label.gii',
            atlas_dir + f'/tpl-fs32k/Icosahedron{res}.R.label.gii']
    elif target[:3] == 'Fus':
        net = nb.load(randy15_dir + '/derivatives/group' +
                      f'/{target}_space-fs32k.pscalar.nii')
    atlas, _ = am.get_atlas(space)
        
    for i, s in enumerate(T.participant_id):
        dest_dir = dset.data_dir.format(s)
        if smooth is not None:
            file_name = f'{dest_dir}/{s}_space-{space}_{ses_id}'\
                        f'_{target+type}_desc-sm{smooth}'\
                        f'{kernel if kernel is not None else ""}'
        else:
            file_name = f'{dest_dir}/{s}_space-{space}_{ses_id}'\
                        f'_{target+type}'
        
        file_name = file_name + '_binarized' if thres is not None else file_name
        file_name = file_name + '_kt' if keeptop else file_name
        file_name += '.dscalar.nii'
        # if os.path.exists(file_name):
        #     print(f'Already extracted {file_name}, skipping...')
        # else:
        print(f"-- Start extracting rsFC {s} {ses_id} {target+type} "
                f"smooth={smooth} kernel={kernel} binarize={thres} keeptop={keeptop} --")
        try:
            start = time.perf_counter()
            # Get the subject's data
            if smooth is None:
                data_cortex_subj, info = dset.get_data(
                    space='fs32k', ses_id=ses_id, type='Tseries', subj=[i])
                data_cortex_subj = data_cortex_subj.squeeze()
            else:
                _, info = dset.get_data(space='fs32k', ses_id=ses_id,
                                        type='Tseries', subj=[i])
                tmp_subjlist = pd.DataFrame({'participant_id': [s]})
                tmp_subjlist.to_csv(f'{dest_dir}/{s}.tsv', sep='\t', index=False)
                data_cortex_subj = smooth_hcp_fs32k(f'{dest_dir}/{s}.tsv', ses_id=ses_id, 
                                                    type='Tseries', smooth=smooth, 
                                                    kernel=kernel, return_data_only=True)
                os.remove(f'{dest_dir}/{s}.tsv')

            if target[:3] == 'Net' or target[:3] == 'Fus':
                names = [f'Network_{i}' for i in range(1, int(res)+1)]
                if target[:3] == 'Fus':
                    icos = [atlas_dir + f'/tpl-fs32k/Icosahedron1002.L.label.gii',
                        atlas_dir + f'/tpl-fs32k/Icosahedron1002.R.label.gii']
                    # Average the subject's cortical data within each icosahedron 
                    # if using the Fusion connectivity maps (which are given at 
                    # icosahedron1002 resolution)
                    data_cortex_subj, _ = conn.average_within_Icos(
                        icos, data_cortex_subj)
                    names = net.header.get_axis(0).name.tolist()
                
                # Regress each network into the fs32k cortical 
                # data to get a run-specific network timecourse
                network_timecourse = conn.regress_networks(
                    net.get_fdata(), data_cortex_subj)
                        
            elif target[:3] == 'Ico':
                # Average within each parcel
                network_timecourse, names = conn.average_within_Icos(
                    net, data_cortex_subj)
                network_timecourse = network_timecourse.T
                sides = np.repeat(['L', 'R'], len(names) / 2)
                names = [f'Ico_{sides[i]}{name}' for i,name in enumerate(names)]  
                        
            elif target[:3] == 'ICA':
                assert smooth is None, "ICA components were calculated on unsmoothed data, \
                                cannot apply functional connectvity profile on smoothed data!"
                ses_dic = {'ses-rest1': 0, 'ses-rest2': 1}
                ica_dir = dset.base_dir + f'/derivatives/group/node_timeseries/3T_HCP1200_MSMAll_d{res}_ts2'
                network_timecourse = np.loadtxt(ica_dir + f'/{s}.txt').T
                network_timecourse = np.hsplit(network_timecourse,2)[ses_dic[ses_id]]
                names = [f'Network_{i}' for i in range(1, int(res)+1)]

            # Calculate the connectivity fingerprint
            coef = conn.connectivity_fingerprint(data_cortex_subj, network_timecourse, info,
                                                type, threshold=thres, keeptop=keeptop)
            # Make info
            runs = np.repeat([info.run.unique()], len(names))
            net_id = np.tile(np.arange(len(names)),
                            int(coef.shape[0] / len(names))) + 1
            info = pd.DataFrame({'sn': [s] * coef.shape[0],
                                'sess': [ses_id] * coef.shape[0],
                                'run': runs,
                                'half': 2 - (runs < runs[-1]),
                                'net_id': net_id,
                                'names': names * int(coef.shape[0] / len(names))})

            # Save the data
            C = atlas.data_to_cifti(coef, info.names)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)

            nb.save(C,  file_name)
            info.to_csv(f'{dest_dir}/{s}_{ses_id}_info-{target+type}.tsv', 
                        sep='\t', index=False)
            
            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - time {elapse}")
        except:
            print(f"   Incomplete - {s} {ses_id} {target+type}!")

def smooth_hcp_fs32k(bulk, ses_id='ses-s1', type='Tseries', smooth=1, kernel=None,
                     return_data_only=True):
    hcp_dataset = ds.DataSetUkbResting(randy15_dir)
    T = pd.read_csv(bulk, sep='\t')

    # get the surfaces for smoothing
    surf_L = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-L_midthickness.surf.gii'
    surf_R = atlas_dir + f'/tpl-fs32k/tpl-fs32k_hemi-R_midthickness.surf.gii'

    for s in T.participant_id.drop_duplicates():
        # Make smoothed file name
        dest_dir = hcp_dataset.data_dir.format(s)
        cifti_out = dest_dir + f'/{s}_space-fs32k_{ses_id}_{type}_desc' \
                    f'-sm{smooth}{kernel if kernel is not None else ""}.dscalar.nii'
        
        if os.path.exists(cifti_out):
            print(f"Already smoothed for {s} fs32k {ses_id} {smooth}")
            if return_data_only:
                data = nb.load(cifti_out).get_fdata()
                return data
        else:
            print(f'- Smoothing data for {s} fs32k {ses_id} in {smooth}mm ...')
            contain_nan = False
            # Load the unsmoothed data
            input_file = hcp_dataset.data_dir.format(s) \
                        + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
            C = nb.load(input_file)

            # fill nan with zeros if unsmoothed data contains any
            if np.isnan(C.get_fdata()).any():
                contain_nan = True
                mask = np.isnan(C.get_fdata())
                C = nb.Cifti2Image(dataobj=np.nan_to_num(C.get_fdata()), header=C.header)
                nb.save(C, f'{s}_tmp.dscalar.nii')
                input_file = f'{s}_tmp.dscalar.nii'

            # Write in smoothed surface data (filled with 0)
            smooth_cmd = f"wb_command -cifti-smoothing {input_file} " \
                        f"{smooth} {smooth} COLUMN {cifti_out} " \
                        f"{f'-{kernel} ' if kernel is not None else ''}" \
                        f"-left-surface {surf_L} -right-surface {surf_R} " \
                        f"-fix-zeros-surface"
            subprocess.run(smooth_cmd, shell=True)
            
            C = nb.load(cifti_out)
            data = C.get_fdata()
            
            if contain_nan:
                os.remove(f'{s}_tmp.dscalar.nii')
                # Replace 0s back to NaN (we don't want the 0s impact model learning)
                data[mask] = np.nan
                C = nb.Cifti2Image(dataobj=data, header=C.header)
                nb.save(C, cifti_out)
            
            if return_data_only:
                os.remove(cifti_out)
                return data

def binarize_rsfc_fs32k(bulk, ses_id='ses-s1', type='Tseries', smooth=None, kernel=None,
                        thres=None):
    randy_dataset = ds.DataSetRANDY15(data_dir)
    T = pd.read_csv(bulk, sep='\t')

    for s in T.participant_id:
        # Make smoothed file name
        dest_dir = randy_dataset.data_dir.format(s)
        
        if smooth is not None:
            file_name = f'{dest_dir}/{s}_space-fs32k_{ses_id}'\
                        f'_{type}_desc-sm{smooth}'\
                        f'{kernel if kernel is not None else ""}'
        else:
            file_name = f'{dest_dir}/{s}_space-fs32k_{ses_id}'\
                        f'_{type}'

        print(f'- Banarizing data for {s} fs32k {ses_id} {smooth} ...')
        try:
            # Load the data to be binarized
            C = nb.load(file_name + '.dscalar.nii')
            coef = conn.binarize_top_percent(C.get_fdata(), percent=thres, keep_top=False)
            C = nb.Cifti2Image(dataobj=coef, header=C.header)
            nb.save(C, file_name + f'_binarized-{thres}.dscalar.nii')
            print(f'- Done.')
        except:
            print(f'Missing file subject {s} for {ses_id}, skip..')

def fsaverage6_to_fs32k(src_file, hemi, outfile, type='label', return_data_only=False):
    mesh_dir = f'{atlas_dir}/standard_mesh_atlases/resample_fsaverage'
    fs6 = mesh_dir + f'/fsaverage6_std_sphere.{hemi}.41k_fsavg_{hemi}.surf.gii'
    fs32k = mesh_dir + f'/fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii'

    fs6_mid = mesh_dir + f'/fsaverage6.{hemi}.midthickness_va_avg.41k_fsavg_{hemi}.shape.gii'
    fs32k_mid = mesh_dir + f'/fs_LR.{hemi}.midthickness_va_avg.32k_fs_LR.shape.gii'

    cmd = f'wb_command -{type}-resample {src_file} {fs6} {fs32k} ' \
          f'ADAP_BARY_AREA {outfile} -area-metrics {fs6_mid} {fs32k_mid}'
    subprocess.run(cmd, shell=True)

    if return_data_only:
        fs32k_data = nb.load(outfile)
        data = np.stack([dat.data for dat in fs32k_data.darrays])
        os.remove(outfile)
        return data

def fs32k_to_fsaverage6(src_file, hemi, outfile, type='label', return_data_only=False):
    mesh_dir = f'{atlas_dir}/standard_mesh_atlases/resample_fsaverage'
    fs6 = mesh_dir + f'/fsaverage6_std_sphere.{hemi}.41k_fsavg_{hemi}.surf.gii'
    fs32k = mesh_dir + f'/fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii'

    fs6_mid = mesh_dir + f'/fsaverage6.{hemi}.midthickness_va_avg.41k_fsavg_{hemi}.shape.gii'
    fs32k_mid = mesh_dir + f'/fs_LR.{hemi}.midthickness_va_avg.32k_fs_LR.shape.gii'

    cmd = f'wb_command -{type}-resample {src_file} {fs32k} {fs6} ' \
          f'ADAP_BARY_AREA {outfile} -area-metrics {fs32k_mid} {fs6_mid}'
    subprocess.run(cmd, shell=True)

    if return_data_only:
        fs6_data = nb.load(outfile)
        data = np.stack([dat.data for dat in fs6_data.darrays])
        os.remove(outfile)
        return data
    
def get_DU15_parcellation(file_name='DU15NET_Prior', atlas_space='fs32k'):
    atlas, _ = am.get_atlas(atlas_space)
    DU15_dir = '/data/tge/dzhi/workspace/DU15NET'
    file = nb.load(DU15_dir + f'/HCP/fsLR_32k/{file_name}_fsLR_32k.dlabel.nii')
    DU15 = atlas.cifti_to_data(file).reshape(-1)

    info = pd.read_csv(DU15_dir + '/DU15NET_ColorLUT.csv')
    network_names = list(info['Abbreviation'])
    colors = info[["R","G","B","A"]].to_numpy().astype(float)
    colors[:, :3] = colors[:, :3] / 255

    return DU15, network_names, colors

def convert_fs6_indivPar_to_fs32k(src_file):
    atlas, _ = am.get_atlas('fs32k')
    atlas.calculate_symmetry()
    data = nb.load(src_file).get_fdata()

    lh_data = data.T[:40962, :]
    rh_data = data.T[40962:, :]
            
    print(f'   Converting data to fsaverage6 func.gii')
    lh_gii = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft', 
                                column_names=[f'time_{i+1}' for i in range(lh_data.shape[1])])
    rh_gii = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight', 
                                column_names=[f'time_{i+1}' for i in range(rh_data.shape[1])])
    tmp_file_L = Path(src_file).parent / 'tmp.L.func.gii'
    tmp_file_R = Path(src_file).parent / 'tmp.R.func.gii'
    nb.save(lh_gii, tmp_file_L)
    nb.save(rh_gii, tmp_file_R)

    print(f'   Mapping data from fsaverage6 space to fs32k space ...')
    lh_fs32k_gii = Path(src_file).parent / 'tmp_fs32k.L.func.gii'
    rh_fs32k_gii = Path(src_file).parent / 'tmp_fs32k.R.func.gii'
    imgL = fsaverage6_to_fs32k(tmp_file_L, 'L', lh_fs32k_gii, type='metric', return_data_only=True)
    imgR = fsaverage6_to_fs32k(tmp_file_R, 'R', rh_fs32k_gii, type='metric', return_data_only=True)
    print(f'   Remove tmp files, Done.')
    os.remove(tmp_file_L)
    os.remove(tmp_file_R)

    print(f'   Combine Left and Right hemisphere into a single cifti ...')
    dat_L = imgL.darrays[0].data[atlas.vertex_mask[0]]
    dat_R = imgR.darrays[0].data[atlas.vertex_mask[1]]

    return np.concat([dat_L, dat_R])

def extract_timeseries(ses_id='ses-rest', type='Tseries', atlas='MNISymC3'):
    randy15_dataset = DataSetRANDY15(data_dir)
    ses = ses_id.split('-')[1].upper()

    # Extract the data for each participant
    T = randy15_dataset.get_participants()
    for row in T.itertuples(index=False):
        sub_id = row.participant_num
        s = row.participant_id
        ses_dir = f'{randy15_dir}/{sub_id}/{ses}'

        print(f'Extract {s}')
        file_list_L = [f for f in os.listdir(ses_dir) if f.startswith("lh.")]
        file_list_R = [f for f in os.listdir(ses_dir) if f.startswith("rh.")]
        assert len(file_list_L) == len(file_list_R), \
            "file number of left / right doesn't match!"
        
        for file_L, file_R in zip(file_list_L, file_list_R):
            assert file_L.split('.')[1] == file_R.split('.')[1], \
                    "run id for L / R hemisphere doesn't match!"
            run_id = int(re.search(r'\d+$', file_L.split('.')[1]).group())
            print(f'-- Start processing {s} {ses_id} run {run_id} ...')

            # Option 2 - us wb command resample
            lh_data = nb.load(f'{ses_dir}/{file_L}').get_fdata()
            lh_data = lh_data.reshape(-1, lh_data.shape[-1], order='F')
            rh_data = nb.load(f'{ses_dir}/{file_R}').get_fdata()
            rh_data = rh_data.reshape(-1, rh_data.shape[-1], order='F')

            if np.equal(lh_data, rh_data).all():
                print(f'{sub_id} ses {file_L} and {file_R} data identical!')
                break
            
            print(f'   Converting data to fsaverage6 func.gii')
            lh_gii = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft', 
                                        column_names=[f'time_{i+1}' for i in range(lh_data.shape[1])])
            rh_gii = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight', 
                                        column_names=[f'time_{i+1}' for i in range(rh_data.shape[1])])
            nb.save(lh_gii, ses_dir + f'/{sub_id}_tmp_{file_L}.L.func.gii')
            nb.save(rh_gii, ses_dir + f'/{sub_id}_tmp_{file_R}.R.func.gii')

            print(f'   Mapping data from fsaverage6 space to fs32k space ...')
            lh_fs32k_gii = ses_dir + f'/P1_REST{run_id}_wb.L.func.gii'
            rh_fs32k_gii = ses_dir + f'/P1_REST{run_id}_wb.R.func.gii'
            fsaverage6_to_fs32k(ses_dir + f'/{sub_id}_tmp_{file_L}.L.func.gii', 'L', lh_fs32k_gii)
            fsaverage6_to_fs32k(ses_dir + f'/{sub_id}_tmp_{file_R}.R.func.gii', 'R', rh_fs32k_gii)
            print(f'   Remove tmp files, Done.')
            os.remove(ses_dir + f'/{sub_id}_tmp_{file_L}.L.func.gii')
            os.remove(ses_dir + f'/{sub_id}_tmp_{file_R}.R.func.gii')

            print(f'   Combine Left and Right hemisphere into a single cifti ...')
            dest_dir = randy15_dataset.estimates_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)

            try:
                cmd = (f'wb_command -cifti-create-dense-timeseries '
                            f'{dest_dir}/{s}_space-fs32k_{ses_id}_run{run_id}.dtseries.nii '
                            f'-left-metric {lh_fs32k_gii} -right-metric {rh_fs32k_gii} '
                            f'-timestep 1.0 -timestart 0')
                subprocess.run(cmd, shell=True, check=True)
                os.remove(lh_fs32k_gii)
                os.remove(rh_fs32k_gii)

            except:
                print(f"Failed to combine in single CIFTI for {s} {ses_id} run {run_id}!")
        
            print('-- Done!')

def import_rs_timeseries(run_id=1, type='Tseries', atlas='MNISymC3'):
    myatlas, _ = am.get_atlas('fs32k')
    myatlas.calculate_symmetry()
    randy15_dataset = DataSetRANDY15(data_dir)

    T = pd.read_csv(data_dir + '/participants.tsv', delimiter='\t')
    for s in T.participant_id:
        dest_dir = randy15_dataset.data_dir.format(s)
        # if os.path.exists(dest_dir + f'/{s}_space-fs32k_{ses_id}_Tseries.dscalar.nii'):
        #     print(f"Already imported subject {s} {ses_id} Tseries, skipping...")
        # else:
        source_dir = randy15_dataset.estimates_dir.format(s)

        try:
            print(f"-- Start importing subject {s} {run_id} Tseries --")
            
            start = time.perf_counter()
            info = pd.DataFrame()
            data_file = source_dir + f'/{s}_space-fs32k_ses-rest_run{run_id}.dtseries.nii'
            data = myatlas.cifti_to_data(data_file)

            num_tpoint = [f'T{i+1:04}' for i in range(data.shape[0])]
            info = pd.DataFrame({'sn': [s] * data.shape[0],
                                        'run': [run_id] * data.shape[0],
                                        'timepoint': num_tpoint,
                                        'task': 'rest',
                                        'time_id': [i+1 for i in range(data.shape[0])],
                                        'names': num_tpoint})

            C = myatlas.data_to_cifti(data, num_tpoint)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)
            
            nb.save(C, dest_dir +
                        f'/{s}_space-fs32k_ses-rest{run_id}_Tseries.dscalar.nii')
            info.to_csv(dest_dir + f'/{s}_ses-rest{run_id}_Tseries.tsv', 
                        sep='\t', index=False)
            
            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - time {elapse}")
        except:
            print(f"-- Failed importing subject {s} {run_id} Tseries --")


def concatenate_rs_timeseries(run_id_list, space='fs32k', start_point=0, duration=600, tr=1):
    myatlas, _ = am.get_atlas('fs32k')
    myatlas.calculate_symmetry()
    randy15_dataset = DataSetRANDY15(data_dir)

    T = pd.read_csv(data_dir + '/participants.tsv', delimiter='\t')
    data_all, info_all = [],[]
    for s in T.participant_id:
        source_dir = randy15_dataset.data_dir.format(s)

        data = []
        for run_id in run_id_list:
            print(f"-- Concatenating subject {s} run{run_id} Tseries --")
            start = time.perf_counter()
            data_file = source_dir + f'/{s}_space-{space}_ses-rest{run_id}_Tseries.dscalar.nii'
            this_data = myatlas.cifti_to_data(data_file)

            finish = time.perf_counter()
            elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
            print(f"   Done - time {elapse}")
            data.append(ut.zstandarize_ts(this_data))

        data = np.vstack(data)
        end_point = int(start_point + duration / tr)
        data = data[start_point:end_point, :]
        this_num_tpoint = [f'T{i+1:04}' for i in range(end_point-start_point)]

        ## Make new data info
        info = pd.DataFrame({'sn': [s] * data.shape[0],
                            'run': [1] * data.shape[0],
                            'timepoint': this_num_tpoint,
                            'task': ['rest'] * data.shape[0],
                            'time_id': [i+1 for i in range(data.shape[0])],
                            'names': this_num_tpoint})
        
        C = myatlas.data_to_cifti(np.vstack(data), this_num_tpoint)
        nb.save(C, source_dir +
                f'/{s}_space-{space}_ses-rest{duration}s2_Tseries.dscalar.nii')
        info.to_csv(source_dir + f'/{s}_ses-rest{duration}s2_Tseries.tsv',
                    sep='\t', index=False)

        info_all.append(info)
        data_all.append(data)

    assert all(info_all[0].shape == df.shape for df in info_all[1:]), \
        "Not all subjects have same length Tseries!"
    info_com = info_all[0]
    data_all = np.stack(data_all)

    return data_all, info_com


def concatenate_rest_fc(rest_runs=None, space='fs32k', type='Ico642Run',
                        ses_out='ses-rest', binarized='0.1',
                        subj=None, skip_missing=True):
    """Concatenate run-wise resting-state FC files into one ses-rest file.

    The input files are expected in each subject's data directory as:
        sub-*_space-{space}_ses-rest{run}_{type}_binarized-{binarized}.dscalar.nii
        sub-*_ses-rest{run}_{type}.tsv

    The output files are:
        sub-*_space-{space}_{ses_out}_{type}_binarized-{binarized}.dscalar.nii
        sub-*_{ses_out}_{type}.tsv
    """
    randy15_dataset = DataSetRANDY15(data_dir)
    T = randy15_dataset.get_participants()

    if subj is None:
        subjects = T.participant_id
    elif isinstance(subj, str):
        subjects = [subj]
    else:
        subjects = subj

    summary = []
    for s in subjects:
        dest_dir = randy15_dataset.data_dir.format(s)
        if not os.path.exists(dest_dir):
            print(f'Missing data directory for {s}, skip ...')
            continue

        if rest_runs is None:
            pattern = re.compile(
                rf'{re.escape(s)}_space-{re.escape(space)}_ses-rest(\d+)'
                rf'_{re.escape(type)}'
                rf'_binarized-{re.escape(str(binarized))}\.dscalar\.nii$'
            )
            found_runs = []
            for file in os.listdir(dest_dir):
                match = pattern.match(file)
                if match:
                    found_runs.append(int(match.group(1)))
            run_list = sorted(found_runs)
        else:
            run_list = list(rest_runs)

        if not run_list:
            print(f'No resting FC files found for {s}, skip ...')
            continue

        data_list, info_list, row_names = [], [], []
        brain_axis = None
        included_runs = []
        print(f'Concatenating resting FC for {s}: runs {run_list}')

        for run in run_list:
            run = int(run)
            ses_id = f'ses-rest{run}'
            cifti_file = (f'{dest_dir}/{s}_space-{space}_{ses_id}_{type}'
                          f'_binarized-{binarized}.dscalar.nii')
            tsv_file = f'{dest_dir}/{s}_{ses_id}_{type}.tsv'

            if not os.path.exists(cifti_file) or not os.path.exists(tsv_file):
                msg = f'Missing resting FC files for {s} {ses_id}'
                if skip_missing:
                    print(f'   {msg}, skip run ...')
                    continue
                raise FileNotFoundError(msg)

            C = nb.load(cifti_file)
            data = np.asarray(C.dataobj)
            info = pd.read_csv(tsv_file, sep='\t')
            if data.shape[0] != info.shape[0]:
                raise ValueError(
                    f'{s} {ses_id}: CIFTI rows ({data.shape[0]}) do not '
                    f'match TSV rows ({info.shape[0]})'
                )

            if brain_axis is None:
                brain_axis = C.header.get_axis(1)
            elif C.shape[1] != brain_axis.size:
                raise ValueError(
                    f'{s} {ses_id}: brain axis length ({C.shape[1]}) does '
                    f'not match first run ({brain_axis.size})'
                )

            info = info.copy()
            if 'sess' in info.columns:
                info['orig_sess'] = info['sess']
            if 'run' in info.columns:
                info['orig_run'] = info['run']
            if 'names' in info.columns:
                info['orig_names'] = info['names']

            info['sess'] = ses_out
            info['rest_run'] = run
            info['run'] = run
            if 'half' in info.columns:
                info['half'] = 2 - (run % 2)

            if 'names' in info.columns:
                base_names = info['orig_names'].astype(str).str.replace(
                    r'-run\d+$', '', regex=True)
            else:
                base_names = pd.Series(
                    C.header.get_axis(0).name.astype(str), index=info.index)
            info['names'] = [
                f'{name}-run{run:02d}' for name in base_names
            ]

            data_list.append(data)
            info_list.append(info)
            row_names.extend(info['names'].tolist())
            included_runs.append(run)

        if not data_list:
            print(f'No complete resting FC inputs for {s}, skip ...')
            continue

        data_all = np.vstack(data_list)
        info_all = pd.concat(info_list, ignore_index=True)
        row_axis = nb.cifti2.ScalarAxis(row_names)
        header = nb.Cifti2Header.from_axes((row_axis, brain_axis))
        cifti_out = nb.Cifti2Image(dataobj=data_all, header=header)

        nb.save(
            cifti_out,
            f'{dest_dir}/{s}_space-{space}_{ses_out}_{type}'
            f'_binarized-{binarized}.dscalar.nii'
        )
        info_all.to_csv(
            f'{dest_dir}/{s}_{ses_out}_{type}.tsv',
            sep='\t',
            index=False
        )

        summary.append({
            'participant_id': s,
            'ses_id': ses_out,
            'type': type,
            'n_runs': len(included_runs),
            'runs': ','.join(map(str, included_runs)),
            'n_rows': data_all.shape[0]
        })

    return pd.DataFrame(summary)


def extract_residual_timeseries(source_dir, ses_id='EPROJ', type='Tseries',
                                space='MNISymC3', subj=None,
                                skip_existing=False):
    randy15_dataset = DataSetRANDY15(data_dir)
    myatlas, _ = am.get_atlas(space)

    # Extract the data for each participant
    T = randy15_dataset.get_participants()
    T = T.loc[T["complete_task"] == 1]
    if subj is not None:
        if isinstance(subj, str):
            subj = [subj]
        T = T.loc[T.participant_id.isin(subj)]

    for row in T.itertuples(index=False):
        sub_id = row.participant_num
        s = row.participant_id
        ses_dir = f'{source_dir}/{sub_id}/{ses_id}'
        dest_dir = randy15_dataset.data_dir.format(s)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)
        dest_cifti = (f'{dest_dir}/{s}_space-{space}_ses-{ses_id}'
                      f'_Residuals.dscalar.nii')
        dest_tsv = f'{dest_dir}/{s}_ses-{ses_id}_Residuals.tsv'

        if not os.path.exists(ses_dir):
            print(f"{s} doesn't have session {ses_id} data, skip ...")
            continue

        if skip_existing and os.path.exists(dest_cifti) and os.path.exists(dest_tsv):
            print(f'{s} ses-{ses_id} residuals already exist, skip ...')
            continue

        print(f'Extract {s}')
        file_list_L = sorted([f for f in os.listdir(ses_dir) if f.startswith(f"lh.{ses_id}")],
                        key=lambda x: int(re.search(rf'{ses_id}(\d+)', x).group(1)))
        file_list_R = sorted([f for f in os.listdir(ses_dir) if f.startswith(f"rh.{ses_id}")],
                        key=lambda x: int(re.search(rf'{ses_id}(\d+)', x).group(1)))
        assert len(file_list_L) == len(file_list_R), \
            "file number of left / right doesn't match!"
        
        data, num_tpoint = [],[]
        info = pd.DataFrame()
        for file_L, file_R in zip(file_list_L, file_list_R):
            assert file_L.split('.')[1] == file_R.split('.')[1], \
                    "run id for L / R hemisphere doesn't match!"
            run_id = int(re.search(r'\d+$', file_L.split('.')[1]).group())
            print(f'-- Start processing {s} {ses_id} run {run_id} ...')

            # Option 2 - us wb command resample
            lh_data = nb.load(f'{ses_dir}/{file_L}').get_fdata()
            lh_data = lh_data.reshape(-1, lh_data.shape[-1], order='F')
            rh_data = nb.load(f'{ses_dir}/{file_R}').get_fdata()
            rh_data = rh_data.reshape(-1, rh_data.shape[-1], order='F')

            if np.equal(lh_data, rh_data).all():
                print(f'{sub_id} ses {file_L} and {file_R} data identical!')
                continue
            
            print(f'   Converting data to fsaverage6 func.gii')
            lh_gii = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft', 
                                        column_names=[f'time_{i+1}' for i in range(lh_data.shape[1])])
            rh_gii = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight', 
                                        column_names=[f'time_{i+1}' for i in range(rh_data.shape[1])])
            nb.save(lh_gii, ses_dir + f'/{sub_id}_tmp_{ses_id}{run_id}.L.func.gii')
            nb.save(rh_gii, ses_dir + f'/{sub_id}_tmp_{ses_id}{run_id}.R.func.gii')

            print(f'   Mapping data from fsaverage6 space to fs32k space ...')
            lh_fs32k_gii = ses_dir + f'/{sub_id}_{ses_id}{run_id}_wb.L.func.gii'
            rh_fs32k_gii = ses_dir + f'/{sub_id}_{ses_id}{run_id}_wb.R.func.gii'
            lh_fs32k_data = fsaverage6_to_fs32k(ses_dir + f'/{sub_id}_tmp_{ses_id}{run_id}.L.func.gii', 'L', lh_fs32k_gii,
                                type='metric', return_data_only=True)
            rh_fs32k_data = fsaverage6_to_fs32k(ses_dir + f'/{sub_id}_tmp_{ses_id}{run_id}.R.func.gii', 'R', rh_fs32k_gii,
                                type='metric', return_data_only=True)
            print(f'   Remove tmp files, Done.')
            os.remove(ses_dir + f'/{sub_id}_tmp_{ses_id}{run_id}.L.func.gii')
            os.remove(ses_dir + f'/{sub_id}_tmp_{ses_id}{run_id}.R.func.gii')

            print(f'   Combine Left and Right hemisphere / write info ...')
            this_data = np.hstack([lh_fs32k_data[:,myatlas.vertex_mask[0]], 
                              rh_fs32k_data[:,myatlas.vertex_mask[1]]])
            this_num_tpoint = [f'T{i+1:04}' for i in range(this_data.shape[0])]

            this_info = pd.DataFrame({'sn': [s] * this_data.shape[0],
                                        'run': [run_id] * this_data.shape[0],
                                        'timepoint': this_num_tpoint,
                                        'task': ses_id,
                                        'time_id': [i+1 for i in range(this_data.shape[0])],
                                        'names': this_num_tpoint,
                                        'type':'residuals'})
            info = pd.concat([info, this_info], ignore_index=True)
            data.append(this_data)
            num_tpoint += this_num_tpoint

        C = myatlas.data_to_cifti(np.vstack(data), num_tpoint)
        nb.save(C, dest_cifti)
        info.to_csv(dest_tsv, sep='\t', index=False)
    
        print('-- Done!')


def _fix_block_to_bool(fix_block):
    """Convert a fix_block column to a boolean mask."""
    if pd.api.types.is_bool_dtype(fix_block):
        return fix_block.fillna(False).to_numpy(dtype=bool)

    if pd.api.types.is_numeric_dtype(fix_block):
        return fix_block.fillna(0).to_numpy() != 0

    values = fix_block.fillna('').astype(str).str.strip().str.lower()
    true_values = {'1', 'true', 't', 'yes', 'y', 'fix', 'fixation'}
    false_values = {'0', 'false', 'f', 'no', 'n', 'nonfix', 'non-fix',
                    'non_fix', 'task', 'none', ''}
    valid_values = true_values | false_values
    unknown = sorted(set(values.unique()) - valid_values)
    if unknown:
        raise ValueError(f'Unknown fix_block values: {unknown}')

    return values.isin(true_values).to_numpy()


def _make_residual_split_info(info, source_indices, output_type,
                              ses_id, run_length):
    info_out = info.iloc[source_indices].copy().reset_index(drop=True)
    n_timepoints = len(info_out)
    n_runs = int(n_timepoints / run_length)

    for col in ['run', 'time_id', 'timepoint', 'names', 'task', 'type']:
        if col in info_out.columns and f'orig_{col}' not in info_out.columns:
            info_out[f'orig_{col}'] = info_out[col].to_numpy()

    time_id = np.arange(1, n_timepoints + 1)
    info_out['run'] = np.repeat(np.arange(1, n_runs + 1), run_length)
    info_out['time_id'] = time_id
    info_out['timepoint'] = [f'T{i:04d}' for i in time_id]
    info_out['task'] = ses_id
    info_out['type'] = output_type
    info_out['names'] = [
        f'{output_type}_run-{run:02d}_T{time:04d}'
        for run, time in zip(
            info_out['run'],
            np.tile(np.arange(1, run_length + 1), n_runs)
        )
    ]
    info_out['source_index'] = source_indices + 1

    return info_out


def _format_run_label(run):
    try:
        run_num = float(run)
        if run_num.is_integer():
            return f'{int(run_num):02d}'
        return str(run)
    except (TypeError, ValueError):
        return str(run)


def _make_trimmed_session_info(info, source_indices, output_type, run_col):
    info_out = info.iloc[source_indices].copy().reset_index(drop=True)

    for col in ['time_id', 'timepoint', 'names', 'type']:
        if col in info_out.columns and f'orig_{col}' not in info_out.columns:
            info_out[f'orig_{col}'] = info_out[col].to_numpy()

    time_id = info_out.groupby(run_col).cumcount() + 1
    info_out['time_id'] = time_id
    info_out['timepoint'] = [f'T{i:04d}' for i in time_id]
    info_out['type'] = output_type
    info_out['names'] = [
        f'run-{_format_run_label(run)}_T{time:04d}'
        for run, time in zip(info_out[run_col], time_id)
    ]
    info_out['source_index'] = source_indices + 1

    return info_out


def _zstandardize_split_source_within_run(data, info, source_indices,
                                          run_col='run'):
    """Z-standardize selected rows within each original run."""
    source_info = info.iloc[source_indices]
    standardized = np.empty(
        (len(source_indices), data.shape[1]),
        dtype=np.float32
    )

    for run in source_info[run_col].drop_duplicates():
        source_pos = np.flatnonzero(source_info[run_col].to_numpy() == run)
        run_idx = source_indices[source_pos]
        run_data = data[run_idx, :].astype(np.float32, copy=False)
        run_mean = np.nanmean(run_data, axis=0, keepdims=True)
        run_std = np.nanstd(run_data, axis=0, keepdims=True)
        run_std[run_std == 0] = 1
        standardized[source_pos, :] = (run_data - run_mean) / run_std

    return standardized


def trim_session_by_run_length(ses_id='MOTOR', trimmed_length=410,
                               space='fs32k', type='Residuals',
                               output_type=None, run_col='run',
                               subj=None, overwrite=False):
    """Trim each run in a concatenated session file to a fixed length.

    For each subject, this reads:
        sub-*_space-{space}_ses-{ses_id}_{type}.dscalar.nii
        sub-*_ses-{ses_id}_{type}.tsv

    It keeps the first ``trimmed_length`` time points from every run in
    ``run_col`` and writes a paired CIFTI/TSV. By default the output type is
    ``{type}Trim{trimmed_length}``, for example ``ResidualsTrim410``.
    """
    randy15_dataset = DataSetRANDY15(data_dir)
    ses_id = ses_id.replace('ses-', '')
    output_type = output_type or f'{type}Trim{trimmed_length}'
    if output_type == type and not overwrite:
        raise ValueError(
            'output_type matches input type. Set overwrite=True if you want '
            'to replace the original files.'
        )

    T = randy15_dataset.get_participants()
    if subj is None:
        subjects = T.participant_id
    elif isinstance(subj, str):
        subjects = [subj]
    else:
        subjects = subj

    summary = []
    for s in subjects:
        dest_dir = randy15_dataset.data_dir.format(s)
        cifti_file = (f'{dest_dir}/{s}_space-{space}_ses-{ses_id}_'
                      f'{type}.dscalar.nii')
        tsv_file = f'{dest_dir}/{s}_ses-{ses_id}_{type}.tsv'

        if not os.path.exists(cifti_file) or not os.path.exists(tsv_file):
            print(f'Missing {type} files for {s} ses-{ses_id}, skip ...')
            continue

        print(f'Trimming {type} for {s} ses-{ses_id} to '
              f'{trimmed_length} time points per run')
        C = nb.load(cifti_file)
        data = np.asarray(C.dataobj)
        info = pd.read_csv(tsv_file, sep='\t')

        if run_col not in info.columns:
            raise KeyError(f'{tsv_file} does not contain column {run_col}')
        if data.shape[0] != info.shape[0]:
            raise ValueError(
                f'{s} ses-{ses_id}: CIFTI rows ({data.shape[0]}) do not '
                f'match TSV rows ({info.shape[0]})'
            )

        source_indices = []
        run_lengths = []
        for run in info[run_col].drop_duplicates():
            run_idx = np.flatnonzero(info[run_col].to_numpy() == run)
            run_lengths.append(len(run_idx))
            if len(run_idx) < trimmed_length:
                raise ValueError(
                    f'{s} ses-{ses_id} run {run} only has '
                    f'{len(run_idx)} time points, shorter than '
                    f'trimmed_length={trimmed_length}'
                )
            source_indices.extend(run_idx[:trimmed_length])

        source_indices = np.asarray(source_indices)
        trim_info = _make_trimmed_session_info(
            info, source_indices, output_type, run_col)
        row_axis = nb.cifti2.ScalarAxis(trim_info['names'].tolist())
        header = nb.Cifti2Header.from_axes((row_axis, C.header.get_axis(1)))
        trim_cifti = nb.Cifti2Image(
            dataobj=data[source_indices, :], header=header)

        nb.save(
            trim_cifti,
            f'{dest_dir}/{s}_space-{space}_ses-{ses_id}_'
            f'{output_type}.dscalar.nii'
        )
        trim_info.to_csv(
            f'{dest_dir}/{s}_ses-{ses_id}_{output_type}.tsv',
            sep='\t',
            index=False
        )

        summary.append({
            'participant_id': s,
            'ses_id': ses_id,
            'type': type,
            'output_type': output_type,
            'n_runs': len(run_lengths),
            'original_timepoints': int(np.sum(run_lengths)),
            'kept_timepoints': int(len(source_indices)),
            'trimmed_timepoints': int(np.sum(run_lengths) - len(source_indices)),
            'min_run_length': int(np.min(run_lengths)),
            'max_run_length': int(np.max(run_lengths))
        })

    return pd.DataFrame(summary)


def split_residuals_by_fix_block(ses_id='MOTOR', space='fs32k',
                                 run_length=410, fix_col='fix_block',
                                 input_type='Residuals', subj=None,
                                 run_col='run',
                                 standardize_within_run=True):
    """Split RANDY15 concatenated residuals into matched FIX/non-FIX series.

    The input files are expected in each subject's data directory:
        sub-*_space-{space}_ses-{ses_id}_{input_type}.dscalar.nii
        sub-*_ses-{ses_id}_{input_type}.tsv

    The function uses ``fix_col`` to split the residual rows, z-standardizes
    each split source within each original run, finds the smaller split
    length, keeps only complete ``run_length`` chunks from the beginning of
    both splits, and writes:
        sub-*_space-{space}_ses-{ses_id}_FixResiduals.dscalar.nii
        sub-*_space-{space}_ses-{ses_id}_nonFixResiduals.dscalar.nii
    plus matching TSV files.
    """
    randy15_dataset = DataSetRANDY15(data_dir)
    ses_id = ses_id.replace('ses-', '')
    T = randy15_dataset.get_participants()

    if subj is None:
        subjects = T.participant_id
    elif isinstance(subj, str):
        subjects = [subj]
    else:
        subjects = subj

    summary = []
    for s in subjects:
        dest_dir = randy15_dataset.data_dir.format(s)
        cifti_file = (f'{dest_dir}/{s}_space-{space}_ses-{ses_id}'
                      f'_{input_type}.dscalar.nii')
        tsv_file = f'{dest_dir}/{s}_ses-{ses_id}_{input_type}.tsv'

        if not os.path.exists(cifti_file) or not os.path.exists(tsv_file):
            print(f'Missing {input_type} files for {s} ses-{ses_id}, skip ...')
            continue

        print(f'Splitting {input_type} for {s} ses-{ses_id}')
        C = nb.load(cifti_file)
        data = np.asarray(C.dataobj)
        info = pd.read_csv(tsv_file, sep='\t')

        if fix_col not in info.columns:
            raise KeyError(f'{tsv_file} does not contain column {fix_col}')
        if standardize_within_run and run_col not in info.columns:
            raise KeyError(f'{tsv_file} does not contain column {run_col}')
        if data.shape[0] != info.shape[0]:
            raise ValueError(
                f'{s} ses-{ses_id}: CIFTI rows ({data.shape[0]}) do not '
                f'match TSV rows ({info.shape[0]})'
            )

        fix_mask = _fix_block_to_bool(info[fix_col])
        fix_idx = np.flatnonzero(fix_mask)
        nonfix_idx = np.flatnonzero(~fix_mask)
        matched_runs = min(len(fix_idx), len(nonfix_idx)) // run_length
        keep_timepoints = matched_runs * run_length

        print(f'   FIX={len(fix_idx)}, non-FIX={len(nonfix_idx)}, '
              f'keep={keep_timepoints} ({matched_runs} runs)')
        if keep_timepoints == 0:
            print(f'   No complete {run_length}-timepoint matched run, skip ...')
            continue

        split_specs = [
            ('FixResiduals', fix_idx),
            ('nonFixResiduals', nonfix_idx)
        ]
        brain_axis = C.header.get_axis(1)
        for output_type, all_source_idx in split_specs:
            if standardize_within_run:
                print(f'   Z-standardizing {output_type} within each {run_col}')
                source_data = _zstandardize_split_source_within_run(
                    data, info, all_source_idx, run_col=run_col)
            else:
                source_data = data[all_source_idx, :]

            source_idx = all_source_idx[:keep_timepoints]
            split_data = source_data[:keep_timepoints, :]
            split_info = _make_residual_split_info(
                info, source_idx, output_type, ses_id, run_length)
            split_info['standardized_within_run'] = standardize_within_run
            split_info['standardized_run_col'] = run_col
            row_axis = nb.cifti2.ScalarAxis(split_info['names'].tolist())
            header = nb.Cifti2Header.from_axes((row_axis, brain_axis))
            split_cifti = nb.Cifti2Image(dataobj=split_data, header=header)

            nb.save(
                split_cifti,
                f'{dest_dir}/{s}_space-{space}_ses-{ses_id}_'
                f'{output_type}.dscalar.nii'
            )
            split_info.to_csv(
                f'{dest_dir}/{s}_ses-{ses_id}_{output_type}.tsv',
                sep='\t',
                index=False
            )

        summary.append({
            'participant_id': s,
            'ses_id': ses_id,
            'input_type': input_type,
            'standardized_within_run': standardize_within_run,
            'fix_timepoints': len(fix_idx),
            'nonfix_timepoints': len(nonfix_idx),
            'matched_runs': matched_runs,
            'kept_timepoints': keep_timepoints,
            'trimmed_fix_timepoints': len(fix_idx) - keep_timepoints,
            'trimmed_nonfix_timepoints': len(nonfix_idx) - keep_timepoints
        })

    return pd.DataFrame(summary)


def extract_betas(ses_id='ses-task', type='Tseries', atlas='MNISymC3'):
    randy15_dataset = DataSetRANDY15(data_dir)
    ses = ses_id.split('-')[1].upper()

    # Extract the data for each participant
    T = randy15_dataset.get_participants()
    for row in T.itertuples(index=False):
        sub_id = row.participant_num
        s = row.participant_id
        ses_dir = f'{randy15_dir}/{sub_id}/{ses}'

        print(f'Extract {s}')
        file_list_L = [f for f in os.listdir(ses_dir) if f.startswith("lh.")]
        file_list_R = [f for f in os.listdir(ses_dir) if f.startswith("rh.")]
        assert len(file_list_L) == len(file_list_R), \
            "file number of left / right doesn't match!"
        
        for file_L, file_R in zip(file_list_L, file_list_R):
            assert file_L.split('.')[1] == file_R.split('.')[1], \
                    "run id for L / R hemisphere doesn't match!"
            run_id = int(re.search(r'\d+$', file_L.split('.')[1]).group())
            print(f'-- Start processing {s} {ses_id} run {run_id} ...')

            # Option 2 - us wb command resample
            lh_data = nb.load(f'{ses_dir}/{file_L}').get_fdata()
            lh_data = lh_data.reshape(-1, lh_data.shape[-1], order='F')
            rh_data = nb.load(f'{ses_dir}/{file_R}').get_fdata()
            rh_data = rh_data.reshape(-1, rh_data.shape[-1], order='F')

            if np.equal(lh_data, rh_data).all():
                print(f'{sub_id} ses {file_L} and {file_R} data identical!')
                break
            
            print(f'   Converting data to fsaverage6 func.gii')
            lh_gii = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft', 
                                        column_names=[f'time_{i+1}' for i in range(lh_data.shape[1])])
            rh_gii = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight', 
                                        column_names=[f'time_{i+1}' for i in range(rh_data.shape[1])])
            nb.save(lh_gii, ses_dir + f'/{sub_id}_tmp_{file_L}.L.func.gii')
            nb.save(rh_gii, ses_dir + f'/{sub_id}_tmp_{file_R}.R.func.gii')

            print(f'   Mapping data from fsaverage6 space to fs32k space ...')
            lh_fs32k_gii = ses_dir + f'/P1_REST{run_id}_wb.L.func.gii'
            rh_fs32k_gii = ses_dir + f'/P1_REST{run_id}_wb.R.func.gii'
            fsaverage6_to_fs32k(ses_dir + f'/{sub_id}_tmp_{file_L}.L.func.gii', 'L', lh_fs32k_gii)
            fsaverage6_to_fs32k(ses_dir + f'/{sub_id}_tmp_{file_R}.R.func.gii', 'R', rh_fs32k_gii)
            print(f'   Remove tmp files, Done.')
            os.remove(ses_dir + f'/{sub_id}_tmp_{file_L}.L.func.gii')
            os.remove(ses_dir + f'/{sub_id}_tmp_{file_R}.R.func.gii')

            print(f'   Combine Left and Right hemisphere into a single cifti ...')
            dest_dir = randy15_dataset.estimates_dir.format(s)
            Path(dest_dir).mkdir(parents=True, exist_ok=True)

            try:
                cmd = (f'wb_command -cifti-create-dense-timeseries '
                            f'{dest_dir}/{s}_space-fs32k_{ses_id}_run{run_id}.dtseries.nii '
                            f'-left-metric {lh_fs32k_gii} -right-metric {rh_fs32k_gii} '
                            f'-timestep 1.0 -timestart 0')
                subprocess.run(cmd, shell=True, check=True)
                os.remove(lh_fs32k_gii)
                os.remove(rh_fs32k_gii)

            except:
                print(f"Failed to combine in single CIFTI for {s} {ses_id} run {run_id}!")
        
            print('-- Done!')

def import_contrasts(source_dir, session_info=['EPROJ','LANG','MOTOR','NBACK','TOM','VISME','VODDK'],
                     space='fs32k'):
    """This is to import all run-wise task contrasts
    """
    myatlas, _ = am.get_atlas(space)
    dataset = DataSetRANDY15(data_dir)
    T = dataset.get_participants()

    for i, s in enumerate(T.participant_id):
        orig_subid = T.participant_num[i]
        # Gather all session directories
        betas_dir = Path(source_dir) / 'Task_GLM'

        dest_folder = Path(dataset.contrast_dir.format(s))
        for session_index, session_name in enumerate(session_info):
            ses_dir = betas_dir / orig_subid / session_name
            beta_files = sorted(list(ses_dir.glob(f"{orig_subid}_*_concat[Zz]stats.dtseries.nii")))

            for reg_num, file in enumerate(beta_files):
                # 1. load data and make func.gii for both L/R hemisphere
                match = re.search(f"sm2_(.*?)_concat[Zz]stats", file.stem)
                cond_name = match.group(1) if match else file.stem

                print(f'   Converting data to fsaverage6 func.gii')
                data = nb.load(file).get_fdata()
                num_run = data.shape[0]
                lh_data = np.split(data, 2, axis=1)[0].T
                rh_data = np.split(data, 2, axis=1)[1].T
                
                lh_gii = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft', 
                                            column_names=[f'run_{i+1}' for i in range(lh_data.shape[1])])
                rh_gii = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight', 
                                            column_names=[f'run_{i+1}' for i in range(rh_data.shape[1])])
                nb.save(lh_gii, ses_dir / f'{orig_subid}_tmp_{cond_name}.L.func.gii')
                nb.save(rh_gii, ses_dir / f'{orig_subid}_tmp_{cond_name}.R.func.gii')

                # 2. Convert fsaverage6 func.gii to fslr32k space
                print(f'   Mapping data from fsaverage6 space to fs32k space ...')
                lh_data_fs32k = fsaverage6_to_fs32k(ses_dir / f'{orig_subid}_tmp_{cond_name}.L.func.gii', 'L', 
                                    ses_dir / f'{orig_subid}_{session_name}_{cond_name}_wb.L.func.gii',
                                        type='metric', return_data_only=True)
                rh_data_fs32k = fsaverage6_to_fs32k(ses_dir / f'{orig_subid}_tmp_{cond_name}.R.func.gii', 'R', 
                                    ses_dir / f'{orig_subid}_{session_name}_{cond_name}_wb.R.func.gii',
                                        type='metric', return_data_only=True)
                print(f'   Remove tmp files, Done.')
                os.remove(ses_dir / f'{orig_subid}_tmp_{cond_name}.L.func.gii')
                os.remove(ses_dir / f'{orig_subid}_tmp_{cond_name}.R.func.gii')

                # 3. Combine L/R fslr32k data
                this_data = np.hstack([lh_data_fs32k[:,myatlas.vertex_mask[0]], 
                                        rh_data_fs32k[:,myatlas.vertex_mask[1]]])
                con_nams = [f'{cond_name}-run{r+1:02d}' for r in range(this_data.shape[0])]
                C = myatlas.data_to_cifti(this_data, con_nams)

                # make reg info file
                info = pd.DataFrame({"sn": [s] * len(con_nams),
                                     'run': np.arange(len(con_nams)) + 1,
                                     "cond_id": [reg_num + 1] * len(con_nams) ,
                                     "task_name": [session_name] * len(con_nams),
                                     "task_uni_num": [session_index + 1] * len(con_nams),
                                     "contrast_name": [cond_name] * len(con_nams),
                                     "domain": [session_name] * len(con_nams),
                                     "domain_abbr": [session_name] * len(con_nams)})

                # Write-in beta files in destination folder
                Path(dest_folder).mkdir(parents=True, exist_ok=True)
                info.to_csv(dest_folder / f"{s}_{session_name}_{cond_name}_Contrasts.tsv", sep="\t", index=False)
                dest_file = dest_folder / f"{s}_{session_name}_{cond_name}_fs32k_sm2_Zmap.dscalar.nii"
                nb.save(C, dest_file)
                print(f'Successfully import run-wise contrasts for subject {s}: {dest_file}!')


def import_betas(source_dir, space='fs32k'):
    myatlas, _ = am.get_atlas(space)
    dataset = DataSetRANDY15(data_dir)
    reg_info = pd.read_csv(dataset.base_dir + '/regressor_info.tsv', sep='\t')

    T = dataset.get_participants()
    for i, s in enumerate(T.participant_id):
        orig_subid = T.participant_num[i]
        # Gather all session directories
        betas_dir = Path(source_dir) / 'Task_GLM_PE'
        session_info = ['EPROJ','LANG','MOTOR','NBACK','TOM','VISME','VODDK']

        dest_folder = Path(dataset.estimates_dir.format(s) + "/ses-task")
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder, exist_ok=True)

        reginfo_data = []  # To store reginfo entries
        # Each task domain
        for session_index, session_name in enumerate(session_info):
            ses_dir = betas_dir / orig_subid / session_name
            if not ses_dir.exists():
                print(f'Subject {s} domain {session_name} does not exist, skip ...')
                continue

            run_folders = sorted([f for f in ses_dir.iterdir() if f.is_dir()], 
                                    key=lambda f: int(f.name))
            session_run_start = int(np.unique(reg_info.loc[reg_info['task_domain'] 
                                    == session_name]['start_run_num'])[0])

            # Each run
            for run_num, run_dir in enumerate(run_folders):
                lh_pe = list((run_dir/'lh').glob(f"pe*.nii.gz"))
                rh_pe = list((run_dir/'rh').glob(f"pe*.nii.gz"))
                # Filter and sort by the number in the filename
                odd_lh_pe = sorted([f for f in lh_pe if int(f.name.removeprefix('pe').removesuffix('.nii.gz')) % 2 == 1],
                                    key=lambda f: int(f.name.removeprefix('pe').removesuffix('.nii.gz')))
                odd_rh_pe = sorted([f for f in rh_pe if int(f.name.removeprefix('pe').removesuffix('.nii.gz')) % 2 == 1],
                                    key=lambda f: int(f.name.removeprefix('pe').removesuffix('.nii.gz')))
                assert len(odd_lh_pe) == len(odd_rh_pe), "PE number doesn't match for L/R hemisphere!"
                global_run_counter = session_run_start + run_num

                # Import beta weights
                for reg_num, (lh_file, rh_file) in enumerate(zip(odd_lh_pe, odd_rh_pe)):
                    this_reg_info = reg_info.loc[(reg_info['task_domain'] == session_name) & 
                                             (reg_info['reg_num'] == reg_num + 1)]
                    
                    # 1. load data and make func.gii for both L/R in fsaverage space
                    print(f'   Importing subject {s}, domain {session_name}, run {run_num+1}, PE {reg_num+1}')
                    lh_data = nb.load(lh_file).get_fdata().reshape(-1, 1, order='F')
                    rh_data = nb.load(rh_file).get_fdata().reshape(-1, 1, order='F')
                    
                    if np.equal(lh_data, rh_data).all():
                        print(f'{s} {session_name} run {run_num+1} PE {reg_num+1} data are identical!')
                        break

                    lh_gii = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft')
                    rh_gii = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight')
                    tmp_fname = f'tmp_{orig_subid}_{session_name}{run_num+1}_PE{reg_num+1}'
                    nb.save(lh_gii, run_dir / f'{tmp_fname}.L.func.gii')
                    nb.save(rh_gii, run_dir / f'{tmp_fname}.R.func.gii')

                    # 2. Convert fsaverage func.gii to fslr32k space
                    lh_data_fs32k = fsaverage6_to_fs32k(run_dir / f'{tmp_fname}.L.func.gii', 'L', 
                                        run_dir / f'{tmp_fname}_wb.L.func.gii', type='metric', return_data_only=True)
                    rh_data_fs32k = fsaverage6_to_fs32k(run_dir / f'{tmp_fname}.R.func.gii', 'R', 
                                        run_dir / f'{tmp_fname}_wb.R.func.gii', type='metric', return_data_only=True)
                    print(f'   Remove tmp files, Done.')
                    os.remove(run_dir / f'{tmp_fname}.L.func.gii')
                    os.remove(run_dir / f'{tmp_fname}.R.func.gii')

                    # 3. Combine L/R fslr32k data
                    this_data = np.hstack([lh_data_fs32k[:,myatlas.vertex_mask[0]], 
                                           rh_data_fs32k[:,myatlas.vertex_mask[1]]])
                    cond_name = this_reg_info['cond_name'].iloc[0]
                    C = myatlas.data_to_cifti(this_data, [cond_name])

                    # make reg info file
                    reg_id = int(this_reg_info['reg_id'].iloc[0])
                    reginfo_data.append({"sn": run_num+1,
                                        "run": global_run_counter,
                                        "task_name": session_name,
                                        "cond_name": cond_name,
                                        "task_num": session_index+1,
                                        "cond_num": reg_id,
                                        "reg_id": reg_id,
                                        "reg_num": reg_num + 1,
                                        "half": 2 if (run_num+1) % 2 == 0 else 1})

                    # Write-in beta files in destination folder
                    dest_file = dest_folder / f"{s}_ses-task_run-{global_run_counter:02d}_reg-{reg_id:02d}_beta.dscalar.nii"
                    Path(dest_folder).mkdir(parents=True, exist_ok=True)
                    nb.save(C, dest_file)

            print(f'Successfully import betas for subject {s}, task domain {session_name}, total run: {len(run_folders)}')

        # Write-in reg info file
        reginfo_df = pd.DataFrame(reginfo_data, columns=["sn", "run", "task_name", "cond_name", 
                                                         "task_num", "reg_id","reg_num", "half"])
        reginfo_df.to_csv(dest_folder / f"{s}_ses-task_reginfo.tsv", sep="\t", index=False)
        print(f"Saved task reginfo file for subject {s}")


def integrate_betas(space='fs32k'):
    myatlas, _ = am.get_atlas(space)
    dataset = DataSetRANDY15(data_dir)

    T = dataset.get_participants()
    for i, s in enumerate(T.participant_id):
        reg_info = pd.read_csv(dataset.estimates_dir.format(s) + 
                               f'/ses-task/{s}_ses-task_reginfo.tsv', sep='\t')

        dest_folder = Path(dataset.estimates_dir.format(s) + "/ses-task")
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder, exist_ok=True)

        data, conds = [],[]
        for index, row in reg_info.iterrows():
            file = f'/{s}_ses-task_run-{row.run:02d}_reg-{row.reg_id:02d}_beta.dscalar.nii'
            pe_data = nb.load(dataset.estimates_dir.format(s) + 
                              f'/ses-task/{file}').get_fdata()
            data.append(pe_data)
            conds.append(f'{row.task_name}{row.run}_{row.cond_name}')
        
        data = np.vstack(data)
        # Write-in beta files in destination folder
        C = myatlas.data_to_cifti(data, conds)
        nb.save(C, dest_folder / f'{s}_ses-task_beta.dscalar.nii')
        print(f'Successfully integrated subject {s} betas into a single cifti.')


def import_resms(source_dir, space='fs32k'):
    myatlas, _ = am.get_atlas(space)
    dataset = DataSetRANDY15(data_dir)
    reg_info = pd.read_csv(dataset.base_dir + '/regressor_info.tsv', sep='\t')
    T = dataset.get_participants()

    for i, participant in enumerate(T.participant_id):
        # Setup destination file
        orig_subid = T.participant_num[i]
        res_dir = Path(f'{source_dir}/Task_GLM_PE/{orig_subid}')
        dest_folder = dataset.estimates_dir.format(participant) + "/ses-task"
        dest_file = f'{dest_folder}/{participant}_ses-task_resms.dscalar.nii'
        if not Path(dest_folder).exists():
            os.makedirs(dest_folder, exist_ok=True)

        lh_resms_files = [f for f in res_dir.rglob('sigmasquareds.nii.gz')
                        if f.parent.name == 'lh']
        rh_resms_files = [f for f in res_dir.rglob('sigmasquareds.nii.gz')
                        if f.parent.name == 'rh']

        lh_resms_data, rh_resms_data = [], []
        for (lh_file, rh_file) in zip(lh_resms_files, rh_resms_files):
            lh_data = nb.load(lh_file).get_fdata().reshape(-1, 1, order='F')
            rh_data = nb.load(rh_file).get_fdata().reshape(-1, 1, order='F')
            lh_resms_data.append(lh_data)
            rh_resms_data.append(rh_data)

        lh_resms_data = np.mean(lh_resms_data, axis=0)
        rh_resms_data = np.mean(rh_resms_data, axis=0)

        if np.equal(lh_data, rh_data).all():
            break

        print(f'   Making gifti in fsaverage space')
        lh_gii = nt.make_func_gifti(lh_resms_data, anatomical_struct='CortexLeft')
        rh_gii = nt.make_func_gifti(rh_resms_data, anatomical_struct='CortexRight')
        nb.save(lh_gii, res_dir / 'tmp_resms.L.func.gii')
        nb.save(rh_gii, res_dir / 'tmp_resms.R.func.gii')

        # 2. Convert fsaverage func.gii to fslr32k space
        lh_data_fs32k = fsaverage6_to_fs32k(res_dir / 'tmp_resms.L.func.gii', 'L', 
                            res_dir / 'tmp_resms_wb.L.func.gii', type='metric', return_data_only=True)
        rh_data_fs32k = fsaverage6_to_fs32k(res_dir / 'tmp_resms.R.func.gii', 'R', 
                            res_dir / 'tmp_resms_wb.R.func.gii', type='metric', return_data_only=True)
        print(f'   Remove tmp files, Done.')
        os.remove(res_dir / 'tmp_resms.L.func.gii')
        os.remove(res_dir / 'tmp_resms.R.func.gii')

        # 3. Combine L/R fslr32k data
        this_data = np.hstack([lh_data_fs32k[:,myatlas.vertex_mask[0]], 
                                rh_data_fs32k[:,myatlas.vertex_mask[1]]])
        C = myatlas.data_to_cifti(this_data, ['resms'])
        
        nb.save(C, dest_file)
        print(f'Copied resms file to {dest_file} for participant {participant}')


def extract_task_contrasts(src_dir, space='fs32k'):
    task_domains = ["Episodic_Projection","Motor","N-Back","Oddball",
                    "Sentence_Processing","Theory_of_Mind","Visual"]
    atlas, _ = am.get_atlas(space)
    atlas.calculate_symmetry()
    randy_dataset = DataSetRANDY15(data_dir)
    T = randy_dataset.get_participants()

    for i, sub in enumerate([f'P{i+1:02}' for i in range(15)]):
        subj_num = T.iloc[i].participant_id
        contrasts_dir = src_dir.format(sub)
        all_files = [f for f in Path(contrasts_dir).iterdir() if f.is_file()]
        
        dest_dir = randy_dataset.task_contrasts_dir.format(subj_num)
        Path(dest_dir).mkdir(parents=True, exist_ok=True)

        contrasts, con_nams, domains = [], [], []
        for file in all_files:
            data = convert_fs6_indivPar_to_fs32k(str(file))
            match = re.search(f"{sub}_(.*?)_fsaverage6", file.stem)
            file_name = match.group(1) if match else file.stem
            this_task_domain = next((t for t in task_domains 
                                     if file_name.startswith(t + "_")), None)
            this_contrast = file_name[len(this_task_domain) + 1:] if this_task_domain else None

            contrasts.append(data)
            con_nams.append(this_contrast)
            domains.append(this_task_domain)

        info = pd.DataFrame({'sn': [subj_num] * len(con_nams),
                             'cond_id': np.arange(len(con_nams))+1,
                             'contrast_name': con_nams,
                             'domain': domains})
        
        C = atlas.data_to_cifti(np.stack(contrasts), con_nams)
        nb.save(C, dest_dir + f'/{subj_num}_AllContrasts_fs32k_sm2_Zmap.dscalar.nii')
        info.to_csv(dest_dir + f'/{subj_num}_AllContrasts.tsv', sep='\t', index=False)
        print(f'Done subject {subj_num} all contrast')


def smooth_randytask_fs32k(ses_id='ses-task', type='CondRun', smooth=None, kernel='fwhm'):
    dataset = DataSetRANDY15(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Smooth data for {s} fs32k {ses_id} in {smooth} {kernel} ...')

        start = time.perf_counter()
        file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'
        ut.smooth_fs32k_data(file, smooth=smooth, kernel=kernel)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")


def mask_randytask_fs32k(ses_id='ses-s1', type='CondHalf', high_percent=0.1, low_percent=0.1,
                    smooth=None, z_transfer=False, binarized=False):
    dataset = DataSetRANDY15(data_dir)
    T = dataset.get_participants()

    for s in T.participant_id:
        print(f'Mask data for {s} fs32k {ses_id} in high {high_percent} low {low_percent} ...')

        start = time.perf_counter()
        if smooth is not None:
            file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}_desc-sm{smooth}.dscalar.nii'
        else:
            file = dataset.data_dir.format(s) + f'/{s}_space-fs32k_{ses_id}_{type}.dscalar.nii'

        ut.mask_fs32k_data(file, high_percent=high_percent, low_percent=low_percent,
                           z_transfer=z_transfer, binarized=binarized)
        finish = time.perf_counter()
        elapse = time.strftime('%H:%M:%S', time.gmtime(finish - start))
        print(f"- Done subject {s} - time {elapse}.")

if __name__ == "__main__":
    dataset = DataSetRANDY15(data_dir)

    # Extract the data for each participant
    T = dataset.get_participants()
    T = T.loc[T["complete_task"] == 1]
    # Make func.gii
    # for ses in range(1,18):
    #     # Option 1 - use CBIG conversion
    #     lh_data = spio.loadmat(randy15_dir + f'/P1_REST{ses}.mat')['left_ts']
    #     rh_data = spio.loadmat(randy15_dir + f'/P1_REST{ses}.mat')['right_ts']
    #     gifti_L = nt.make_func_gifti(lh_data, anatomical_struct='CortexLeft', 
    #                                  column_names=[f'time_{i+1}' for i in range(lh_data.shape[1])])
    #     gifti_R = nt.make_func_gifti(rh_data, anatomical_struct='CortexRight',
    #                                  column_names=[f'time_{i+1}' for i in range(rh_data.shape[1])])
    #     nb.save(gifti_L, randy15_dir + f'/P1_REST{ses}_yeo.L.func.gii')
    #     nb.save(gifti_R, randy15_dir + f'/P1_REST{ses}_yeo.R.func.gii')
    # extract_timeseries(ses_id='ses-rest', type='Tseries', atlas='fs32k')
    
    # for i in range(24):
    #     import_rs_timeseries(run_id=i+1, type='Tseries', atlas='fs32k')
    
    ## Extract Randy15 task residuals and functional connectivity
    session_info = ['NBACK','TOM','VISME','VODDK']
    for ses in session_info:
        extract_residual_timeseries(randy15_dir + '/Task_Residuals',
                                    ses_id=ses, type='Residuals', space='fs32k')
        # trim_session_by_run_length(ses_id=ses, trimmed_length=410,
        #                            type='Residuals', space='fs32k',
        #                            subj=T.participant_id)
        # split_residuals_by_fix_block(ses_id=ses, space='fs32k',
        #                              run_length=410, subj=T.participant_id)
        # for conn_type in ['Ico642ResTrim410Run', 'Ico642ResRun']:
        #     conn.get_connectivity_fingerprint('RANDY15', type=conn_type,
        #                                       space='fs32k',
        #                                       ses_id=f'ses-{ses}',
        #                                       subj=T.participant_id)

    ## Extract Randy15 task contrasts
    # import_contrasts('/data/tge/dzhi/projects/RANDY15')
    
    ## Import Randy15 task data (beta maps)
    # import_betas('/data/tge/dzhi/projects/RANDY15')
    # import_resms('/data/tge/dzhi/projects/RANDY15')
    # integrate_betas(space='fs32k')

    ## Extract Randy15 task data
    # dataset.extract_all(ses_id='ses-task', type='CondRun', atlas='fs32k')
    # dataset.extract_all(ses_id='ses-task', type='CondAll', atlas='fs32k')

    ## Test load task betas
    # dat, info, ds = ds.get_dataset(base_dir, 'RANDY15', atlas='fs32k', sess='ses-task',
    #                                 type='CondRun', subj=['sub-10','sub-11','sub-12'], smooth=None)
    
    ## Group average data
    # dataset.group_average_data(ses_id='ses-task', type='CondAll', atlas='fs32k')

    ## Concatenate resting runs
    # concatenate_rs_timeseries([1,2,3,4,14,15,13,12,11], start_point=0, duration=3387)
    concatenate_rest_fc(rest_runs=None, space='fs32k', type='Ico642Run',
                        ses_out='ses-rest', binarized='0.1',
                        subj=T.participant_id)

    # Extract Randy15 rest data rsFC
    # for i in ["3387s2"]:
    #     conn.get_connectivity_fingerprint('RANDY15', type='Ico642Run', space='fs32k', ses_id=f'ses-rest{i}',
    #                                 subj=None)

    # -- Get connectivity fingerprint --
    dname = 'HCP'
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest1')
    # conn.get_connectivity_fingerprint(dname,
    #                                   type='Net67Run', space='MNISymC2', ses_id='ses-rest2')

    # for t in [162]:
    #     get_hcp_fs32k_rsfc(type=f'Ico{t}Run', space='fs32k', ses_id='ses-rest1', 
    #                     subj_list='/subj_list/HCP40_training_set.tsv',
    #                     smooth=None, kernel='fwhm', thres=None, keeptop=False)
    #     get_hcp_fs32k_rsfc(type=f'Ico{t}Run', space='fs32k', ses_id='ses-rest2',
    #                     subj_list='/subj_list/test_split/HCP923_test_set_split_1.tsv',
    #                     smooth=4, kernel='fwhm', thres=0.1, keeptop=False)


    #  -- fs32k smoothing (cortex)
    # smooth_hcp_fs32k(hcp_dir + '/subj_list/HCP40_validation_set.tsv', ses_id='ses-rest1',
    #                 type=f'Tseries', smooth=4, kernel='fwhm', return_data_only=False)
    # smooth_hcp_fs32k(hcp_dir + '/subj_list/HCP40_validation_set.tsv', ses_id='ses-rest2',
    #                 type=f'Tseries', smooth=4, kernel='fwhm', return_data_only=False)
    
    #  -- fs32k binarizing (cortex)
    for i in ["3387s2"]:
        for t in [0.1]:
            binarize_rsfc_fs32k(data_dir + '/participants.tsv', ses_id=f'ses-rest{i}',
                                type='Ico642Run', smooth=None, kernel='fwhm', thres=t)
    # name = 'HCP'

    # smooth_randytask_fs32k(ses_id='ses-task', type='CondRun', smooth=4, kernel='fwhm')
    mask_randytask_fs32k(ses_id='ses-task', type='CondRun', high_percent=0.1, low_percent=0.1,
                smooth='4fwhm', z_transfer=True, binarized=False)

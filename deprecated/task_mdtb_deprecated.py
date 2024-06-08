
def make_tinfo_file():
    """
        Creates info file for time series data in the following format:

        run	timepoint	task	time_id
        1	T0001	rest	1
        1	T0002	rest	2
        1	T0003	rest	3
        1	T0004	rest	4
        1	T0005	rest	5
        1	T0006	rest	6
        1	T0007	rest	7
    """

    T = pd.read_csv(f'{fusion_dir}/participants.tsv', delimiter='\t')
    for subject in T.iterrows():
        subject = subject[1].participant_id
        for session in sessions:
            tinfo_file = f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_tinfo.tsv"
            # Load last run
            run = runs[-1]
            img = nib.load(f"{fusion_dir}/derivatives/{subject}/estimates/ses-s{session}/{subject}_ses-s{session}_run-{run}.nii")
            # Get the number of timepoints
            timepoints = img.shape[-1]

            with open(tinfo_file, 'w') as f:
                f.write('run\ttimepoint\ttask\ttime_id\n')
                for run in runs:
                    for timepoint in range(1, timepoints+1):
                        # make time_id continuous across runs
                        time_id = (int(run) - 1) * timepoints + timepoint
                        f.write(f"{int(run)}\tT{time_id:04d}\ttask\t{time_id}\n")
                
            print(f'Created {tinfo_file}')
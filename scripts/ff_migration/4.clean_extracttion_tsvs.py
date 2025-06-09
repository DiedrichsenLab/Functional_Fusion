import pandas as pd
from pathlib import Path
import importlib
import re
# Base directory resolution
base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = '/cifs/diedrichsen/data/FunctionalFusion_new'
if not Path(base_dir).exists():
    base_dir = 'Y:/data/FunctionalFusion_new'

# Dataset names
dataset_names = [
    "MDTB",
    "Pontine",
    "Nishimoto",
    "IBC",
    "WMFS",
    "Demand",
    "Somatotopic",
    "DMCC",
    "Language",
    "HCPur100",
    "Social"
]
dataset_names = [
    "Pontine"
]

valid_types = ["TaskRun", "TaskHalf", "TaskAll", "CondRun", "CondHalf", "CondAll"]


# Load dataset description
desc_path = Path(base_dir) / 'dataset_description.tsv'
df = pd.read_csv(desc_path, sep='\t')



# Import general module
dataset_module = importlib.import_module("Functional_Fusion.dataset")


for name in dataset_names:
    print(f"Processing dataset: {name}")
    row = df[df['name'] == name]
    class_name = row['class_name'].values[0]
    dataset_dir = Path(base_dir) / name

    # initialize dataset class
    DatasetClass = getattr(dataset_module, class_name)
    dataset_object = DatasetClass(str(dataset_dir))

    extract_dir = dataset_dir / 'derivatives' / 'ffextract'
    for subdir in extract_dir.iterdir():
        print(f"Processing subject: {subdir.name}")
        if not subdir.is_dir() or not subdir.name.startswith('sub-') or subdir.name.startswith('group'):
            print(f"Skipping {subdir.name} as it is not a valid subject directory.")
            continue
        sub = subdir.name

        # Loop through exisiting extarcted .tsv files
        for file in subdir.glob('*.tsv'):
            fname = file.name

            if any(tag in fname for tag in valid_types):
                # Extract session and type from filename
                ses_match = re.search(r'ses-[a-zA-Z0-9]+', fname)
                type_match = next((tag for tag in valid_types if tag in fname), None)

                if not ses_match or not type_match:
                    print(f"[SKIP] Invalid file name format: {fname}")
                    continue
                ses = ses_match.group()

                print(f"Extracting: {name} | {sub} | {ses} | {type_match}")

                participants = dataset_object.get_participants()
                sub_index = participants.index[participants['participant_id'] == sub].tolist()

                if not sub_index:
                    print(f"WARNING Subject {sub} not found in participants.tsv")
                    continue

                try:
                    dataset_object.extract_all(
                        ses_id=ses,
                        type=type_match,
                        atlas="MNISymC3",
                        smooth=None,
                        subj=sub_index
                    )
                except Exception as e:
                    print(f"[ERROR] extract_all failed for {fname}: {e}")


#!/usr/bin/env python3
"""Run and profile HCP fs32k rsFC extraction.

The profiler records wall time and resident memory (RSS) for the complete
function and for each major processing step.  By default, RSS includes child
processes, which is important when ``wb_command`` is used for smoothing.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import threading
import time
from contextlib import contextmanager
from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
import psutil

import Functional_Fusion.atlas_map as am
import Functional_Fusion.connectivity as conn
import Functional_Fusion.dataset as ds
from Functional_Fusion.dataset import DataSetHcpResting


DEFAULT_HCP_DIR = "/data/tge/Tian/HCP_img"
DEFAULT_SUBJECT_LIST = "/subj_list/HCP40_training_KONG2019.tsv"


class StepProfiler:
    """Continuously sample RSS and summarize nested timed code sections."""

    def __init__(self, interval: float = 0.1, include_children: bool = True):
        if interval <= 0:
            raise ValueError("The sampling interval must be greater than zero.")
        self.interval = interval
        self.include_children = include_children
        self.process = psutil.Process(os.getpid())
        self.samples: list[tuple[float, int]] = []
        self.records: list[dict[str, object]] = []
        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None

    def _rss_bytes(self) -> int:
        try:
            total = self.process.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return 0

        if self.include_children:
            for child in self.process.children(recursive=True):
                try:
                    total += child.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
        return total

    def _take_sample(self) -> tuple[float, int]:
        sample = (time.perf_counter(), self._rss_bytes())
        with self._lock:
            self.samples.append(sample)
        return sample

    def _sample_loop(self) -> None:
        while not self._stop_event.wait(self.interval):
            self._take_sample()

    def start(self) -> None:
        if self._thread is not None:
            raise RuntimeError("Profiler has already been started.")
        self._take_sample()
        self._thread = threading.Thread(target=self._sample_loop, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._take_sample()
        self._stop_event.set()
        self._thread.join()
        self._thread = None

    @contextmanager
    def measure(self, step: str, subject: str = "", scope: str = "step"):
        """Measure one section while using the shared memory-sampling thread."""
        started_at, start_rss = self._take_sample()
        status = "ok"
        error = ""
        try:
            yield
        except BaseException as exc:
            status = "error"
            error = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            finished_at, end_rss = self._take_sample()
            with self._lock:
                values = [
                    rss
                    for timestamp, rss in self.samples
                    if started_at <= timestamp <= finished_at
                ]
            peak_rss = max(values, default=max(start_rss, end_rss))
            mib = 1024**2
            self.records.append(
                {
                    "sequence": len(self.records) + 1,
                    "scope": scope,
                    "subject": subject,
                    "step": step,
                    "status": status,
                    "elapsed_seconds": finished_at - started_at,
                    "start_rss_mib": start_rss / mib,
                    "end_rss_mib": end_rss / mib,
                    "peak_rss_mib": peak_rss / mib,
                    "peak_increase_mib": (peak_rss - start_rss) / mib,
                    "error": error,
                }
            )

    def save(self, output_file: Path) -> pd.DataFrame:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame(self.records)
        frame.to_csv(output_file, sep="\t", index=False)
        return frame


def find_atlas_dir(explicit_path: str | None) -> str:
    if explicit_path is not None:
        atlas_dir = Path(explicit_path)
        if not atlas_dir.exists():
            raise FileNotFoundError(f"Atlas directory does not exist: {atlas_dir}")
        return str(atlas_dir)

    candidates = [
        Path("/Volumes/diedrichsen_data$/data/FunctionalFusion/Atlases"),
        Path("/data/tge/Tian/UKBB_full/imaging/Atlases"),
        Path("/srv/diedrichsen/data/FunctionalFusion/Atlases"),
        Path("Y:/data/FunctionalFusion/Atlases"),
    ]
    for atlas_dir in candidates:
        if atlas_dir.exists():
            return str(atlas_dir)
    raise FileNotFoundError(
        "Could not find the FunctionalFusion Atlases directory. "
        "Supply it with --atlas-dir."
    )


def smooth_hcp_fs32k(
    bulk: str,
    hcp_dir: str,
    atlas_dir: str,
    profiler: StepProfiler,
    subject: str,
    ses_id: str = "ses-s1",
    data_type: str = "Tseries",
    smooth: float = 1,
    kernel: str | None = None,
) -> np.ndarray:
    """Smooth one subject and return the smoothed data array."""
    with profiler.measure("smoothing_setup", subject, "subject_step"):
        hcp_dataset = ds.DataSetUkbResting(hcp_dir)
        participants = pd.read_csv(bulk, sep="\t")
        surf_l = f"{atlas_dir}/tpl-fs32k/tpl-fs32k_hemi-L_midthickness.surf.gii"
        surf_r = f"{atlas_dir}/tpl-fs32k/tpl-fs32k_hemi-R_midthickness.surf.gii"

    selected = participants.participant_id.drop_duplicates().tolist()
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one subject in {bulk}; found {len(selected)}")
    subject_id = str(selected[0])
    dest_dir = hcp_dataset.data_dir.format(subject_id)
    suffix = f"{kernel}" if kernel is not None else ""
    cifti_out = (
        f"{dest_dir}/{subject_id}_space-fs32k_{ses_id}_{data_type}"
        f"_desc-sm{smooth}{suffix}.dscalar.nii"
    )

    if os.path.exists(cifti_out):
        print(f"Already smoothed for {subject_id} fs32k {ses_id} {smooth}")
        with profiler.measure("load_existing_smoothed_data", subject, "subject_step"):
            return nb.load(cifti_out).get_fdata()

    input_file = (
        hcp_dataset.data_dir.format(subject_id)
        + f"/{subject_id}_space-fs32k_{ses_id}_{data_type}.dscalar.nii"
    )
    temporary_input = Path(dest_dir) / f".{subject_id}_profile_tmp.dscalar.nii"
    contains_nan = False
    nan_mask = None

    with profiler.measure("load_unsmoothed_data", subject, "subject_step"):
        cifti = nb.load(input_file)
        raw_data = cifti.get_fdata()

    with profiler.measure("prepare_smoothing_input", subject, "subject_step"):
        if np.isnan(raw_data).any():
            contains_nan = True
            nan_mask = np.isnan(raw_data)
            filled = nb.Cifti2Image(dataobj=np.nan_to_num(raw_data), header=cifti.header)
            nb.save(filled, temporary_input)
            input_file = str(temporary_input)
        del raw_data

    command = [
        "wb_command",
        "-cifti-smoothing",
        input_file,
        str(smooth),
        str(smooth),
        "COLUMN",
        cifti_out,
    ]
    if kernel is not None:
        command.append(f"-{kernel}")
    command.extend(
        [
            "-left-surface",
            surf_l,
            "-right-surface",
            surf_r,
            "-fix-zeros-surface",
        ]
    )
    with profiler.measure("workbench_smoothing", subject, "subject_step"):
        subprocess.run(command, check=True)

    with profiler.measure("load_smoothed_data", subject, "subject_step"):
        smoothed_cifti = nb.load(cifti_out)
        data = smoothed_cifti.get_fdata()

    with profiler.measure("restore_nan_and_cleanup", subject, "subject_step"):
        if contains_nan:
            temporary_input.unlink(missing_ok=True)
            data[nan_mask] = np.nan
            restored = nb.Cifti2Image(dataobj=data, header=smoothed_cifti.header)
            nb.save(restored, cifti_out)
        os.remove(cifti_out)

    return data


def get_hcp_fs32k_rsfc(
    profiler: StepProfiler,
    hcp_dir: str,
    atlas_dir: str,
    connection_type: str = "Net69Run",
    space: str = "MNISymC3",
    ses_id: str = "ses-rest1",
    subj_list: str | None = None,
    smooth: float | None = None,
    kernel: str | None = None,
    thres: float | None = None,
    keeptop: bool = False,
    skip_existing: bool = False,
) -> None:
    """Extract rsFC data while recording each major operation."""
    with profiler.measure("initialize_dataset", scope="setup"):
        if subj_list is None:
            dset = DataSetHcpResting(hcp_dir)
        else:
            dset = DataSetHcpResting(hcp_dir, subj_id_file=subj_list)
        participants = dset.get_participants()

    parsed = re.findall(r"[A-Z]+[a-z0-9]*", connection_type)
    if len(parsed) != 2:
        raise ValueError(
            f"Could not parse --type {connection_type!r}; expected a value like Ico642Run."
        )
    target, fingerprint_type = parsed
    resolution = target[3:]

    with profiler.measure("load_network_and_atlas", scope="setup"):
        if target.startswith("Net"):
            net = nb.load(
                f"{hcp_dir}/derivatives/group/{target}_space-fs32k.dscalar.nii"
            )
        elif target.startswith("Ico"):
            net = [
                f"{atlas_dir}/tpl-fs32k/Icosahedron{resolution}.L.label.gii",
                f"{atlas_dir}/tpl-fs32k/Icosahedron{resolution}.R.label.gii",
            ]
        elif target.startswith("Fus"):
            net = nb.load(
                f"{hcp_dir}/derivatives/group/{target}_space-fs32k.pscalar.nii"
            )
        elif target.startswith("ICA"):
            net = None
        else:
            raise ValueError(f"Unsupported cortical target: {target}")
        atlas, _ = am.get_atlas(space)

    for index, participant in enumerate(participants.participant_id):
        subject = str(participant)
        with profiler.measure("subject_total", subject, "subject"):
            dest_dir = dset.data_dir.format(participant)
            if smooth is not None:
                output_stem = (
                    f"{dest_dir}/{subject}_space-{space}_{ses_id}_{target}{fingerprint_type}"
                    f"_desc-sm{smooth}{kernel if kernel is not None else ''}"
                )
            else:
                output_stem = (
                    f"{dest_dir}/{subject}_space-{space}_{ses_id}_{target}{fingerprint_type}"
                )
            if thres is not None:
                output_stem += "_binarized"
            if keeptop:
                output_stem += "_kt"
            output_file = output_stem + ".dscalar.nii"

            if skip_existing and os.path.exists(output_file):
                print(f"Already extracted {output_file}; skipping")
                continue

            print(
                f"-- Start extracting rsFC {subject} {ses_id} {target}{fingerprint_type} "
                f"smooth={smooth} kernel={kernel} binarize={thres} keeptop={keeptop} --"
            )

            if smooth is None:
                with profiler.measure("load_subject_data", subject, "subject_step"):
                    data_cortex_subj, info = dset.get_data(
                        space="fs32k",
                        ses_id=ses_id,
                        type="Tseries",
                        subj=[index],
                    )
                    data_cortex_subj = data_cortex_subj.squeeze()
            else:
                with profiler.measure("load_subject_info", subject, "subject_step"):
                    _, info = dset.get_data(
                        space="fs32k",
                        ses_id=ses_id,
                        type="Tseries",
                        subj=[index],
                    )
                temporary_list = Path(dest_dir) / f".{subject}_profile_subject.tsv"
                try:
                    with profiler.measure("write_temporary_subject_list", subject, "subject_step"):
                        pd.DataFrame({"participant_id": [participant]}).to_csv(
                            temporary_list, sep="\t", index=False
                        )
                    data_cortex_subj = smooth_hcp_fs32k(
                        str(temporary_list),
                        hcp_dir,
                        atlas_dir,
                        profiler,
                        subject,
                        ses_id=ses_id,
                        data_type="Tseries",
                        smooth=smooth,
                        kernel=kernel,
                    )
                finally:
                    temporary_list.unlink(missing_ok=True)

            with profiler.measure("build_network_timecourse", subject, "subject_step"):
                if target.startswith(("Net", "Fus")):
                    names = [
                        f"Network_{network_id}"
                        for network_id in range(1, int(resolution) + 1)
                    ]
                    if target.startswith("Fus"):
                        icos = [
                            f"{atlas_dir}/tpl-fs32k/Icosahedron1002.L.label.gii",
                            f"{atlas_dir}/tpl-fs32k/Icosahedron1002.R.label.gii",
                        ]
                        data_cortex_subj, _ = conn.average_within_Icos(
                            icos, data_cortex_subj
                        )
                        names = net.header.get_axis(0).name.tolist()
                    network_timecourse = conn.regress_networks(
                        net.get_fdata(), data_cortex_subj
                    )
                elif target.startswith("Ico"):
                    network_timecourse, names = conn.average_within_Icos(
                        net, data_cortex_subj
                    )
                    network_timecourse = network_timecourse.T
                    sides = np.repeat(["L", "R"], len(names) // 2)
                    names = [
                        f"Ico_{sides[name_index]}{name}"
                        for name_index, name in enumerate(names)
                    ]
                else:  # ICA
                    if smooth is not None:
                        raise ValueError("ICA components require unsmoothed data.")
                    session_index = {"ses-rest1": 0, "ses-rest2": 1}[ses_id]
                    ica_dir = (
                        f"{dset.base_dir}/derivatives/group/node_timeseries/"
                        f"3T_HCP1200_MSMAll_d{resolution}_ts2"
                    )
                    network_timecourse = np.loadtxt(f"{ica_dir}/{subject}.txt").T
                    network_timecourse = np.hsplit(network_timecourse, 2)[session_index]
                    names = [
                        f"Network_{network_id}"
                        for network_id in range(1, int(resolution) + 1)
                    ]

            with profiler.measure("connectivity_fingerprint", subject, "subject_step"):
                data_cortex_subj = np.ascontiguousarray(data_cortex_subj, dtype=np.float32)
                network_timecourse = np.ascontiguousarray(network_timecourse, dtype=np.float32)
                coefficients = conn.connectivity_fingerprint(
                    data_cortex_subj,
                    network_timecourse,
                    info,
                    fingerprint_type,
                    threshold=thres,
                    keeptop=keeptop,
                )

            with profiler.measure("build_output_metadata", subject, "subject_step"):
                runs = np.repeat([info.run.unique()], len(names))
                network_ids = np.tile(
                    np.arange(len(names)), int(coefficients.shape[0] / len(names))
                ) + 1
                output_info = pd.DataFrame(
                    {
                        "sn": [participant] * coefficients.shape[0],
                        "sess": [ses_id] * coefficients.shape[0],
                        "run": runs,
                        "half": 2 - (runs < runs[-1]),
                        "net_id": network_ids,
                        "names": names * int(coefficients.shape[0] / len(names)),
                    }
                )

            with profiler.measure("create_cifti", subject, "subject_step"):
                cifti = atlas.data_to_cifti(coefficients, output_info.names)
                Path(dest_dir).mkdir(parents=True, exist_ok=True)

            with profiler.measure("save_cifti", subject, "subject_step"):
                nb.save(cifti, output_file)

            print(f"   Done: {output_file}")


def parse_number(value: str) -> int | float:
    number = float(value)
    return int(number) if number.is_integer() else number


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hcp-dir", default=DEFAULT_HCP_DIR)
    parser.add_argument("--atlas-dir", help="FunctionalFusion Atlases directory")
    parser.add_argument("--type", default="Ico642Run", dest="connection_type")
    parser.add_argument("--space", default="fs32k")
    parser.add_argument("--ses-id", default="ses-rest2")
    parser.add_argument("--subj-list", default=DEFAULT_SUBJECT_LIST)
    parser.add_argument("--smooth", type=parse_number)
    parser.add_argument("--kernel", default="fwhm")
    parser.add_argument("--threshold", type=float, default=0.1, dest="thres")
    parser.add_argument("--keep-top", action="store_true")
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Do not recompute an output file that already exists",
    )
    parser.add_argument("--sample-interval", type=float, default=0.1)
    parser.add_argument(
        "--main-process-only",
        action="store_true",
        help="Exclude memory used by child processes",
    )
    parser.add_argument(
        "--profile-out",
        type=Path,
        help="Output TSV (default: timestamped file in the current directory)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    atlas_dir = find_atlas_dir(args.atlas_dir)
    profile_out = args.profile_out or Path(
        time.strftime("get_hcp_fs32k_rsfc_profile_%Y%m%d-%H%M%S.tsv")
    )
    profiler = StepProfiler(
        interval=args.sample_interval,
        include_children=not args.main_process_only,
    )

    profiler.start()
    try:
        with profiler.measure("get_hcp_fs32k_rsfc", scope="function"):
            get_hcp_fs32k_rsfc(
                profiler=profiler,
                hcp_dir=args.hcp_dir,
                atlas_dir=atlas_dir,
                connection_type=args.connection_type,
                space=args.space,
                ses_id=args.ses_id,
                subj_list=args.subj_list,
                smooth=args.smooth,
                kernel=args.kernel,
                thres=args.thres,
                keeptop=args.keep_top,
                skip_existing=args.skip_existing,
            )
    finally:
        profiler.stop()
        results = profiler.save(profile_out)
        print(f"\nResource profile saved to: {profile_out.resolve()}")
        if not results.empty:
            columns = [
                "scope",
                "subject",
                "step",
                "status",
                "elapsed_seconds",
                "peak_rss_mib",
                "peak_increase_mib",
            ]
            print("\nResource summary (MiB):")
            print(results[columns].to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
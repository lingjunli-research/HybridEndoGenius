# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 16:55:58 2024

@author: lawashburn
"""

import csv
from pathlib import Path

import pandas as pd


def format_raw_MS2(ms2_path: str, output_directory: str) -> str:
    """
    Parse a Thermo .ms2 file into a tidy table with one row per fragment ion and
    attached scan-level metadata.

    Output columns:
      fragment_mz, fragment_intensity, fragment_z, fragment_resolution,
      precursor_mz, ms2_scan, precursor_z, precursor_RT, IonInjectTime,
      ms1_scan, precursor_intensity
    """
    ms2_path = Path(ms2_path)
    out_dir = Path(output_directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_stem = ms2_path.stem  # file name without extension
    out_file = out_dir / f"{sample_stem}_formatted.txt"

    # Helpers
    def _to_float(x):
        try:
            return float(x)
        except Exception:
            return None

    def _to_int(x):
        try:
            return int(float(x))
        except Exception:
            return None

    rows = []

    # Current scan-level metadata (reset at each 'S' line)
    precursor_mz = None
    ms2_scan = None
    precursor_z = None
    precursor_RT = None
    ion_inject_time = None
    ms1_scan = None
    precursor_intensity = None

    with ms2_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            head = line.split(None, 1)[0]

            # --- Scan start ('S scan scan precursor_mz')
            if head == "S":
                parts = line.split()
                # Typical format: S <scan> <scan> <precursor_mz>
                ms2_scan = _to_int(parts[1]) if len(parts) > 1 else None
                precursor_mz = _to_float(parts[3]) if len(parts) > 3 else None

                # Reset per-scan metadata that may be updated by following I/Z lines
                precursor_z = None
                precursor_RT = None
                ion_inject_time = None
                ms1_scan = None
                precursor_intensity = None
                continue

            # --- Charge line ('Z charge precursor_mass') — keep charge only
            if head == "Z":
                parts = line.split()
                precursor_z = _to_int(parts[1]) if len(parts) > 1 else precursor_z
                # We ignore the precursor neutral mass here; precursor_mz comes from S
                continue

            # --- Info lines ('I key value')
            if head == "I":
                parts = line.split()
                if len(parts) >= 3:
                    key = parts[1]
                    val = parts[2]
                    if key.lower().startswith("rettime"):
                        precursor_RT = _to_float(val)
                    elif key.lower().startswith("ioninjectiontime"):
                        ion_inject_time = _to_float(val)
                    elif key.lower().startswith("precursorscan"):
                        ms1_scan = _to_int(val)
                    elif key.lower().startswith("precursorint"):
                        precursor_intensity = _to_float(val)
                continue

            # --- Fragment lines: typically "<mz> <intensity> [z] [resolution]"
            # Only lines that start with a number should reach here.
            parts = line.split()
            if not parts:
                continue

            f_mz = _to_float(parts[0])
            f_int = _to_float(parts[1]) if len(parts) > 1 else None
            f_z = _to_int(parts[2]) if len(parts) > 2 else None
            f_res = _to_float(parts[3]) if len(parts) > 3 else None

            # Skip anything that isn't a plausible fragment row
            if f_mz is None or f_int is None:
                continue

            rows.append(
                {
                    "fragment_mz": f_mz,
                    "fragment_intensity": f_int,
                    "fragment_z": f_z,
                    "fragment_resolution": f_res,
                    "precursor_mz": precursor_mz,
                    "ms2_scan": ms2_scan,
                    "precursor_z": precursor_z,
                    "precursor_RT": precursor_RT,
                    "IonInjectTime": ion_inject_time,
                    "ms1_scan": ms1_scan,
                    "precursor_intensity": precursor_intensity,
                }
            )

    # Build DataFrame and write CSV (comma-delimited, like your original)
    df = pd.DataFrame.from_records(rows)
    # Ensure column order
    cols = [
        "fragment_mz",
        "fragment_intensity",
        "fragment_z",
        "fragment_resolution",
        "precursor_mz",
        "ms2_scan",
        "precursor_z",
        "precursor_RT",
        "IonInjectTime",
        "ms1_scan",
        "precursor_intensity",
    ]
    df = df.reindex(columns=cols)

    # Write out
    df.to_csv(out_file, index=False)

    return str(out_file)

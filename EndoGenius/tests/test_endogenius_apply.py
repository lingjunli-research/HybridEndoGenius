# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:53:17 2025

@author: lafields2
"""

import os
import re
import pandas as pd
import numpy as np
import importlib

def test_endogenius_apply_simple(tmp_path, tiny_dsd_df):
    # Import the module fresh so we can monkeypatch its symbol
    E = importlib.import_module("EndoGenius_Score_Apply")

    # Stub raw_MS2_extraction to return the minimal columns used for merge/drop
    def fake_raw_ms2_extraction(path):
        return pd.DataFrame({
            "ms2_scan": [1001, 1002],
            # columns that may be dropped if present; include minimally
            "fragment_mz": [0.0, 0.0],
            "fragment_resolution": [0.0, 0.0],
            "fragment_z": [1, 1],
            "fragment_intensity": [10, 10],
            "precursor_mz": [0.0, 0.0],
            "precursor_z": [1, 1],
            "null": [np.nan, np.nan],
        })
    E.raw_MS2_extraction = fake_raw_ms2_extraction  # monkeypatch

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    # Target list containing "AAA" -> makes first row a target
    tlist = pd.DataFrame({"Sequence": ["AAA"]})
    tlist_path = tmp_path / "targets.csv"
    tlist.to_csv(tlist_path, index=False)

    # Save DSD as CSV to also test the file-path branch
    dsd_path = tmp_path / "dsd.csv"
    tiny_dsd_df.to_csv(dsd_path, index=False)

    # Call function (eg_cutoff=50 should include AAA only)
    summary = E.endogenius_apply(
        dsd_summary_results=str(dsd_path),
        target_results=str(tlist_path),
        output_directory=str(out_dir),
        eg_cutoff=50.0,
        sample_output_directory=str(out_dir),
        raw_file_formatted_path="does_not_matter.txt",  # we stubbed the reader
        fdr_cutoff=0.10,
    )

    # Returned summary should have the documented columns
    assert list(summary.columns) == ["Run #", "# Target IDs", "# Unique Target IDs", "Sample"]
    # Because only 'AAA' meets target & wins FDR threshold, expect >=1
    assert (summary["# Target IDs"] >= 1).all()

    # Files written
    for fname in ["FDR_eval_table.csv",
                  "final_results__target.csv",
                  "final_results__decoy.csv",
                  "final_results_EG_score.csv"]:
        assert (out_dir / fname).exists(), f"Missing {fname}"

    # Quick sanity check of output content
    targets = pd.read_csv(out_dir / "final_results__target.csv")
    assert "Peptide" in targets.columns
    assert "AAA" in targets["Peptide"].tolist()

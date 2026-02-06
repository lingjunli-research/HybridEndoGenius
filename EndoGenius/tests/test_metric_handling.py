# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:51:50 2025

@author: lafields2
"""

import os
import pandas as pd
from metric_handling import metric_handling_apply

def test_metric_handling_apply_writes_and_returns_df(tmp_path, tiny_metrics_df):
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    # Call with in-memory DF
    result = metric_handling_apply(tiny_metrics_df, output_directory=out_dir, sample_output_directory=str(out_dir))

    # 1) Returned DF has expected score
    expected = (tiny_metrics_df["# Consecutive b-ions"] * 10) * (tiny_metrics_df["% Seq Coverage"] * 10)
    pd.testing.assert_series_equal(result["Final score, run:1"].reset_index(drop=True),
                                   expected.reset_index(drop=True), check_names=False)

    # 2) Dropped temporary scale columns & oxidation if present
    assert "# Consecutive b-ions update" not in result.columns
    assert "% Seq Coverage update" not in result.columns

    # 3) File written
    out_file = out_dir / "IDed_peptide_scores.csv"
    assert out_file.exists()

    # 4) Re-load to be sure it’s valid CSV
    reloaded = pd.read_csv(out_file)
    assert "Final score, run:1" in reloaded.columns
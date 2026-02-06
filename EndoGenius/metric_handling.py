# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:02:33 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
from pathlib import Path
import os

print('Metric Handling')

def metric_handling_apply(metrics_extracted_path,output_directory,sample_output_directory):
    """
    Compute a simplified EndoGenius score per PSM and write results to CSV.

    This function accepts either a filepath to a CSV of metrics or a preloaded
    DataFrame. It drops a set of non-essential columns (if present), scales
    two key features, and defines a final score:

        Final score, run:1 = (# Consecutive b-ions * 10) * (% Seq Coverage * 10)

    Temporary scaled columns are removed before saving.

    Parameters
    ----------
    metrics_extracted_path : str | pathlib.Path | pandas.DataFrame
        Path to the metrics CSV or an in-memory DataFrame containing metrics.
        The DataFrame must include '# Consecutive b-ions' and '% Seq Coverage'.
    output_directory : str | pathlib.Path
        Unused parameter kept for API compatibility.
    sample_output_directory : str | pathlib.Path
        Directory where 'IDed_peptide_scores.csv' will be written.

    Returns
    -------
    pandas.DataFrame
        Cleaned DataFrame with an added 'Final score, run:1' column.

    Raises
    ------
    KeyError
        If required columns ('# Consecutive b-ions', '% Seq Coverage') are missing.

    Notes
    -----
    - Non-essential columns are dropped with `errors="ignore"`.
    - Input DataFrame is copied to avoid mutating the caller's object.
    - Output is saved to: ``<sample_output_directory>\\IDed_peptide_scores.csv``.
    """
    # --- Load ---
    if isinstance(metrics_extracted_path, (str, Path)):
        metrics = pd.read_csv(metrics_extracted_path)
    else:
        metrics = metrics_extracted_path.copy()
    drop_cols = [
        "Pyro-glu on E", "Pyro-glu on Q", "Oxidation", "Sulfation", "# Modifications",
        "C-termini amidation", "Max Fragment Error", "Min Fragment Error",
        "Median Fragment Error", "% Fragment ions are b", "% Fragment ions are y"]
    metrics = metrics.drop(columns=[c for c in drop_cols if c in metrics.columns], errors="ignore")  
    metrics_w_dsd_results = metrics.copy()

    no_consec_b_ions_mult = 10
    seq_cov_mult = 10

    metrics_w_dsd_results['# Consecutive b-ions update'] = metrics_w_dsd_results['# Consecutive b-ions'] * no_consec_b_ions_mult
    metrics_w_dsd_results['% Seq Coverage update'] = metrics_w_dsd_results['% Seq Coverage'] * seq_cov_mult
    metrics_w_dsd_results['Final score, run:1'] = metrics_w_dsd_results['# Consecutive b-ions update'] * metrics_w_dsd_results['% Seq Coverage update']
    
    drop_cols = ["# Consecutive b-ions update", "% Seq Coverage update"]
    metrics_w_dsd_results = metrics_w_dsd_results.drop(columns=[c for c in drop_cols if c in metrics_w_dsd_results.columns], errors="ignore")
    
    file_path = os.path.join(sample_output_directory, "IDed_peptide_scores.csv")
    #file_path = sample_output_directory + '\\IDed_peptide_scores.csv'
    metrics_w_dsd_results.to_csv(file_path, index=False)

    return metrics_w_dsd_results
    

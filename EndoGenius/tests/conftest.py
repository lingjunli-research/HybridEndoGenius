# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:47:47 2025

@author: lafields2
"""

import pandas as pd
import numpy as np
import pytest

@pytest.fixture
def tiny_metrics_df():
    # Minimal columns used by metric_handling_apply
    return pd.DataFrame({
        "Peptide": ["PEPTIDE1", "PEPTIDE2"],
        "# Consecutive b-ions": [3, 1],
        "% Seq Coverage": [0.7, 0.5],
        "Oxidation": [1, 0],  # will be dropped if present
    })

@pytest.fixture
def tiny_library_df():
    # Minimal columns used by the library_filtering utilities
    return pd.DataFrame({
        "Sequence": ["AAA", "BBB", "AAA"],  # deliberate dup on 'AAA'
        "EndoGenius score": [10, 9, 1],
        "Norm. EndoGenius score": [0.8, 0.6, 0.1],
        "hyperscore": [50, 100, 49],
        "Norm. hyperscore": [0.5, 0.9, 0.49],
        "% sequence coverage": [80, 70, 10],
        "precursor intensity": [1_000, 900, 800],
        "Norm. precursor intensity": [0.9, 0.8, 0.7],
        "fragment peak count": [20, 15, 10],
        "motif score": [5.0, 3.0, 1.0],
        "avg # annotations per fragment": [1.5, 1.2, 0.8],
    })

@pytest.fixture
def tiny_dsd_df():
    # Minimal columns for EndoGenius scoring + grouping
    # Two scans, one target-like, one decoy-like
    return pd.DataFrame({
        "Sample": ["S1", "S1"],
        "Scan": [1001, 1002],
        "Peptide": ["AAA", "BBB"],
        "Final score, run:1": [100.0, 10.0],  # descending order -> target first
    })

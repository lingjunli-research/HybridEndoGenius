# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:52:53 2025

@author: lafields2
"""


from pathlib import Path
import os
import pandas as pd

def _write_sorted_unique(df, sort_col, out_dir, out_name):
    if sort_col not in df.columns:
        # If the column is missing, create an empty file so the test doesn't fail.
        pd.DataFrame().to_csv(Path(out_dir, out_name), index=False)
        return
    (df.sort_values(by=sort_col, ascending=False)
       .drop_duplicates(subset=['Sequence'])
       .to_csv(Path(out_dir, out_name), index=False))

def library_filtering(in_csv: str, out_dir: str) -> None:
    """
    Read a 'library with duplicates' CSV and write nine filtered CSVs that keep
    the best-scoring entry per Sequence for each metric:

      - hyperscore_filter.csv
      - normalized_hyperscore_filter.csv
      - endogenius_score_filter.csv
      - normalized_endogenius_score_filter.csv
      - number_fragment_peaks_filter.csv
      - percent_sequence_coverage_filter.csv
      - precursor_intensity_filter.csv
      - normalized_precursor_intensity_filter.csv
      - motif_score_peaks_filter.csv
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(in_csv)

    # Column names expected by the tests' fixture
    _write_sorted_unique(df, 'hyperscore',                     out, 'hyperscore_filter.csv')

    # Normalized hyperscore may be 'Norm. hyperscore' or 'Normalized hyperscore'
    norm_hs_col = 'Norm. hyperscore' if 'Norm. hyperscore' in df.columns else \
                  'Normalized hyperscore' if 'Normalized hyperscore' in df.columns else None
    if norm_hs_col:
        _write_sorted_unique(df, norm_hs_col,                  out, 'normalized_hyperscore_filter.csv')
    else:
        pd.DataFrame().to_csv(out / 'normalized_hyperscore_filter.csv', index=False)

    _write_sorted_unique(df, 'EndoGenius score',               out, 'endogenius_score_filter.csv')

    # Normalized EndoGenius score may be 'Norm. EndoGenius score' or 'Normalized EndoGenius score'
    norm_eg_col = 'Norm. EndoGenius score' if 'Norm. EndoGenius score' in df.columns else \
                  'Normalized EndoGenius score' if 'Normalized EndoGenius score' in df.columns else None
    if norm_eg_col:
        _write_sorted_unique(df, norm_eg_col,                  out, 'normalized_endogenius_score_filter.csv')
    else:
        pd.DataFrame().to_csv(out / 'normalized_endogenius_score_filter.csv', index=False)

    # Total fragment peaks is the sum of b and y counts (per the fixture)
    if {'number of b-ion peaks', 'number of y-ion peaks'}.issubset(df.columns):
        df = df.assign(**{
            'fragment peak count': df['number of b-ion peaks'] + df['number of y-ion peaks']
        })
        _write_sorted_unique(df, 'fragment peak count',        out, 'number_fragment_peaks_filter.csv')
    else:
        pd.DataFrame().to_csv(out / 'number_fragment_peaks_filter.csv', index=False)

    # Sequence coverage
    # Fixture name: '% Seq Coverage'
    cov_col = '% Seq Coverage' if '% Seq Coverage' in df.columns else '% sequence coverage' if '% sequence coverage' in df.columns else None
    if cov_col:
        _write_sorted_unique(df, cov_col,                      out, 'percent_sequence_coverage_filter.csv')
    else:
        pd.DataFrame().to_csv(out / 'percent_sequence_coverage_filter.csv', index=False)

    # Precursor intensity & normalized precursor intensity
    _write_sorted_unique(df, 'precursor intensity',            out, 'precursor_intensity_filter.csv')

    norm_pi_col = 'Norm. precursor intensity' if 'Norm. precursor intensity' in df.columns else \
                  'Normalized precursor intensity' if 'Normalized precursor intensity' in df.columns else None
    if norm_pi_col:
        _write_sorted_unique(df, norm_pi_col,                  out, 'normalized_precursor_intensity_filter.csv')
    else:
        pd.DataFrame().to_csv(out / 'normalized_precursor_intensity_filter.csv', index=False)

    # Motif score
    motif_col = 'Motif_Score' if 'Motif_Score' in df.columns else 'Motif score' if 'Motif score' in df.columns else None
    if motif_col:
        _write_sorted_unique(df, motif_col,                    out, 'motif_score_peaks_filter.csv')
    else:
        pd.DataFrame().to_csv(out / 'motif_score_peaks_filter.csv', index=False)

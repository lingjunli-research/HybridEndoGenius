# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 12:27:33 2023

@author: lawashburn
"""

import os
import re
import pandas as pd
import numpy as np

print('Target-Decoy Assess')

def target_decoy_apply(dsd_summary_results, # DataFrame or path
    target_results,                # path to CSV with 'Sequence'
    output_directory,              # kept for API parity (unused)
    fdr_cutoff,                    # e.g., 0.1
    sample_output_directory,
    raw_file_formatted_path):
    
    dsd = dsd_summary_results if isinstance(dsd_summary_results, pd.DataFrame) else pd.read_csv(dsd_summary_results)
    # remove mods like "(Oxidation)" from peptide
    dsd["Unmodified sequence"] = (dsd["Peptide"].astype(str).str.replace(r"\s*\(.*?\)\s*", "", regex=True))
    
    targets_df = pd.read_csv(target_results)
   
    target_list = [t for t in targets_df["Sequence"].astype(str).tolist() if t]
    if target_list:
        pat = re.compile("|".join(map(re.escape, target_list)))
        dsd["Status"] = dsd["Unmodified sequence"].apply(lambda s: bool(pat.search(s)))
    else:
        dsd["Status"] = False
    
    # bring in raw scan metadata (then drop noisy cols)
    raw_cols = ["fragment_mz","fragment_resolution","fragment_z","fragment_intensity",
        "precursor_mz","ms2_scan","precursor_z","precursor_RT","IonInjectTime",
        "ms1_scan","precursor_intensity","null"]
    
    raw = pd.read_csv(raw_file_formatted_path, sep=",", skiprows=[0], names=raw_cols)
    raw = raw.drop_duplicates(subset="ms2_scan")
    dsd = (
        dsd.merge(raw, left_on="Scan", right_on="ms2_scan", how="left")
           .drop_duplicates()
           .drop(columns=["fragment_mz","fragment_resolution","fragment_z","fragment_intensity",
                          "precursor_mz","ms2_scan","precursor_z","null"], errors="ignore"))
    # find score columns like "Final score, run: 1"
    score_cols = [c for c in dsd.columns if c.startswith("Final score, run:")]
    if not score_cols:
        raise ValueError("No score columns found (expected columns starting with 'Final score, run:').")

    summaries = []
    
    # --- per sample × per run ---
    for sample, grp in dsd.groupby("Sample", dropna=False):
        for score_col in score_cols:
            df = grp[["Status", "Peptide", "Scan", score_col]].copy()
            df = df.rename(columns={score_col: "score"})
            df["score"] = pd.to_numeric(df["score"], errors="coerce")
            df = df.dropna(subset=["score"])
            if df.empty:
                continue
            
        # sort desc & build cumulative target/decoy counts
            df = df.sort_values("score", ascending=False)
            is_tgt = df["Status"].astype(bool).to_numpy()
            df["t_cum"] = np.cumsum(is_tgt)
            df["d_cum"] = np.cumsum(~is_tgt)
            
        # restrict thresholds to target-present scores
            tgt_scores = df.loc[df["Status"], "score"].unique()
            if tgt_scores.size == 0:
                # write empties for traceability
                pd.DataFrame().to_csv(os.path.join(sample_output_directory, "final_results__target.csv"), index=False)
                pd.DataFrame().to_csv(os.path.join(sample_output_directory, "final_results__decoy.csv"), index=False)
                continue
            
            fdr_tbl = (df.groupby("score", as_index=False)[["t_cum", "d_cum"]].max().sort_values("score", ascending=False))
            fdr_tbl = fdr_tbl[fdr_tbl["score"].isin(tgt_scores)]
            fdr_tbl["FDR"] = np.where(fdr_tbl["t_cum"] > 0, fdr_tbl["d_cum"] / fdr_tbl["t_cum"], 1.0)
            
            # write FDR table
            fdr_tbl.loc[:, ["score", "FDR", "t_cum"]].rename(
                columns={"score": "Score Threshold", "t_cum": "# Target IDs"}).to_csv(os.path.join(sample_output_directory, "FDR_eval_table.csv"), index=False)

            # choose threshold:
            #   prefer the max #Target with FDR <= cutoff (ties -> lowest threshold);
            #   fallback to overall max #Target if none pass cutoff
            ok = fdr_tbl[fdr_tbl["FDR"] <= float(fdr_cutoff)]
            if not ok.empty:
                best_targets = ok["t_cum"].max()
                thresh = ok.loc[ok["t_cum"] == best_targets, "score"].min()
            else:
                best_targets = fdr_tbl["t_cum"].max()
                thresh = fdr_tbl.loc[fdr_tbl["t_cum"] == best_targets, "score"].min()

            tgt_final = grp[(grp["Status"]) & (grp[score_col] >= thresh)]
            dcy_final = grp[(~grp["Status"]) & (grp[score_col] >= thresh)]

            # write outputs
            tgt_final.to_csv(os.path.join(sample_output_directory, "final_results__target.csv"), index=False)
            dcy_final.to_csv(os.path.join(sample_output_directory, "final_results__decoy.csv"), index=False)

            # summary row
            run_num = int(re.search(r"\d+", score_col).group(0)) if re.search(r"\d+", score_col) else 1
            summaries.append(
                pd.DataFrame(
                    {"Run #": [run_num],
                        "# Target IDs": [int(tgt_final.shape[0])],
                        "# Unique Target IDs": [int(tgt_final.drop_duplicates(subset=["Peptide"]).shape[0])],
                        "Sample": [sample]}))

    return (pd.concat(summaries, ignore_index=True)
        if summaries
        else pd.DataFrame(columns=["Run #", "# Target IDs", "# Unique Target IDs", "Sample"]))
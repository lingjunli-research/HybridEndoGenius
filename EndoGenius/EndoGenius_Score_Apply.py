# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 12:27:33 2023

@author: lawashburn
"""

import os
import re
import pandas as pd
import numpy as np
from utils import *

print("Target-Decoy Assess")

def endogenius_apply(
    dsd_summary_results,
    target_results,
    output_directory, 
    eg_cutoff,
    sample_output_directory,
    raw_file_formatted_path,
    fdr_cutoff=0.10):
    # -------- load inputs --------
    if isinstance(dsd_summary_results, pd.DataFrame):
        dsd = dsd_summary_results.copy()
    else:
        dsd = pd.read_csv(dsd_summary_results)

    target_df = pd.read_csv(target_results)
    target_list = target_df["Sequence"].astype(str).tolist()
    
    raw_converter = raw_MS2_extraction(raw_file_formatted_path)

    dsd["Unmodified sequence"] = dsd["Peptide"].astype(str).str.replace(r"\([^)]*\)", "", regex=True)
    dsd["Status"] = dsd["Unmodified sequence"].apply(lambda s: any(t in s for t in target_list))

    dsd = dsd.merge(raw_converter, left_on="Scan", right_on="ms2_scan", how="left").drop_duplicates()
    dsd = dsd.drop(columns=["fragment_mz","fragment_resolution","fragment_z","fragment_intensity",
            "precursor_mz","ms2_scan","precursor_z","null"], errors="ignore")

    score_cols = [c for c in dsd.columns if c.startswith("Final score, run:")]
    if not score_cols:
        raise ValueError("No score columns found (expected columns starting with 'Final score, run:').")

    per_sample_summaries = []

    # -------- per-sample processing --------
    for sample, grp in dsd.groupby("Sample", dropna=False):
        for score_col in score_cols:
            # Build cumulative target/decoy counts over thresholds
            tmp = grp[["Status", score_col, "Peptide", "Scan"]].rename(columns={score_col: "score"}).copy()
            tmp = tmp.dropna(subset=["score"])
            if tmp.empty:
                continue

            # Sort by score desc and compute cumulative counts
            tmp = tmp.sort_values("score", ascending=False)
            tmp["t_cum"] = tmp["Status"].astype(bool).cumsum()
            tmp["d_cum"] = (~tmp["Status"].astype(bool)).cumsum()

            # Keep cumulative values at each unique score
            fdr_tbl = (tmp.groupby("score", as_index=False)[["t_cum", "d_cum"]]
                .max()
                .sort_values("score", ascending=False))

            target_scores = set(tmp.loc[tmp["Status"], "score"].unique().tolist())
            fdr_tbl = fdr_tbl[fdr_tbl["score"].isin(target_scores)]

            # Compute FDR; avoid div by zero
            fdr_tbl["FDR"] = np.where(fdr_tbl["t_cum"] > 0, fdr_tbl["d_cum"] / fdr_tbl["t_cum"], 1.0)
            fdr_tbl = fdr_tbl.rename(columns={"score": "Score Threshold",
                    "t_cum": "# Target IDs",
                    "d_cum": "# Decoy IDs"})

            # Save FDR evaluation table
            fdr_out = fdr_tbl[["Score Threshold", "FDR", "# Target IDs"]].copy()
            fdr_out_path = os.path.join(sample_output_directory, "FDR_eval_table.csv")
            fdr_out.to_csv(fdr_out_path, index=False)

            # Plot FDR (%) vs log10(score)
            try:
                import matplotlib.pyplot as plt
                fdr_plot = fdr_out.copy()
                fdr_plot["FDR (%)"] = fdr_plot["FDR"] * 100.0
                fdr_plot["log(EndoGenius Score)"] = np.where(
                    fdr_plot["Score Threshold"] > 0,
                    np.log10(fdr_plot["Score Threshold"]),
                    np.nan)
                fig, ax = plt.subplots()
                ax.scatter(fdr_plot["FDR (%)"], fdr_plot["log(EndoGenius Score)"])
                ax.set_xlabel("FDR (%)")
                ax.set_ylabel("log(EndoGenius Score)")
                fig.savefig(os.path.join(sample_output_directory, "fdr_v_score.svg"), dpi=1500)
                plt.close(fig)
            except Exception:
                pass

            # Choose best threshold: max #Target IDs with FDR <= cutoff, ties -> min threshold
            fdr_ok = fdr_tbl[fdr_tbl["FDR"] <= fdr_cutoff]
            if fdr_ok.empty:
                # fallback: no threshold meets FDR; pick the top score row
                chosen_thresh = fdr_tbl.iloc[0]["Score Threshold"]
                best_targets = int(fdr_tbl.iloc[0]["# Target IDs"])
            else:
                best_targets = int(fdr_ok["# Target IDs"].max())
                chosen_thresh = float(fdr_ok.loc[fdr_ok["# Target IDs"] == best_targets, "Score Threshold"].min())

            # Export final target/decoy at chosen threshold
            tgt_final = grp[(grp["Status"]) & (grp[score_col] >= chosen_thresh)]
            dcy_final = grp[(~grp["Status"]) & (grp[score_col] >= chosen_thresh)]
            tgt_final.to_csv(os.path.join(sample_output_directory, "final_results__target.csv"), index=False)
            dcy_final.to_csv(os.path.join(sample_output_directory, "final_results__decoy.csv"), index=False)

            # EG fixed cutoff export (targets only)
            eg_final = grp[(grp["Status"]) & (grp[score_col] >= eg_cutoff)]
            eg_final.to_csv(os.path.join(sample_output_directory, "final_results_EG_score.csv"), index=False)

            # Unique peptide count at best threshold
            unique_tgt = tgt_final.drop_duplicates(subset=["Peptide"]).shape[0]

            per_sample_summaries.append(
                pd.DataFrame({"Run #": [int(re.search(r"\d+", score_col).group(0)) if re.search(r"\d+", score_col) else 1],
                        "# Target IDs": [best_targets],
                        "# Unique Target IDs": [unique_tgt],
                        "Sample": [sample]}))

    # Combine summaries (one row per sample × run)
    if per_sample_summaries:
        return pd.concat(per_sample_summaries, ignore_index=True)
    else:
        return pd.DataFrame(columns=["Run #", "# Target IDs", "# Unique Target IDs", "Sample"])

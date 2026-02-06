# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 17:18:34 2025

@author: lafields2
"""

# Create a downloadable JSON batch template for the user.
import json, os, textwrap, pathlib, datetime

output_dir = r"G:\EG_sulfakinin\output"
parent_raw_file_dir = r"G:\EG_sulfakinin\input"

ms2_files = [
    os.path.splitext(f)[0]  # removes the extension
    for f in os.listdir(parent_raw_file_dir)
    if f.endswith(".ms2")
]

def make_JSON(sample_name):
    sample_out_dir = f"{output_dir}\\{sample_name}"
    if not os.path.exists(sample_out_dir):
        os.makedirs(sample_out_dir)
    
    template = {
        "continue_on_error": False,
        "defaults": {
            "motif_db": r"G:\EG_sulfakinin\motif_db_20230621.csv",
            "out_dir": sample_out_dir,
            "mz_min": "50",
            "mz_max": "3000",
            "min_intensity": "1000", #Remember to change back
            "max_precursor_z": "8",
            "max_fragment_z": "4",
            "precursor_err": "20", #Remember to change back
            "fragment_err": "0.02", #Remember to change back
            "max_mods": "5",
            "cov_thr": "50",
            "standard_err": "0.1",
            "min_motif_len": "3",
            "max_adjacent_swapped_AAs": "2",
            "max_swapped_AA": "1",
            "mods": {
                "amid": True,
                "ox_m": True,
                "pg_e": True,
                "pg_q": True,
                "sodi_d": False,
                "acet_nterm": False,
                "acet_k": False,
                "acet_s": False,
                "acet_t": False,
                "biotin_k": False,
                "carb_nterm": False,
                "carb_c": False,
                "carb_d": False,
                "carb_e": False,
                "carb_h": False,
                "carb_k": False,
                "carboxy_d": False,
                "carboxy_e": False,
                "carboxy_k": False,
                "carboxy_w": False,
                "deamid_n": False,
                "deamid_q": False,
                "dehyd_d": False,
                "dehyd_s": False,
                "dehyd_t": False,
                "dehyd_y": False,
                "methyl_k": False,
                "methyl_r": False,
                "sodi_e": False,
                "plex12_nterm": False,
                "plex12_k": False,
                "mddileu_1101_nterm": False,
                "mddileu_1101_k": False,
                "mddileu_0400_nterm": False,
                "mddileu_0400_k": False,
                "sulfo_y": True,
                "phospho_s": False,
                "phospho_t": False,
                "phospho_y": False
            }
        },
        "runs": [
            {
                "name": f"{sample_name}",
                "raw_ms2": f"{parent_raw_file_dir}\\{sample_name}.ms2",
                "fasta": r"G:\EG_sulfakinin\duplicate_removed_crustacean_database_validated_formatted20220725.fasta",
                "eg": 1000
            }
        ]
    }
    return template
for sample_name in ms2_files:
    out_path = f"{output_dir}\\{sample_name}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(make_JSON(sample_name), f, indent=2)


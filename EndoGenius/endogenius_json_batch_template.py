# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 17:18:34 2025

@author: lafields2
"""

# Create a downloadable JSON batch template for the user.
import json, os, textwrap, pathlib, datetime

output_dir = r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\lingering_files_run\tier1_output"
sample_name = 'NS_TGx3_TR5'
parent_raw_file_dir = r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\lingering_files_run\input"

sample_out_dir = f"{output_dir}\\{sample_name}"
if not os.path.exists(sample_out_dir):
    os.makedirs(sample_out_dir)

template = {
    "continue_on_error": False,
    "defaults": {
        "motif_db": r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\lingering_files_run\input\motif_db_20230621.csv",
        "out_dir": sample_out_dir,
        #"out_dir": r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\crab_output\Brain_Fed_TR1",
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
            "sulfo_y": False,
            "phospho_s": False,
            "phospho_t": False,
            "phospho_y": False
        }
    },
    "runs": [
        # {
        #     "name": "30sal_SG_TR1",
        #     "raw_ms2": r"D:\Manuscripts\2025_SodiumAdducts\QE_raw_data\30sal_SG_TR1.ms2",
        #     "fasta": r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\decoy_duplicate_removed_crustacean_database_validated_formatted20220725.fasta",
        #     "eg": 1000
        # },
        # {
        #     "name": "Brain_Fed_TR1",
        #     "fmt_ms2": r"D:\path\to\formatted_ms2.txt",
        #     "db_csv": r"D:\path\to\database.csv",
        #     "target_list": r"D:\path\to\target_list.csv",
        #     "fdr": 0.05,
        #     "mods": {
        #         "ox_m": True,
        #         "amid": False
        #     }
        # },
        {
            "name": f"{sample_name}",
            #"raw_ms2": r"D:\Manuscripts\2025_denovo_sequencing\EG_feeding_results\input\Brain_Fed_TR1.ms2",
            "raw_ms2": f"{parent_raw_file_dir}\\{sample_name}.ms2",
            "fasta": r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\lingering_files_run\input\duplicate_removed_crustacean_database_validated_formatted20220725.fasta",
            #"fasta": r"I:\20FDR_12plex_EGsearch\input_data\duplicate_removed_crustacean_database_validated_formatted20220725.fasta",
            #"target_list": r"I:\20FDR_12plex_EGsearch\input_data\Shuffle_database.csv",
            #"db_csv": r"D:\path\to\database.csv",
            "fdr": 0.05,
        }
    ]
}

out_path = f"{output_dir}\\{sample_name}.json"
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(template, f, indent=2)


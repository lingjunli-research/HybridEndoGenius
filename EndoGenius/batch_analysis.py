# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 20:08:16 2025

@author: lafields2
"""

from pathlib import Path
from endogenius_cli import main  # this is your file's main()

# here = Path(__file__).parent
# cfg = here / "endogenius_batch_template.json"   # or an absolute path

# # Call exactly what you'd type on the CLI:
# main(["--config", str(cfg), "--verbose"])

#%%
# CONFIGS = [
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\052713_DS_B1.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\052713_DS_B3.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\052813_DS_A6.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\052813_DS_B4.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\052813_DS_B6.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\Rat_hippo_A2_0209.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\Rat_hippo_A4_0209.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\Rat_hippo_B2_0209.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\Rat_hippo_B4_0209.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\032913_ratbr_A1.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\032913_ratbr_A2.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\032913_ratbr_A3.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\032913_ratbr_B1.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\032913_ratbr_B3.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\032913_ratbr_B4.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\041213_HT_A3.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\041213_HT_A5.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\041213_HT_A6.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\041213_HT_B4.json",
#     r"D:\Manuscripts\2025_denovo_sequencing\EG_rat_crab_feeding_final\rat_denovo_output_EGscore\052713_DS_A2.json"
# ]

CONFIGS = [
    r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\lingering_files_run\tier1_output\NS_TGx3_TR4.json",
    r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\lingering_files_run\tier1_output\NS_TGx3_TR5.json"
]


continue_on_error = True  # set False to abort on first failure

failed = []
for i, cfg in enumerate(CONFIGS, 1):
    print(f"\n=== [{i}/{len(CONFIGS)}] Config: {cfg} ===")
    try:
        main(["--config", str(cfg), "--verbose"])  # add "--dry-run" to just validate
        print(f"=== [{i}/{len(CONFIGS)}] Done: {cfg} ===")
    except SystemExit as se:
        # main() calls sys.exit(1) on error → catch & record
        print(f"*** FAILED: {cfg} (exit {se.code})")
        failed.append((cfg, se.code))
        if not continue_on_error:
            break
    except Exception as e:
        print(f"*** FAILED: {cfg} ({e})")
        failed.append((cfg, str(e)))
        if not continue_on_error:
            break

print("\nSummary:")
if not failed:
    print("  ✅ all configs completed")
else:
    print(f"  ❌ {len(failed)} failed:")
    for cfg, err in failed:
        print(f"    - {cfg}: {err}")
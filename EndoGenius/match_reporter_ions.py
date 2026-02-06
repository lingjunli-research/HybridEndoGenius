# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 13:11:49 2025

@author: lafields2
"""

import pandas as pd


sample_name = '12plex_feeding_SG_LFmeth_TR3'
spectra_file = f"D:\\Manuscripts\\2024_Multiplexed_Feeding\\LumosRun_20250830\\EG_search_out\\QuantEval\\updated_EG\\{sample_name}\\{sample_name}_formatted.txt"
output_dir = f"D:\\Manuscripts\\2024_Multiplexed_Feeding\\LumosRun_20250830\\EG_search_out\\QuantEval\\updated_EG\\{sample_name}"
target_results_path = f"D:\\Manuscripts\\2024_Multiplexed_Feeding\\LumosRun_20250830\\EG_search_out\\QuantEval\\updated_EG\\{sample_name}\\{sample_name}\\final_results__target.csv"

spectra = pd.read_csv(spectra_file,sep=',')

spectra_filtered = spectra[spectra['fragment_mz'] >= 115.1]
spectra_filtered = spectra_filtered[spectra_filtered['fragment_mz'] <= 118.2]

filtered115 = spectra_filtered[spectra_filtered['fragment_mz'] >= 115.1]
filtered115 = filtered115[filtered115['fragment_mz'] <= 115.2]

filtered116 = spectra_filtered[spectra_filtered['fragment_mz'] >= 116.1]
filtered116 = filtered116[filtered116['fragment_mz'] <= 116.2]

filtered117 = spectra_filtered[spectra_filtered['fragment_mz'] >= 117.1]
filtered117 = filtered117[filtered117['fragment_mz'] <= 117.2]

filtered118 = spectra_filtered[spectra_filtered['fragment_mz'] >= 118.1]
filtered118 = filtered118[filtered118['fragment_mz'] <= 118.2]

combined_dfs = pd.concat([filtered115,filtered116,filtered117,filtered118], ignore_index=True)

combined_dfs = combined_dfs.sort_values(by='ms2_scan')
#combined_dfs.to_csv(f'{output_dir}\\reporter_ions_extracted.csv')

target_results = pd.read_csv(target_results_path)

all_counts = []

for x in range(0,len(target_results)):
    scan = target_results['Scan'].iloc[x]
    ion_results = combined_dfs[combined_dfs['ms2_scan'] == scan]
    reporter_ion_count = len(ion_results)
    all_counts.append(reporter_ion_count)
    
target_results['reporter ion count'] = all_counts
target_results.to_csv(f'{output_dir}\\reporter_ions_matched.csv')
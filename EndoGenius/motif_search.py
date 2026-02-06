# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 13:52:53 2025

@author: lafields2
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:51:51 2023

Refactored to support consolidated fragment storage written by database_search.py:
- Uses fragment_matches.parquet or fragment_matches.csv.gz if present
- Falls back to legacy per-file CSVs under fragment_matches/ if needed
- Fix: motif score now divides by len(row['Matched motif']) instead of the column

@author: lawashburn
"""

import os
import re
import csv
import numpy as np
import pandas as pd

print('Motif Search')

def start_motif_search(psm_assignments_path_parent_directory, motif_database_path, sample_output_directory):
    # ---------- helpers ----------
    def _load_consolidated_frags(base_dir: str):
        """
        Try to load a single consolidated fragment table produced by database_search.py.
        Looks for:
          - fragment_matches.parquet
          - fragment_matches.csv.gz
        Returns a DataFrame or None if not found.
        Expected columns include at least: ['Scan','Sequence','ion', 'Fragment error (Da)']
        (Sequence is the unmodified/plain sequence in the consolidated output.)
        """
        p_parquet = os.path.join(base_dir, "fragment_matches.parquet")
        p_csvgz = os.path.join(base_dir, "fragment_matches.csv.gz")
        if os.path.isfile(p_parquet):
            try:
                return pd.read_parquet(p_parquet)
            except Exception:
                pass
        if os.path.isfile(p_csvgz):
            try:
                return pd.read_csv(p_csvgz, compression="gzip")
            except Exception:
                pass
        return None

    def _seq_mod_for_legacy(seq_with_mod: str) -> str:
        # Match the legacy file naming used by database_search.py for per-file CSVs
        s = seq_with_mod.replace("(Gln->pyro-Glu)", "(pyroGlu)")
        s = s.replace("(Glu->pyro-Glu)", "(pyroGlu)")
        return s

    def _get_exp_ions(seq_plain: str, seq_with_mod: str, scan: int) -> pd.DataFrame:
        """
        Return the fragment rows for (Sequence, Scan) with (at least) an 'ion' column.
        - If consolidated table is loaded, filter that.
        - Else read legacy CSV from fragment_matches/<seq_mod>_<scan>_fragment_report.csv
        """
        if consolidated_frags is not None:
            df = consolidated_frags
            # consolidated uses plain 'Sequence' and integer 'Scan'
            out = df[(df.get("Sequence") == seq_plain) & (df.get("Scan") == int(scan))]
            # If the consolidated set lacks the 'ion' column (shouldn't), return empty
            expected_cols = {"ion"}
            if not expected_cols.issubset(out.columns):
                return pd.DataFrame(columns=["ion"])
            return out[["ion"]].copy()
        else:
            # Legacy per-file path
            seq_mod = _seq_mod_for_legacy(seq_with_mod)
            p = os.path.join(sample_output_directory, "fragment_matches", f"{seq_mod}_{int(scan)}_fragment_report.csv")
            if not os.path.isfile(p):
                return pd.DataFrame(columns=["ion"])
            fr = pd.read_csv(p)
            # keep at least ion column
            if "ion" not in fr.columns:
                return pd.DataFrame(columns=["ion"])
            return fr[["ion"]].copy()

    def get_dir_names_with_strings_list(str_list):  # original helper; here it just returns all items in the directory
        full_list = os.listdir(psm_assignments_path_parent_directory)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return final_list

    # ---------- load inputs ----------
    consolidated_frags = _load_consolidated_frags(sample_output_directory)

    query = ''  # search filter (kept for compatibility with your structure)
    parent_dir_list = get_dir_names_with_strings_list([query])

    for directory in parent_dir_list:
        psm_assignments_path = os.path.join(sample_output_directory, 'final_psm_report_out.csv')
        fragment_match_directory = os.path.join(sample_output_directory, 'fragment_matches')

        psm_assignments = pd.read_csv(psm_assignments_path)
        motif_database = pd.read_csv(motif_database_path)
        motif_sequences = motif_database['Sequence'].values.tolist()

        # Ensure we have required columns
        if 'Sequence with mod' not in psm_assignments.columns:
            # If missing, back-fill with the original sequence (best-effort)
            psm_assignments['Sequence with mod'] = psm_assignments['Sequence']

        psm_assignments = psm_assignments.dropna(subset=['Sequence with mod'])

        motif_match_motif_storage = []
        motif_match_sequence_storage = []
        motif_match_sequence_modded_storage = []
        motif_match_scan_storage = []
        motif_match_coverage_storage = []

        for psm in range(len(psm_assignments)):
            sequence_modded = psm_assignments.iloc[psm]['Sequence with mod']
            sequence = psm_assignments.iloc[psm]['Sequence']
            scan = int(psm_assignments.iloc[psm]['Scan'])
            for motif in motif_sequences:
                if motif in sequence:
                    seq_len = len(sequence)
                    # find all motif spans in the (unmodified/plain) sequence
                    for start, end in [m.span() for m in re.finditer(motif, sequence)]:
                        # theoretical b/y for the motif span
                        b_ion_list = [f"b{val+1}" for val in range(start, end)]
                        y_ion_list = [f"y{seq_len - val}" for val in range(start, end)]
                        theo_ion_list = set(b_ion_list + y_ion_list)

                        # experimental ions: from consolidated or legacy CSV
                        exp_ions = _get_exp_ions(sequence, sequence_modded, scan)
                        exp_b_y = exp_ions['ion'].astype(str).tolist()

                        # strip neutral losses to compare only b/y indices
                        exp_b_y_noPTM = []
                        for c in exp_b_y:
                            c2 = c.replace('-H2O', '').replace('-NH3', '')
                            if c2 in theo_ion_list and c2 not in exp_b_y_noPTM:
                                exp_b_y_noPTM.append(c2)

                        motif_seq_cov_calc = len(exp_b_y_noPTM) / max(len(theo_ion_list), 1)

                        motif_match_motif_storage.append(motif)
                        motif_match_sequence_storage.append(sequence)
                        motif_match_sequence_modded_storage.append(sequence_modded)
                        motif_match_scan_storage.append(scan)
                        motif_match_coverage_storage.append(motif_seq_cov_calc)

        motif_matches = pd.DataFrame({
            'Sequence format': motif_match_sequence_storage,   # plain sequence
            'Sequence': motif_match_sequence_modded_storage,   # modded sequence (file naming)
            'Matched motif': motif_match_motif_storage,
            'Scan': motif_match_scan_storage,
            'Motif Sequence Coverage': motif_match_coverage_storage
        })

        # Fix: motif score uses row-wise motif length
        if not motif_matches.empty:
            motif_matches['Motif Score (1)'] = motif_matches.apply(
                lambda row: len(row['Sequence']) / max(len(row['Matched motif']), 1), axis=1
            )
            motif_matches['Sequence Length'] = motif_matches['Sequence'].str.len()
            motif_matches['Norm. Motif Score (2)'] = motif_matches['Motif Score (1)'] * np.sqrt(motif_matches['Sequence Length'])
            motif_matches['Final motif score'] = motif_matches['Norm. Motif Score (2)'] * (motif_matches['Motif Sequence Coverage'])
        else:
            motif_matches['Motif Score (1)'] = []
            motif_matches['Sequence Length'] = []
            motif_matches['Norm. Motif Score (2)'] = []
            motif_matches['Final motif score'] = []

        # Merge back to PSMs (keys must match plain, modded, and scan)
        psm_w_motif_report = pd.merge(
            psm_assignments,
            motif_matches,
            right_on=['Sequence format', 'Sequence', 'Scan'],
            left_on=['Sequence', 'Sequence with mod', 'Scan'],
            how='inner'
        )

        psm_w_motif_report = psm_w_motif_report.sort_values(by='Final motif score', ascending=False).drop_duplicates(subset='Scan')

        file_path = os.path.join(sample_output_directory, 'psm_results_w_motif.csv')
        with open(file_path, 'w', newline='') as filec:
            writerc = csv.writer(filec)
            psm_w_motif_report.to_csv(filec, index=False)

        return psm_w_motif_report

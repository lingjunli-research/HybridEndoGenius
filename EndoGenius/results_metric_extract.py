# -*- coding: utf-8 -*-
"""
Created on Tue May  9 14:30:59 2023

Schema-tolerant version:
- Handles psm_results_w_motif.csv with columns ['Sequence', 'Sequence with mod']
  OR with merge suffixes ['Sequence_x','Sequence_y'] OR ['Sequence format', 'Sequence'].
- Uses consolidated fragment table if available; falls back to legacy per-scan CSVs.
- Precursor error pulled from precursor_AMM_results.csv.
"""

import os
import re
import csv
from statistics import mean
import pandas as pd

print('Results Metric Extract')

def results_metric_extract_start(results_directory, output_directory, sample_output_directory):

    # ---------- helpers ----------
    def get_dir_names_with_strings_list(str_list):
        full_list = [name for name in os.listdir(results_directory)
                     if os.path.isdir(os.path.join(results_directory, name))]
        return [nm for ps in str_list for nm in full_list if ps in nm]

    def get_file_names_with_strings_list(str_list, folder):
        full_list = os.listdir(folder)
        return [nm for ps in str_list for nm in full_list if ps in nm]

    def count_consec(lst):
        if not lst:
            return [0]
        lst = sorted(lst)
        consec = [1]
        for x, y in zip(lst, lst[1:]):
            if x == y - 1:
                consec[-1] += 1
            else:
                consec.append(1)
        return consec

    def _load_consolidated_frags(base_dir: str):
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

    def _legacy_seq_name(seq_with_mod: str) -> str:
        # Legacy file naming uses "(pyroGlu)" replacement in filenames
        s = seq_with_mod.replace('(Gln->pyro-Glu)', '(pyroGlu)')
        s = s.replace('(Glu->pyro-Glu)', '(pyroGlu)')
        return s

    def _get_fragments(seq_with_mod: str, scan: int) -> pd.DataFrame:
        """
        Return DataFrame with at least ['ion','Fragment error (Da)'] for (Sequence with mod, Scan).
        Prefer consolidated table; else load per-file CSV.
        """
        if consolidated_frags is not None:
            df = consolidated_frags
            out = df[(df.get("Sequence") == seq_with_mod) & (df.get("Scan") == int(scan))]
            need = {"ion", "Fragment error (Da)"}
            if not need.issubset(out.columns):
                # Minimal fallback if schema differs
                cols = [c for c in ["ion", "Fragment error (Da)"] if c in out.columns]
                return out[cols].copy() if cols else pd.DataFrame(columns=["ion", "Fragment error (Da)"])
            return out[["ion", "Fragment error (Da)"]].copy()
        else:
            seq_for_file = _legacy_seq_name(seq_with_mod)
            p = os.path.join(sample_output_directory, "fragment_matches", f"{seq_for_file}_{int(scan)}_fragment_report.csv")
            if not os.path.isfile(p):
                return pd.DataFrame(columns=["ion", "Fragment error (Da)"])
            fr = pd.read_csv(p)
            if "ion" not in fr.columns:
                fr["ion"] = []
            if "Fragment error (Da)" not in fr.columns:
                fr["Fragment error (Da)"] = []
            return fr[["ion", "Fragment error (Da)"]].copy()

    def _pick_col(cols, *candidates):
        for c in candidates:
            if c in cols:
                return c
        return None

    # ---------- load inputs ----------
    consolidated_frags = _load_consolidated_frags(sample_output_directory)

    precursor_amm_path = os.path.join(sample_output_directory, "precursor_AMM_results.csv")
    precursor_amm = pd.read_csv(precursor_amm_path) if os.path.isfile(precursor_amm_path) else pd.DataFrame()
    if not precursor_amm.empty and "Scan" in precursor_amm.columns:
        precursor_amm["Scan"] = precursor_amm["Scan"].astype(int, errors="ignore")

    query = ''
    parent_dir_list = get_dir_names_with_strings_list([query])

    # storages
    peptide_storage = []
    scan_storage = []
    sample_storage = []
    avg_fragment_err_storage = []
    precursor_err_storage = []
    motif_score_storage = []
    consec_b_ions_storage = []
    consec_y_ions_storage = []
    c_term_amid_storage = []
    avg_annotate_peak_storage = []
    percent_b_ions_storage = []
    percent_y_ions_storage = []
    avg_ions_per_AA_storage = []
    num_non_neut_ions_storage = []
    hyperscore_storage = []
    e_pyroglu = []
    q_pyroglu = []
    max_fragment_err_storage = []
    min_fragment_err_storage = []
    med_fragment_err_storage = []
    ox_storage = []
    sulfo_storage = []
    num_mods_storage = []
    seq_cov_storage = []

    for a in parent_dir_list:
        a_path = sample_output_directory
        _ = get_file_names_with_strings_list(['.csv'], a_path)  # retained for parity

        final_results_file_path = os.path.join(sample_output_directory, 'psm_results_w_motif.csv')
        final_results = pd.read_csv(final_results_file_path)

        # ---- choose column names robustly ----
        cols = final_results.columns

        # "plain" (no-mod) sequence column from assignment/motif join
        col_plain = _pick_col(cols, 'Sequence', 'Sequence_x', 'Sequence format')
        # "modded" sequence column
        col_mod = _pick_col(cols, 'Sequence with mod', 'Sequence_y', 'Sequence (modded)', 'Sequence_mod')

        # if only one sequence column exists, assume it's the modded one
        if col_mod is None and col_plain is not None and col_plain not in ('Sequence format',):
            # Heuristic: if there is no explicit mod column but we have 'Sequence',
            # treat it as the modded sequence for downstream fragment lookups.
            col_mod = col_plain

        if col_plain is None:
            # still ensure we have something to read (only used for a couple of metrics)
            col_plain = col_mod  # fall back

        # required numeric columns: tolerate exact names from earlier stages
        col_seqcov = _pick_col(cols, 'Sequence coverage', '% Seq Coverage')
        col_corr = _pick_col(cols, 'Correlation value', 'Hyperscore')
        col_motifscore = _pick_col(cols, 'Final motif score', 'Motif_Score')

        # filters/sanity
        if col_seqcov:
            final_results = final_results[final_results[col_seqcov] > 0]
        if "Scan" in cols:
            final_results["Scan"] = final_results["Scan"].astype(int, errors="ignore")

        for ind in final_results.index:
            sequence_modded = final_results.at[ind, col_mod] if col_mod else ""
            sequence_plain = final_results.at[ind, col_plain] if col_plain else ""
            scan = int(final_results.at[ind, 'Scan']) if 'Scan' in final_results.columns else int(final_results.at[ind, 'Scan_y'])
            motif_score = final_results.at[ind, col_motifscore] if col_motifscore else float('nan')
            corr_score = final_results.at[ind, col_corr] if col_corr else float('nan')
            seq_cov = final_results.at[ind, col_seqcov] if col_seqcov else float('nan')

            peptide_storage.append(sequence_modded)
            scan_storage.append(scan)
            motif_score_storage.append(motif_score)
            sample_storage.append(a)
            hyperscore_storage.append(corr_score)
            seq_cov_storage.append(seq_cov)

            num_mods_storage.append(sequence_modded.count('('))
            c_term_amid_storage.append(1 if '(Amidated)' in sequence_modded else 0)
            ox_storage.append(1 if '(Oxidation)' in sequence_modded else 0)
            sulfo_storage.append(1 if 'Sulfation' in sequence_modded else 0)

            if 'Q(Gln->pyro-Glu)' in sequence_modded:
                q_pyroglu.append(1); e_pyroglu.append(0)
            elif 'E(Glu->pyro-Glu)' in sequence_modded:
                e_pyroglu.append(1); q_pyroglu.append(0)
            else:
                e_pyroglu.append(0); q_pyroglu.append(0)

            # --- fragments ---
            fragment_report = _get_fragments(sequence_modded, scan).copy()
            if not fragment_report.empty and 'Fragment error (Da)' in fragment_report.columns:
                frag_err = fragment_report['Fragment error (Da)'].abs().dropna()
                avg_fragment_err_storage.append(frag_err.mean() if not frag_err.empty else float('nan'))
                max_fragment_err_storage.append(frag_err.max() if not frag_err.empty else float('nan'))
                min_fragment_err_storage.append(frag_err.min() if not frag_err.empty else float('nan'))
                med_fragment_err_storage.append(frag_err.median() if not frag_err.empty else float('nan'))
            else:
                avg_fragment_err_storage.append(float('nan'))
                max_fragment_err_storage.append(float('nan'))
                min_fragment_err_storage.append(float('nan'))
                med_fragment_err_storage.append(float('nan'))

            # precursor error (ppm) from precursor_AMM_results.csv by (Sequence modded, Scan)
            prec_err_ppm = float('nan')
            if not precursor_amm.empty and {'Sequence', 'Scan', 'Precursor error (ppm)'}.issubset(precursor_amm.columns):
                hit = precursor_amm[(precursor_amm['Sequence'] == sequence_modded) & (precursor_amm['Scan'] == scan)]
                if not hit.empty:
                    prec_err_ppm = hit.iloc[0]['Precursor error (ppm)']
            precursor_err_storage.append(prec_err_ppm)

            # ion lists and percentages
            fragment_ion_list = fragment_report['ion'].astype(str).tolist() if 'ion' in fragment_report.columns else []
            num_ions = len(fragment_ion_list)

            b_ions = [ion for ion in fragment_ion_list if ion.startswith('b')]
            y_ions = [ion for ion in fragment_ion_list if ion.startswith('y')]

            percent_b_ions_storage.append((len(b_ions) / num_ions) * 100 if num_ions else 0.0)
            percent_y_ions_storage.append((len(y_ions) / num_ions) * 100 if num_ions else 0.0)

            # Average annotations/fragment (legacy logic)
            peptide_length = len(sequence_modded)
            ion_types = []
            for i in range(2, peptide_length + 1): ion_types.append('b' + str(i))
            for i in range(1, peptide_length): ion_types.append('y' + str(i))
            ion_format_list = []
            for aa in fragment_ion_list:
                ion_format_list.append(aa.replace('-H2O', '').replace('-NH3', ''))
            ion_counts = [ion_format_list.count(ion) for ion in ion_types] if ion_types else []
            avg_annotate_peak_storage.append(mean(ion_counts) if ion_counts else 0.0)

            # positions for consecutive b/y
            def _to_int_list(lst, prefix):
                out = []
                for ion in lst:
                    ion2 = ion.replace(prefix, '').replace('-H2O', '').replace('-NH3', '')
                    if ion2.isdigit(): out.append(int(ion2))
                return out

            b_ion_loc = _to_int_list(b_ions, 'b')
            y_ion_loc = _to_int_list(y_ions, 'y')

            peptide_no_mods = re.sub(r"[\(\[].*?[\)\]]", "", sequence_modded)
            num_expected_ions = len(peptide_no_mods) if peptide_no_mods else 1

            consec_b = count_consec(b_ion_loc); max_consec_b = max(consec_b) if consec_b else 0
            consec_b_ions_storage.append(max_consec_b / num_expected_ions)

            consec_y = count_consec(y_ion_loc); max_consec_y = max(consec_y) if consec_y else 0
            consec_y_ions_storage.append(max_consec_y / num_expected_ions)

            # Average number of fragment ions per AA (b + reoriented y)
            y_reorient = [num_expected_ions - loc for loc in y_ion_loc]
            all_ions = b_ion_loc + y_reorient
            all_ions_no_dups = sorted(set(all_ions))
            counts_per_aa = [all_ions.count(k) for k in all_ions_no_dups] if all_ions_no_dups else []
            avg_ions_per_AA_storage.append(mean(counts_per_aa) if counts_per_aa else 0.0)

            # Number of non-neutral-loss fragment ions per AA
            ion_filtered_loc = []
            for ion in fragment_ion_list:
                if ('-H2O' in ion) or ('-NH3' in ion): continue
                if ion.startswith('b'):
                    val = ion.replace('b', '')
                elif ion.startswith('y'):
                    val = ion.replace('y', '')
                else:
                    continue
                if val.isdigit(): ion_filtered_loc.append(val)
            ion_filtered_loc_no_dups = sorted(set(ion_filtered_loc))
            filtered_counts_per_aa = [ion_filtered_loc.count(k) for k in ion_filtered_loc_no_dups] if ion_filtered_loc_no_dups else []
            num_non_neut_ions_storage.append(mean(filtered_counts_per_aa) if filtered_counts_per_aa else 0.0)

    # assemble output
    final_weighting_metrics = pd.DataFrame({
        'Sample': sample_storage,
        'Peptide': peptide_storage,
        'Scan': scan_storage,
        'Average Fragment Error': avg_fragment_err_storage,
        'Precursor Error': precursor_err_storage,
        '# Consecutive b-ions': consec_b_ions_storage,
        '# Consecutive y-ions': consec_y_ions_storage,
        'C-termini amidation': c_term_amid_storage,
        'Motif_Score': motif_score_storage,
        'Average annotations/fragment ': avg_annotate_peak_storage,  # legacy trailing space kept
        '% Fragment ions are b': percent_b_ions_storage,
        '% Fragment ions are y': percent_y_ions_storage,
        'Average number of fragment ions per AA': avg_ions_per_AA_storage,
        'Number of non-neutral-loss fragment ions per AA': num_non_neut_ions_storage,
        'Hyperscore': hyperscore_storage,
        'Pyro-glu on E': e_pyroglu,
        'Pyro-glu on Q': q_pyroglu,
        'Max Fragment Error': max_fragment_err_storage,
        'Min Fragment Error': min_fragment_err_storage,
        'Median Fragment Error': med_fragment_err_storage,
        'Oxidation': ox_storage,
        'Sulfation': sulfo_storage,
        '# Modifications': num_mods_storage,
        '% Seq Coverage': seq_cov_storage
    })

    file_path = os.path.join(sample_output_directory, 'metrics_extracted_motif.csv')
    with open(file_path, 'w', newline='') as filec:
        writerc = csv.writer(filec)
        final_weighting_metrics.to_csv(filec, index=False)

    return final_weighting_metrics

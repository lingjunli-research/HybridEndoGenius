# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 13:43:47 2025

@author: lafields2
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:41:46 2023

Refactored to support consolidated fragment storage written by database_search.py:
- Uses fragment_matches.parquet or fragment_matches.csv.gz if present
- Falls back to legacy per-file CSVs under fragment_matches/ if needed

@author: lawashburn
"""

import os
import re
import csv
import numpy as np
import pandas as pd

print('PSM assignment')

def PSM_assignment_execute(
    standard_err_percent: float,
    confident_seq_cov: float,
    max_adjacent_swapped_AA: int,
    min_motif_len: int,  
    fragment_error_threshold: float,
    num_sub_AAs: int,
    db_search_parent_directory: str, 
    target_path: str,
    motif_path: str,
    sample_output_directory: str,
):
    # ---------- helpers ----------
    def _load_consolidated_frags(base_dir: str):
        """
        Try to load a single consolidated fragment table produced by database_search.py.
        Looks for:
          - fragment_matches.parquet
          - fragment_matches.csv.gz
        Returns a DataFrame or None if not found.
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

    def _seq_mod_for_legacy(seq_mod: str) -> str:
        # File names use "(pyroGlu)"
        return (seq_mod
                .replace("(Glu->pyro-Glu)", "(pyroGlu)")
                .replace("(Gln->pyro-Glu)", "(pyroGlu)"))

    def _get_fragment_report(seq_plain: str, seq_mod: str, scan: int) -> pd.DataFrame:
        """
        Return fragment rows for (Sequence, Scan).
        - If consolidated table is present, match by 'Sequence' and 'Scan'.
          Try modded string first, then plain as fallback.
        - Else read legacy CSV fragments/<seq_modded>_<scan>_fragment_report.csv
        """
        if consolidated_frags is not None:
            df = consolidated_frags
            out = df[(df.get("Sequence") == seq_mod) & (df.get("Scan") == int(scan))]
            if out.empty:
                out = df[(df.get("Sequence") == seq_plain) & (df.get("Scan") == int(scan))]
            return out.copy()
        path = os.path.join(sample_output_directory, "fragment_matches",
                            f"{_seq_mod_for_legacy(seq_mod)}_{int(scan)}_fragment_report.csv")
        if not os.path.isfile(path):
            return pd.DataFrame(columns=["ion", "Fragment error (Da)"])
        fr = pd.read_csv(path)
        # ensure required columns exist
        if "ion" not in fr.columns:
            fr["ion"] = []
        if "Fragment error (Da)" not in fr.columns:
            fr["Fragment error (Da)"] = []
        return fr[["ion", "Fragment error (Da)"]].copy()

    def compare(string1, string2, no_match_c=' ', match_c='|'):
        if len(string2) < len(string1):
            string1, string2 = string2, string1
        result = ''
        n_diff = 0
        for c1, c2 in zip(string1, string2):
            if c1 == c2:
                result += match_c
            else:
                result += no_match_c
                n_diff += 1
        delta = len(string2) - len(string1)
        result += delta * no_match_c
        n_diff += delta
        return (result, n_diff)

    def dif_compare(string1, string2):
        result, n_diff = compare(string1, string2, no_match_c='_')
        return n_diff

    def get_dir_names_with_strings_list(str_list):  # definition for finding a file containing a string in filename in specified directory
        full_list = os.listdir(db_search_parent_directory)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]

        final_final_list = []
        for a in final_list:
            isdir = os.path.isdir(db_search_parent_directory + '\\' + a)
            if isdir:
                final_final_list.append(a)
        return final_final_list
    
    def _min_frag_error(seq_plain: str, seq_with_mod: str, scan: int) -> float:
        fr = _get_fragment_report(seq_plain, seq_with_mod, scan)
        if fr.empty or "Fragment error (Da)" not in fr.columns:
            return float("inf")
        s = pd.to_numeric(fr["Fragment error (Da)"], errors="coerce").abs()
        return float(np.nanmin(s)) if len(s) else float("inf")
    
    # ---------- small utilities ----------
    def _strip_mods(seq: str) -> str:
        return re.sub(r"[\(\[].*?[\)\]]", "", seq)

    def _longest_adjacent_diff_block_len(s1: str, s2: str) -> int:
        # length of the longest run of consecutive indices where s1 != s2
        diffs = [i for i in range(min(len(s1), len(s2))) if s1[i] != s2[i]]
        if not diffs:
            return 0
        # group consecutive integers
        best, cur, prev = 1, 1, diffs[0]
        for x in diffs[1:]:
            if x == prev + 1:
                cur += 1
            else:
                best = max(best, cur)
                cur = 1
            prev = x
        return max(best, cur)

    def _differences_within(s1: str, s2: str, max_diffs: int) -> bool:
        # Hamming-like (same length assumed by callers), stop if exceeds
        diffs = 0
        for a, b in zip(s1, s2):
            if a != b:
                diffs += 1
                if diffs > max_diffs:
                    return False
        return True

    def _motif_union_coverage(seq_plain: str, motifs: list[str]) -> float:
        """
        Fraction of sequence positions covered by any matched motif (union).
        """
        positions = set()
        for m in motifs:
            if not m:
                continue
            for match in re.finditer(re.escape(m), seq_plain):
                positions.update(range(match.start(), match.end()))
        return (len(positions) / max(len(seq_plain), 1)) if positions else 0.0

    def _choose_extreme(df: pd.DataFrame, col: str, how: str) -> pd.DataFrame:
        if df.empty or col not in df.columns:
            return df
        if how == "max":
            val = df[col].max()
        else:
            val = df[col].min()
        return df[df[col] == val]

    def _choose_correlation_rounded(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "Correlation value" not in df.columns:
            return df
        tmp = df.copy()
        tmp["__corr_round__"] = tmp["Correlation value"].round(0)
        mx = tmp["__corr_round__"].max()
        out = tmp[tmp["__corr_round__"] == mx].drop(columns="__corr_round__")
        return out

    def _tie_same_length_and_swaps(df: pd.DataFrame) -> pd.DataFrame:
        """
        If all sequences same length and two candidates, pick the pair whose
        longest adjacent swap block <= threshold. If satisfied, return that pair;
        if exactly two, we then accept both (the original logic later samples one).
        """
        if df.empty:
            return df
        tmp = df.copy()
        tmp["Seq_Len"] = tmp["Sequence"].str.len()
        if tmp["Seq_Len"].nunique() != 1:
            return df
        seqs = tmp["Sequence"].tolist()
        if len(seqs) != 2:
            return df
        s1, s2 = seqs
        if _longest_adjacent_diff_block_len(s1, s2) <= max_adjacent_swapped_AA:
            # keep both; caller decides next step (often labels H(+))
            return df
        # else keep original df; caller will sample/random later
        return df

    def _count_frags_ge_thresh(seq_plain: str, seq_mod: str, scan: int) -> int:
        fr = _get_fragment_report(seq_plain, seq_mod, scan)
        if fr.empty or "Fragment error (Da)" not in fr.columns:
            return 0
        return int((fr["Fragment error (Da)"] >= fragment_error_threshold).sum())

    def _avg_frag_error_le(seq_plain: str, seq_mod: str, scan: int):
        fr = _get_fragment_report(seq_plain, seq_mod, scan)
        if fr.empty or "Fragment error (Da)" not in fr.columns:
            return np.inf
        keep = fr[fr["Fragment error (Da)"] <= fragment_error_threshold]["Fragment error (Da)"].abs()
        return keep.mean() if not keep.empty else np.inf

    # ---------- load inputs ----------
    target = pd.read_csv(target_path)
    target_list = target['Sequence'].values.tolist()

    motif_db = pd.read_csv(motif_path)
    motif_list = motif_db['Sequence'].values.tolist()

    consolidated_frags = _load_consolidated_frags(sample_output_directory)

    # ---------- main ----------
    query = ''  # search for ion list pertaining to the sequence
    parent_dir_list = get_dir_names_with_strings_list([query])  # not really used (kept for structure compatibility)

    all_corr_path = os.path.join(sample_output_directory, "all_correlation_results.csv")
    all_corr = pd.read_csv(all_corr_path)
    all_corr["Sequence with mod"] = all_corr["Sequence"].astype(str)
    all_corr["Sequence"] = all_corr["Sequence"].astype(str).str.replace(r"\([^()]*\)", "", regex=True)

    max_seq_cov = all_corr["Sequence coverage"].max()
    _ = all_corr[all_corr["Sequence coverage"] >= max_seq_cov * (1 - standard_err_percent)]["Correlation value"].std()

    # mark singleton scans
    all_corr["count"] = all_corr.groupby("Scan")["Scan"].transform("size")

    # collect final rows
    assigned_rows = []

        # ---------- precompute motif flags & stats per Sequence (plain) ----------
    motif_map = {}
    for seq in all_corr["Sequence"].unique():
        has_motif = False
        first_motif = np.nan
        matched = []
        for m in motif_list:
            if m and m in seq:
                matched.append(m)
                if not has_motif:
                    first_motif = m
                    has_motif = True
        # metrics used in tie-breakers
        motif_len = len(first_motif) if isinstance(first_motif, str) else np.nan
        motif_seq_ratio = (motif_len * np.sqrt(len(seq))) if isinstance(motif_len, (int, float)) else np.nan
        motif_count = len(matched)
        motif_cov = _motif_union_coverage(seq, matched) if matched else 0.0
        motif_map[seq] = {
            "Motif": first_motif,
            "Motif status": bool(matched),
            "Motif Length": motif_len if matched else np.nan,
            "Motif:Seq Ratio": motif_seq_ratio if matched else np.nan,
            "# Matching motifs": motif_count,
            "Sequence/Motif Coverage": motif_cov,
        }

    motif_df = pd.DataFrame.from_dict(motif_map, orient="index").reset_index().rename(columns={"index": "Sequence"})
    work = all_corr.merge(motif_df, on="Sequence", how="left")
    work["Seq Length"] = work["Sequence"].str.len()
    
    # ---------- per-scan resolver ----------
    def resolve_scan(scan_df: pd.DataFrame) -> pd.DataFrame:
        """
        Return exactly one row with 'Step assigned' filled, using the original
        A/B/C… cascade but without repeated code.
        """
        # Helper to annotate step and stop when single row remains
        def _apply_and_mark(df_in, df_out, step_code):
            if df_out.empty:
                return df_in, False  # no change
            if len(df_out) == 1:
                out = df_out.copy()
                out["Step assigned"] = step_code
                return out, True
            return df_out.copy(), False

        # A: motif present?
        with_motif = scan_df[scan_df["Motif status"] == True]
        if len(with_motif) == 1:
            out = with_motif.copy()
            out["Step assigned"] = "filter A(+)"
            return out

        if len(with_motif) > 1:
            # B: max Motif:Seq Ratio
            tmp = _choose_extreme(with_motif, "Motif:Seq Ratio", "max")
            tmp, done = _apply_and_mark(with_motif, tmp, "filter A(++) B(+)")
            if done:
                return tmp

            # C: max # Matching motifs
            tmp2 = _choose_extreme(tmp, "# Matching motifs", "max")
            tmp2, done = _apply_and_mark(tmp, tmp2, "filter A(++) B(++) C(+)")
            if done:
                return tmp2

            # D: max Sequence/Motif Coverage
            tmp3 = _choose_extreme(tmp2, "Sequence/Motif Coverage", "max")
            tmp3, done = _apply_and_mark(tmp2, tmp3, "filter A(++) B(++) C(++) D(+)")
            if done:
                return tmp3

            # E: require Sequence coverage >= confident_seq_cov
            tmp4 = tmp3[tmp3["Sequence coverage"] >= confident_seq_cov]
            if len(tmp4) == 1:
                out = tmp4.copy(); out["Step assigned"] = "filter A(++) B(++) C(++) D(++) E(+)"
                return out
            if len(tmp4) > 1:
                tmp3 = tmp4  # continue cascade on this subset

            # G: max Correlation value
            tmp5 = _choose_extreme(tmp3, "Correlation value", "max")
            tmp5, done = _apply_and_mark(tmp3, tmp5, "filter A(++) B(++) C(++) D(++) E(++) G(+)")
            if done:
                return tmp5

            # I: max Sequence coverage
            tmp6 = _choose_extreme(tmp5, "Sequence coverage", "max")
            tmp6, done = _apply_and_mark(tmp5, tmp6, "filter A(++) B(++) C(++) D(++) E(++) G(++) I(+)")
            if done:
                return tmp6

            # H: same-length & swapped-AA block <= threshold?
            tmp7 = _tie_same_length_and_swaps(tmp6)
            if len(tmp7) == 2 and tmp7["Seq Length"].nunique() == 1:
                # both acceptable; mark H(+) and keep both—sample one for determinism
                out = tmp7.sample(n=1, random_state=1).copy()
                out["Step assigned"] = "filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(+)"
                return out
            # J: tiebreak by number of frags with error >= threshold
            if len(tmp6) > 1:
                tmp6 = tmp6.copy()
                tmp6["__frag_len__"] = [
                    _count_frags_ge_thresh(r["Sequence"], r["Sequence with mod"], r["Scan"])
                    for _, r in tmp6.iterrows()
                ]
                tmp7 = _choose_extreme(tmp6, "__frag_len__", "max").drop(columns="__frag_len__", errors="ignore")
                if len(tmp7) == 1:
                    out = tmp7.copy(); out["Step assigned"] = "filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(-) J(+)"
                    return out

            # fallback: sample 1 (keeps behavior consistent with old code)
            out = tmp6.sample(n=1, random_state=1).copy()
            out["Step assigned"] = "filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(-?)"
            return out

        # ---- No motif path ----
        # E: Sequence coverage >= confident_seq_cov
        filter_E = scan_df[scan_df["Sequence coverage"] >= confident_seq_cov]
        if len(filter_E) == 1:
            out = filter_E.copy(); out["Step assigned"] = "filter A(-) E(+)"
            return out
        if len(filter_E) > 1:
            # I: max Sequence coverage
            tmp = _choose_extreme(filter_E, "Sequence coverage", "max")
            if len(tmp) == 1:
                out = tmp.copy(); out["Step assigned"] = "filter A(-) E(+) I(+)"
                return out
            # K: correlation rounded
            tmp2 = _choose_correlation_rounded(tmp)
            if len(tmp2) == 1:
                out = tmp2.copy(); out["Step assigned"] = "filter A(-) E(++) I(++) K(+)"
                return out
            # H: same-length & swaps allowed
            tmp3 = _tie_same_length_and_swaps(tmp2)
            if len(tmp3) == 2 and tmp3["Seq Length"].nunique() == 1:
                out = tmp3.sample(n=1, random_state=1).copy()
                out["Step assigned"] = "filter A(-)E(++)I(++)K(++)H(+)"
                return out
            out = tmp2.sample(n=1, random_state=1).copy()
            out["Step assigned"] = "filter A(-)E(++)I(++)K(++)H(-?)"
            return out

        # I: max Sequence coverage among all
        tmp = _choose_extreme(scan_df, "Sequence coverage", "max")
        if len(tmp) == 1:
            out = tmp.copy(); out["Step assigned"] = "filter A(-) E(-) I(+)"
            return out

        # F: keep candidates with Correlation value > 0
        filter_F = tmp[tmp["Correlation value"] > 0] if "Correlation value" in tmp.columns else tmp
        if len(filter_F) == 1:
            out = filter_F.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(+)"
            return out
        if len(filter_F) > 1:
            # J: max count of frags with error >= threshold
            tmpJ = filter_F.copy()
            tmpJ["__frag_len__"] = [
                _count_frags_ge_thresh(r["Sequence"], r["Sequence with mod"], r["Scan"])
                for _, r in tmpJ.iterrows()
            ]
            bestJ = _choose_extreme(tmpJ, "__frag_len__", "max").drop(columns="__frag_len__", errors="ignore")
            if len(bestJ) == 1:
                out = bestJ.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(++) J(+)"
                return out
            # G: max correlation
            bestG = _choose_extreme(bestJ, "Correlation value", "max")
            if len(bestG) == 1:
                out = bestG.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(++) J(++) G(+)"
                return out
            # H: same-length & swaps
            h = _tie_same_length_and_swaps(bestG)
            if len(h) == 2 and h["Seq Length"].nunique() == 1:
                out = h.sample(n=1, random_state=1).copy()
                out["Step assigned"] = "filter A(-) E(-) I(++) F(++) J(++) G(++) H(+)"
                return out
            out = bestG.sample(n=1, random_state=1).copy()
            out["Step assigned"] = "filter A(-) E(-) I(++) F(++) J(++) G(++) H(-?)"
            return out

        # L: choose min avg fragment error (<= threshold), else O: min fragment error
        cand = tmp
        cand = cand.copy()
        cand["__avg_err__"] = [
            _avg_frag_error_le(r["Sequence"], r["Sequence with mod"], r["Scan"]) for _, r in cand.iterrows()
        ]
        if np.isfinite(cand["__avg_err__"]).any():
            bestL = _choose_extreme(cand[np.isfinite(cand["__avg_err__"])], "__avg_err__", "min").drop(columns="__avg_err__", errors="ignore")
            if len(bestL) == 1:
                out = bestL.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(+)"
                return out
            # N/H fallbacks
            bestL = bestL.drop_duplicates(subset="Sequence")
            if len(bestL) == 1:
                out = bestL.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(++) (N+)"
                return out
            # H on bestL (same-length & swaps)
            h2 = _tie_same_length_and_swaps(bestL)
            if len(h2) == 2 and h2["Seq Length"].nunique() == 1:
                out = h2.sample(n=1, random_state=1).copy()
                out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(++) (N++) H(+)"
                return out
            # O: min fragment error among all candidates
        # O branch: compute minimal single-frag error
        tmpO = scan_df.copy()
        tmpO["__min_err__"] = tmpO.apply(lambda r: _min_frag_error(r["Sequence"], r["Sequence with mod"], int(r["Scan"])),axis=1)
        for _, r in tmpO.iterrows():
            fr = _get_fragment_report(r["Sequence"], r["Sequence with mod"], r["Scan"])
            tmpO.at[_, "__min_err__"] = fr["Fragment error (Da)"].min() if ("Fragment error (Da)" in fr.columns and not fr.empty) else np.inf
        bestO = _choose_extreme(tmpO, "__min_err__", "min").drop(columns="__min_err__", errors="ignore")
        if len(bestO) == 1:
            out = bestO.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(-) O(+)"
            return out
        # N/P: dedupe, then if 2 left and small diffs, accept; else sample
        bestN = bestO.drop_duplicates(subset="Sequence")
        if len(bestN) == 1:
            out = bestN.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(-) O(++) N(+)"
            return out
        if len(bestN) == 2:
            s1, s2 = bestN["Sequence"].tolist()
            if _differences_within(s1, s2, num_sub_AAs):
                out = bestN.copy().sample(n=1, random_state=1)
                out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(-) O(++) N(+) P(+)"
                return out
        # G: final resort max correlation
        bestG2 = _choose_extreme(bestN, "Correlation value", "max")
        if len(bestG2) == 1:
            out = bestG2.copy(); out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(-) O(++) N(-) G(+)"
            return out
        out = bestN.sample(n=1, random_state=1).copy()
        out["Step assigned"] = "filter A(-) E(-) I(++) F(-) L(-) O(++) N(-) G(-)"
        return out

    # ---------- run per-scan ----------
    singles = work[work["count"] == 1].copy()
    singles["Step assigned"] = "Only peptide match for scan"

    multis = work[work["count"] > 1].copy()
    chosen = []
    for scan_val, scan_df in multis.groupby("Scan", sort=False):
        chosen.append(resolve_scan(scan_df))

    final_psm = pd.concat([singles] + chosen, ignore_index=True) if chosen else singles

    # ---------- target/decoy and output ----------
    final_psm_target = final_psm[final_psm["Sequence"].isin(target_list)].copy()
    final_psm_decoy = final_psm[~final_psm["Sequence"].isin(target_list)].copy()
    final_psm_target["Status"] = "Target"
    final_psm_decoy["Status"] = "Decoy"
    final_psm_output = pd.concat([final_psm_target, final_psm_decoy], ignore_index=True)

    out_csv = os.path.join(sample_output_directory, "final_psm_report_out.csv")
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        final_psm_output.to_csv(fh, index=False)

    return final_psm_output

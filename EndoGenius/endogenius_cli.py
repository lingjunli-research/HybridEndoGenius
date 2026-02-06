# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 17:17:40 2025

@author: lafields2
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
import argparse
import os
import sys
from typing import Callable, Iterable, Tuple, Dict, List, Optional
import pandas as pd
import json
from copy import deepcopy

# --- local imports (unchanged dependencies in your project) ---
from utils import fasta_to_df  # you used this in make_target_list

def _is_missing(x) -> bool:
    """True if x is None/empty string/empty DataFrame/Series/container."""
    if x is None:
        return True
    if isinstance(x, str):
        return len(x.strip()) == 0
    try:
        import pandas as _pd
        if isinstance(x, (_pd.DataFrame, _pd.Series)):
            return x.empty
    except Exception:
        pass
    if isinstance(x, (list, tuple, dict, set)):
        return len(x) == 0
    return False

def _debug_shape(name: str, obj):
    """Compact debug print for stage outputs."""
    import pandas as pd
    if isinstance(obj, pd.DataFrame):
        print(f"[DEBUG] {name}: DataFrame shape={obj.shape}")
    elif isinstance(obj, pd.Series):
        print(f"[DEBUG] {name}: Series len={len(obj)}")
    elif isinstance(obj, (list, tuple, dict, set)):
        print(f"[DEBUG] {name}: {type(obj).__name__} size={len(obj)}")
    else:
        print(f"[DEBUG] {name}: {type(obj).__name__} -> {obj!r}")

# ----------------------------
# Small utilities (ported)
# ----------------------------

def save_csv(df: pd.DataFrame, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

def sibling_mzml(path_like: str) -> Optional[str]:
    if not isinstance(path_like, str):
        return None
    # do replacement case-insensitively on the tail
    p = Path(path_like)
    stem = p.name
    low = stem.lower()
    repls = [("_formatted.ms2", ".mzML"), ("_formatted.txt", ".mzML"),
             (".txt", ".mzML"), (".ms2", ".mzML")]
    for old, new in repls:
        if low.endswith(old):
            candidate = str(p.with_name(stem[:len(stem)-len(old)] + new))
            return candidate if os.path.isfile(candidate) else None
    # last resort: same stem with .mzML
    candidate = str(p.with_suffix(".mzML"))
    return candidate if os.path.isfile(candidate) else None

def require_truthy(items: Iterable[Tuple[str, str, str]]) -> Optional[str]:
    for _, val, msg in items:
        if not str(val).strip():
            return msg
    return None

# ----------------------------
# Mod configuration
# ----------------------------
# We mirror your TK IntVars as boolean CLI flags. The make_db() path uses the
# same variable_mod_dict schema your db_generator.make_a_DB expects.

def build_variable_mod_dict(args: argparse.Namespace) -> Dict[str, object]:
    d: Dict[str, object] = {}

    def set_true(key: str) -> None:
        d[key] = True

    def append_list(key: str, value: str) -> None:
        d.setdefault(key, []).append(value)

    # Terminal / fixed-site flags
    if args.amid:                 set_true('-Amidated')
    if args.acet_nterm:           set_true('Acetylation-')
    if args.carb_nterm:           set_true('Carbamidomethylation-')
    if args.plex12_nterm:         set_true('12PlexDiLeu-')
    if args.mddileu_1101_nterm:   set_true('mdDiLeu1101-')
    if args.mddileu_0400_nterm:   set_true('mdDiLeu0400-')
    if args.pg_e:                 d['(Glu->pyro-Glu)'] = ['ntermE']
    if args.pg_q:                 d['(Gln->pyro-Glu)'] = ['ntermQ']

    # Residue-scoped variable mods
    if args.acet_k:               append_list('(Acetylation)', 'K')
    if args.acet_s:               append_list('(Acetylation)', 'S')
    if args.acet_t:               append_list('(Acetylation)', 'T')

    if args.biotin_k:             append_list('(Biotinylation)', 'K')

    if args.carb_c:               append_list('(Carbamidomethylation)', 'C')
    if args.carb_d:               append_list('(Carbamidomethylation)', 'D')
    if args.carb_e:               append_list('(Carbamidomethylation)', 'E')
    if args.carb_h:               append_list('(Carbamidomethylation)', 'H')
    if args.carb_k:               append_list('(Carbamidomethylation)', 'K')

    if args.carboxy_d:            append_list('(Carboxylation)', 'D')
    if args.carboxy_e:            append_list('(Carboxylation)', 'E')
    if args.carboxy_k:            append_list('(Carboxylation)', 'K')
    if args.carboxy_w:            append_list('(Carboxylation)', 'W')

    if args.deamid_n:             append_list('(Deamidation)', 'N')
    if args.deamid_q:             append_list('(Deamidation)', 'Q')

    if args.dehyd_d:              append_list('(Dehydration)', 'D')
    if args.dehyd_s:              append_list('(Dehydration)', 'S')
    if args.dehyd_t:              append_list('(Dehydration)', 'T')
    if args.dehyd_y:              append_list('(Dehydration)', 'Y')

    if args.methyl_k:             append_list('(Methylation)', 'K')
    if args.methyl_r:             append_list('(Methylation)', 'R')

    if args.ox_m:                 append_list('(Oxidation)', 'M')

    if args.sodi_d:               append_list('(SodiumAdduct)', 'D')
    if args.sodi_e:               append_list('(SodiumAdduct)', 'E')

    if args.plex12_k:             append_list('(12PlexDiLeu)', 'K')
    if args.mddileu_1101_k:       append_list('(mdDiLeu1101)', 'K')
    if args.mddileu_0400_k:       append_list('(mdDiLeu0400)', 'K')

    if args.sulfo_y:              append_list('(Sulfation)', 'Y')
    if args.phospho_s:            append_list('(Phosphorylation)', 'S')
    if args.phospho_t:            append_list('(Phosphorylation)', 'T')
    if args.phospho_y:            append_list('(Phosphorylation)', 'Y')

    return d

# ----------------------------
# Core actions (ported)
# ----------------------------

def make_target_list(fasta_path: str, output_folder: str) -> str:
    target_fasta = fasta_to_df(fasta_path)
    target_list = pd.DataFrame({'Sequence': target_fasta})
    file_path = os.path.join(output_folder, 'target_list.csv')
    save_csv(target_list, file_path)
    return file_path

def make_db(args: argparse.Namespace) -> str:
    from db_generator import make_a_DB
    variable_mod_dict = build_variable_mod_dict(args)
    if getattr(args, "verbose", False):
        print("Variable mod dict:", variable_mod_dict)
    max_mods_number = int(args.max_mods)
    out = make_a_DB(variable_mod_dict, args.fasta, args.out_dir, max_mods_number)
    if not out:
        raise RuntimeError("make_a_DB() returned no path (None/empty).")
    return out

def checked_clear_begin_search(
    predefined_db_path: str,
    output_parent_directory: str,
    raw_file_formatted_path: str,
    target_path: str,
    args: argparse.Namespace,
) -> None:
    """
    Run the EndoGenius first-pass search pipeline (CLI).
    Mirrors GUI order; uses args instead of Tk variables.
    """
    # Lazy imports to match project layout
    from database_search import raw_file_detail_extraction, launch_db_search_pt1
    from PSM_assignment import PSM_assignment_execute
    from motif_search import start_motif_search
    from results_metric_extract import results_metric_extract_start
    from metric_handling import metric_handling_apply
    from target_decoy_assess import target_decoy_apply
    from EndoGenius_Score_Apply import endogenius_apply

    # Resolve per-sample output dir
    details = raw_file_detail_extraction(raw_file_formatted_path, output_parent_directory)
    sample_name, sample_output_directory = details[0], details[1]

    # Numeric knobs
    precursor_error_cutoff = float(args.precursor_err)
    fragment_error_cutoff  = float(args.fragment_err)
    min_mz                 = float(args.mz_min)
    min_intensity_pass     = int(args.min_intensity)
    standard_err_percent   = float(args.standard_err)

    # Mod flags (same order as GUI)
    mod_flags = [
        args.amid, args.acet_nterm, args.acet_k, args.acet_s, args.acet_t,
        args.biotin_k, args.carb_nterm, args.carb_c, args.carb_d, args.carb_e,
        args.carb_h, args.carb_k, args.carboxy_d, args.carboxy_e, args.carboxy_k,
        args.carboxy_w, args.deamid_n, args.deamid_q, args.dehyd_d, args.dehyd_s,
        args.dehyd_t, args.dehyd_y, args.methyl_k, args.methyl_r, args.ox_m,
        args.sodi_d, args.sodi_e, args.plex12_nterm, args.plex12_k,
        args.mddileu_1101_nterm, args.mddileu_1101_k, args.mddileu_0400_nterm, args.mddileu_0400_k,
        args.pg_e, args.pg_q, args.sulfo_y, args.phospho_s, args.phospho_t, args.phospho_y
    ]
    mod_values = [int(b) for b in mod_flags]

    max_modifications  = int(args.max_mods)
    choose_mzml_source = args.fmt_ms2 if args.fmt_ms2 else args.raw_ms2

    if getattr(args, "verbose", False):
        print("Calling launch_db_search_pt1 with:")
        print(f"  predefined_db_path      = {predefined_db_path}")
        print(f"  output_parent_directory = {output_parent_directory}")
        print(f"  choose_mzml_source      = {choose_mzml_source}")
        print(f"  raw_file_formatted_path = {raw_file_formatted_path}")
        print(f"  precursor_error_cutoff  = {precursor_error_cutoff}")
        print(f"  fragment_error_cutoff   = {fragment_error_cutoff}")
        print(f"  min_mz                  = {min_mz}")
        print(f"  min_intensity_pass      = {min_intensity_pass}")
        print(f"  standard_err_percent    = {standard_err_percent}")
        print(f"  max_modifications       = {max_modifications}")
        print(f"  sample_output_directory = {sample_output_directory}")

    # ---- DB search (pt1)
    first_pass_db = launch_db_search_pt1(
        predefined_db_path,
        output_parent_directory,
        choose_mzml_source,
        raw_file_formatted_path,
        precursor_error_cutoff,
        fragment_error_cutoff,
        min_mz,
        min_intensity_pass,
        standard_err_percent,
        *mod_values,
        max_modifications,
        sample_output_directory
    )
    if getattr(args, "verbose", False):
        _debug_shape("launch_db_search_pt1", first_pass_db)
    if _is_missing(first_pass_db):
        raise RuntimeError("launch_db_search_pt1() produced no results.")
    print('First pass DB complete')

    # ---- PSM assignment
    confident_seq_cov       = float(args.cov_thr)
    max_adjacent_swapped_AA = int(args.max_adjacent_swapped_AAs)
    min_motif_len           = int(args.min_motif_len)
    num_sub_AAs             = int(args.max_swapped_AA)
    motif_path              = args.motif_db

    first_pass_PSM_assign = PSM_assignment_execute(
        standard_err_percent, confident_seq_cov, max_adjacent_swapped_AA, min_motif_len,
        fragment_error_cutoff, num_sub_AAs, output_parent_directory, target_path, motif_path, sample_output_directory)
    if getattr(args, "verbose", False):
        _debug_shape("PSM_assignment_execute", first_pass_PSM_assign)
    print('First pass PSM assign complete')

    # ---- Motif search
    first_pass_motif_search = start_motif_search(output_parent_directory, motif_path, sample_output_directory)
    if getattr(args, "verbose", False):
        _debug_shape("start_motif_search", first_pass_motif_search)
    print('First pass motif search complete')

    # ---- Metric extract
    first_pass_weighting_extract = results_metric_extract_start(
        output_parent_directory, output_parent_directory, sample_output_directory)
    if getattr(args, "verbose", False):
        _debug_shape("results_metric_extract_start", first_pass_weighting_extract)
    if _is_missing(first_pass_weighting_extract):
        raise RuntimeError("results_metric_extract_start() produced no results.")
    print('First pass weighting extract complete')

    # ---- Metric apply
    first_pass_metric_apply = metric_handling_apply(
        first_pass_weighting_extract, output_parent_directory, sample_output_directory)
    if getattr(args, "verbose", False):
        _debug_shape("metric_handling_apply", first_pass_metric_apply)
    if _is_missing(first_pass_metric_apply):
        raise RuntimeError("metric_handling_apply() produced no results.")
    print('First pass metric apply complete')

    # ---- Final filtering: FDR or EndoGenius
    if args.fdr is not None and args.eg is not None:
        raise ValueError('Must input either an FDR cutoff OR EndoGenius Score cutoff, not both.')

    if args.fdr is not None:
        target_decoy_apply(first_pass_metric_apply, target_path, output_parent_directory,
                           float(args.fdr), sample_output_directory, raw_file_formatted_path)
    elif args.eg is not None:
        endogenius_apply(first_pass_metric_apply, target_path, output_parent_directory,
                         float(args.eg), sample_output_directory, raw_file_formatted_path)
    else:
        raise ValueError('Must input either an FDR cutoff OR EndoGenius Score cutoff.')
    print('First pass target-decoy complete')

def begin_search_confirmed(args: argparse.Namespace) -> None:
    """
    Decide branches (prebuilt DB vs build; raw vs formatted) and
    call checked_clear_begin_search(...).
    """
    out_dir = args.out_dir

    # If using a prebuilt DB, a target list is required
    if args.db_csv:
        if not args.target_list:
            raise ValueError("When using --db-csv, you must supply --target-list.")
        predefined_db_path = args.db_csv
        target_path = args.target_list

        # Prefer formatted, otherwise format the raw
        if args.fmt_ms2:
            raw_file_formatted_path = args.fmt_ms2
        else:
            from format_MS2_file_RT_IIT import format_raw_MS2
            raw_file_formatted_path = format_raw_MS2(args.raw_ms2, out_dir)
            if _is_missing(raw_file_formatted_path):
                raise RuntimeError("format_raw_MS2() returned no path.")
        checked_clear_begin_search(predefined_db_path, out_dir, raw_file_formatted_path, target_path, args)
        print("Your database search is complete.")
        return

    # Build DB & target list from FASTA
    from utils import fasta_to_df
    from db_generator import make_a_DB

    variable_mod_dict = build_variable_mod_dict(args)  # you already have this helper in the CLI
    max_mods_number   = int(args.max_mods)
    new_db_path       = make_a_DB(variable_mod_dict, args.fasta, out_dir, max_mods_number)

    # target list
    target_fasta = fasta_to_df(args.fasta)
    import pandas as _pd, os as _os
    target_path  = _os.path.join(out_dir, 'target_list.csv')
    _pd.DataFrame({'Sequence': target_fasta}).to_csv(target_path, index=False)

    # format or use formatted
    if args.fmt_ms2:
        raw_file_formatted_path = args.fmt_ms2
    else:
        from format_MS2_file_RT_IIT import format_raw_MS2
        raw_file_formatted_path = format_raw_MS2(args.raw_ms2, out_dir)
        if _is_missing(raw_file_formatted_path):
            raise RuntimeError("format_raw_MS2() returned no path.")

    checked_clear_begin_search(new_db_path, out_dir, raw_file_formatted_path, target_path, args)
    print("Your database search is complete.")


def _print_args_summary(args):
    # Only the most relevant knobs + any that are present
    def _get(name): return getattr(args, name, None)
    fields = {
        "name": _get("name") or "<unnamed>",
        "raw_ms2": _get("raw_ms2"),
        "fmt_ms2": _get("fmt_ms2"),
        "fasta": _get("fasta"),
        "db_csv": _get("db_csv"),
        "target_list": _get("target_list"),
        "motif_db": _get("motif_db"),
        "out_dir": _get("out_dir"),
        "fdr": _get("fdr"),
        "eg": _get("eg"),
        "mz_min": _get("mz_min"),
        "mz_max": _get("mz_max"),
        "precursor_err": _get("precursor_err"),
        "fragment_err": _get("fragment_err"),
        "max_mods": _get("max_mods"),
    }
    print("Resolved run arguments:")
    for k, v in fields.items():
        print(f"  - {k}: {v}")

def _must_exist(path: str, label: str):
    if not isinstance(path, str) or not path.strip():
        raise ValueError(f"{label} is missing or empty.")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")

def _ensure_dir(path: str, label: str):
    if not isinstance(path, str) or not path.strip():
        raise ValueError(f"{label} is missing or empty.")
    Path(path).mkdir(parents=True, exist_ok=True)

def _preflight(args: argparse.Namespace) -> None:
    """Strict validation + existence checks before running anything heavy."""
    if getattr(args, "verbose", False):
        _print_args_summary(args)

    # Required text presence (as you had)
    err = require_truthy([
        ("mz_min", args.mz_min,            'Input minimum m/z value'),
        ("mz_max", args.mz_max,            'Input maximum m/z value'),
        ("min_int", args.min_intensity,    'Input minimum intensity value'),
        ("max_prec_z", args.max_precursor_z,'Input maximum precursor charge value'),
        ("max_frag_z", args.max_fragment_z, 'Input maximum fragment charge value'),
        ("prec_err", args.precursor_err,   'Input maximum precursor error value'),
        ("frag_err", args.fragment_err,    'Input maximum fragment error value'),
        ("max_mods", args.max_mods,        'Input maximum # modifications per peptide (0 if none)'),
        ("motif_db", args.motif_db,        'Select motif database'),
        ("cov_thr", args.cov_thr,          'Input confident coverage threshold'),
        ("out_dir", args.out_dir,          'Select output folder'),
    ])
    if err:
        raise ValueError(err)

    # XORs
    if not bool(args.raw_ms2) ^ bool(args.fmt_ms2):
        raise ValueError('Provide either --raw-ms2 OR --fmt-ms2 (exactly one).')
    if not bool(args.fasta) ^ bool(args.db_csv):
        raise ValueError('Provide either --fasta OR --db-csv (exactly one).')
    if (args.fdr is None) == (args.eg is None):
        raise ValueError('Provide either --fdr OR --eg (exactly one).')
    if args.db_csv and not args.target_list:
        raise ValueError('When using --db-csv, you must also supply --target-list.')

    # Existence checks
    _ensure_dir(args.out_dir, "Output directory")
    _must_exist(args.motif_db, "Motif DB")
    if args.raw_ms2:
        _must_exist(args.raw_ms2, "Raw MS2")
    if args.fmt_ms2:
        _must_exist(args.fmt_ms2, "Formatted MS2")
    if args.fasta:
        _must_exist(args.fasta, "FASTA")
    if args.db_csv:
        _must_exist(args.db_csv, "Prebuilt DB CSV")
        _must_exist(args.target_list, "Target list CSV")

    # Sibling .mzML
    ms2_path = args.raw_ms2 or args.fmt_ms2
    mzml = sibling_mzml(ms2_path)
    if mzml is None:
        raise FileNotFoundError(
            f"Corresponding .mzML file does not exist for {ms2_path}. "
            f"Expected sibling with '.mzML' extension."
        )


# ----------------------------
# begin_search (CLI)
# ----------------------------

def begin_search(args: argparse.Namespace) -> None:
    """
    Validate inputs (mirror GUI rules) then delegate to begin_search_confirmed(args).
    """
    # Presence checks (mirrors GUI validation from the Tk version)
    err = require_truthy([
        ("mz_min",       args.mz_min,          'Input minimum m/z value'),
        ("mz_max",       args.mz_max,          'Input maximum m/z value'),
        ("min_int",      args.min_intensity,   'Input minimum intensity value'),
        ("max_prec_z",   args.max_precursor_z, 'Input maximum precursor charge value'),
        ("max_frag_z",   args.max_fragment_z,  'Input maximum fragment charge value'),
        ("prec_err",     args.precursor_err,   'Input maximum precursor error value'),
        ("frag_err",     args.fragment_err,    'Input maximum fragment error value'),
        ("max_mods",     args.max_mods,        'Input maximum # modifications per peptide (0 if none)'),
        ("motif_db",     args.motif_db,        'Select motif database'),
        ("cov_thr",      args.cov_thr,         'Input confident coverage threshold'),
        ("out_dir",      args.out_dir,         'Select output folder'),
    ])
    if err:
        raise ValueError(err)

    # XOR: spectra
    if not bool(args.raw_ms2) ^ bool(args.fmt_ms2):
        raise ValueError('Provide either --raw-ms2 OR --fmt-ms2 (exactly one).')

    # XOR: DB source
    if not bool(args.fasta) ^ bool(args.db_csv):
        raise ValueError('Provide either --fasta OR --db-csv (exactly one).')

    # mzML sibling must exist for whichever spectra we process
    ms2_path = args.raw_ms2 or args.fmt_ms2
    if sibling_mzml(ms2_path) is None:
        raise FileNotFoundError(f"Corresponding .mzML file does not exist next to: {ms2_path}")

    # XOR: threshold
    if (args.fdr is None) == (args.eg is None):
        raise ValueError('Provide either --fdr OR --eg (exactly one).')

    # If using a prebuilt DB, we need target list now
    if args.db_csv and not args.target_list:
        raise ValueError('When using --db-csv, you must also supply --target-list.')

    # Sanity: ensure motif DB exists and out_dir is creatable
    _must_exist(args.motif_db, "Motif DB")
    _ensure_dir(args.out_dir,  "Output directory")

    # (optional) verbose argument echo
    if getattr(args, "verbose", False):
        _print_args_summary(args)

    # Delegate
    begin_search_confirmed(args)

# ----------------------------
# Argument parser
# ----------------------------

def build_parser(required_single_run: bool = True) -> argparse.ArgumentParser:
    """
    If required_single_run=True: enforce required flags (single-run CLI).
    If False: no required flags; used for batch JSON where each run supplies values.
    """
    p = argparse.ArgumentParser(
        description="EndoGenius CLI — run the first-pass search pipeline (formerly GUI).",
        add_help=True
    )

    # Batch config (exists in both modes; not required)
    p.add_argument("--config", help="JSON config with defaults and a list of runs")

    # Spectra inputs (mutually exclusive)
    g_ms2 = p.add_mutually_exclusive_group(required=required_single_run)
    g_ms2.add_argument("--raw-ms2", help="Path to raw .ms2 file")
    g_ms2.add_argument("--fmt-ms2", help="Path to formatted MS2 .txt")

    # DB inputs (mutually exclusive)
    g_db = p.add_mutually_exclusive_group(required=required_single_run)
    g_db.add_argument("--fasta", help="Path to FASTA (if building DB)")
    g_db.add_argument("--db-csv", help="Path to prebuilt database CSV")

    # If using prebuilt DB, target list is required (we validate in code, not parser)
    p.add_argument("--target-list", help="Target peptide list CSV (required with --db-csv)")

    # General paths / params
    p.add_argument("--out-dir", required=required_single_run, help="Output parent directory")
    p.add_argument("--motif-db", required=required_single_run, help="Path to motif database CSV")

    # Numeric knobs (defaults mirror your GUI)
    p.add_argument("--mz-min", default="50")
    p.add_argument("--mz-max", default="3000")
    p.add_argument("--min-intensity", default="1000")
    p.add_argument("--max-precursor-z", default="8")
    p.add_argument("--max-fragment-z", default="4")
    p.add_argument("--precursor-err", default="30")     # ppm?
    p.add_argument("--fragment-err",  default="0.05")   # Da
    p.add_argument("--max-mods", default="3")
    p.add_argument("--cov-thr", default="50")
    p.add_argument("--standard-err", default="0.1")
    p.add_argument("--min-motif-len", default="3")
    p.add_argument("--max-adjacent-swapped-AAs", default="2")
    p.add_argument("--max-swapped-AA", default="1")

    # Thresholds (mutually exclusive)
    g_thr = p.add_mutually_exclusive_group(required=required_single_run)
    g_thr.add_argument("--fdr", type=float, help="FDR cutoff (e.g., 0.05)")
    g_thr.add_argument("--eg",  type=float, help="EndoGenius score cutoff (e.g., 1000)")

    # Mod flags
    p.add_argument("--amid",            action="store_true", default=True)
    p.add_argument("--ox-m",            action="store_true", default=True)
    p.add_argument("--pg-e",            action="store_true", default=True)
    p.add_argument("--pg-q",            action="store_true", default=True)
    p.add_argument("--sodi-d",          action="store_true", default=True)

    p.add_argument("--acet-nterm",      action="store_true")
    p.add_argument("--acet-k",          action="store_true")
    p.add_argument("--acet-s",          action="store_true")
    p.add_argument("--acet-t",          action="store_true")
    p.add_argument("--biotin-k",        action="store_true")
    p.add_argument("--carb-nterm",      action="store_true")
    p.add_argument("--carb-c",          action="store_true")
    p.add_argument("--carb-d",          action="store_true")
    p.add_argument("--carb-e",          action="store_true")
    p.add_argument("--carb-h",          action="store_true")
    p.add_argument("--carb-k",          action="store_true")
    p.add_argument("--carboxy-d",       action="store_true")
    p.add_argument("--carboxy-e",       action="store_true")
    p.add_argument("--carboxy-k",       action="store_true")
    p.add_argument("--carboxy-w",       action="store_true")
    p.add_argument("--deamid-n",        action="store_true")
    p.add_argument("--deamid-q",        action="store_true")
    p.add_argument("--dehyd-d",         action="store_true")
    p.add_argument("--dehyd-s",         action="store_true")
    p.add_argument("--dehyd-t",         action="store_true")
    p.add_argument("--dehyd-y",         action="store_true")
    p.add_argument("--methyl-k",        action="store_true")
    p.add_argument("--methyl-r",        action="store_true")
    p.add_argument("--sodi-e",          action="store_true")
    p.add_argument("--plex12-nterm",    action="store_true")
    p.add_argument("--plex12-k",        action="store_true")
    p.add_argument("--mddileu-1101-nterm", action="store_true")
    p.add_argument("--mddileu-1101-k",  action="store_true")
    p.add_argument("--mddileu-0400-nterm", action="store_true")
    p.add_argument("--mddileu-0400-k",  action="store_true")
    p.add_argument("--sulfo-y",         action="store_true")
    p.add_argument("--phospho-s",       action="store_true")
    p.add_argument("--phospho-t",       action="store_true")
    p.add_argument("--phospho-y",       action="store_true")
    p.add_argument("--verbose", action="store_true", help="Print resolved arguments before each run")
    p.add_argument("--dry-run", action="store_true", help="Validate inputs and exit without running the pipeline")
    return p

# ----------------------------
# Entrypoint
# ----------------------------
def _boolify_mods(args: argparse.Namespace, mods: Dict[str, bool]) -> None:
    """
    Apply the boolean 'mods' dict to argparse Namespace in-place.
    Unknown keys are ignored; known ones set to True/False.
    """
    for k, v in mods.items():
        if hasattr(args, k.replace("-", "_")):
            setattr(args, k.replace("-", "_"), bool(v))

def _ns_from_defaults_and_run(parser: argparse.ArgumentParser,
                              defaults: Dict,
                              run: Dict) -> argparse.Namespace:
    """
    Build an argparse.Namespace for a single run by:
    1) starting with parser defaults,
    2) overlaying JSON defaults,
    3) overlaying the run object.
    """
    # base from parser defaults
    base = vars(parser.parse_args([]))  # no argv -> defaults only
    merged = deepcopy(base)

    # overlay defaults (non-mods first)
    for k, v in defaults.items():
        if k == "mods":
            continue
        merged[k] = v

    # overlay run (non-mods first)
    for k, v in run.items():
        if k in ("mods", "name"):
            continue
        merged[k] = v

    # turn back into Namespace
    ns = argparse.Namespace(**merged)

    # apply mod booleans: defaults.mods then run.mods (run overrides)
    if "mods" in defaults and isinstance(defaults["mods"], dict):
        _boolify_mods(ns, defaults["mods"])
    if "mods" in run and isinstance(run["mods"], dict):
        _boolify_mods(ns, run["mods"])

    return ns

def run_config_batch(config_path: str, parser: argparse.ArgumentParser) -> None:
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    continue_on_error = bool(cfg.get("continue_on_error", False))
    defaults = cfg.get("defaults", {})
    runs = cfg.get("runs", [])
    if not isinstance(runs, list) or not runs:
        raise ValueError("Config has no 'runs' array with at least one item.")

    for i, run in enumerate(runs, start=1):
        name = run.get("name") or f"run_{i}"
        print(f"\n=== [{i}/{len(runs)}] Starting: {name} ===")
        try:
            args = _ns_from_defaults_and_run(parser, defaults, run)
            args.name = name
            if getattr(args, "verbose", False):
                _print_args_summary(args)
            begin_search(args)
            print(f"=== [{i}/{len(runs)}] Finished: {name} ===")
        except Exception as e:
            print(f"=== [{i}/{len(runs)}] FAILED: {name} :: {e}", file=sys.stderr)
            if not continue_on_error:
                raise
def main(argv: List[str] | None = None) -> None:
    # lightweight parser to just see if --config was supplied
    peek = argparse.ArgumentParser(add_help=False)
    peek.add_argument("--config")
    known, _ = peek.parse_known_args(argv)

    # batch -> relaxed parser; single-run -> strict parser
    parser = build_parser(required_single_run=not bool(known.config))
    args = parser.parse_args(argv)

    try:
        if args.config:
            # Use a relaxed parser for defaults inside the batch builder
            batch_parser = build_parser(required_single_run=False)
            run_config_batch(args.config, batch_parser)
        else:
            begin_search(args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 15:56:23 2025

@author: lafields2
"""

from pathlib import Path
from tkinter.filedialog import askopenfilename
from tkinter import filedialog
import pymsgbox
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage, Checkbutton, messagebox, Menu, Toplevel, Listbox, StringVar, IntVar, END
import pandas as pd
import csv
import webbrowser
import os
import subprocess
from typing import Callable, Iterable, Tuple, Dict, List, Optional
from utils import *
import tkinter as tk
from tkinter import ttk

def save_csv(df: pd.DataFrame, path: str) -> None:
    """
    Save a DataFrame to CSV, creating parent directories if needed.

    Parameters
    ----------
    df : pandas.DataFrame
        Data to write.
    path : str
        Destination CSV path.

    Notes
    -----
    Overwrites the file if it already exists.
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

def pick_file(title: str, patterns: Tuple[str, Tuple[str, ...]]) -> str:
    """
    Open a file-selection dialog and return the chosen path.

    Parameters
    ----------
    title : str
        Window title for the dialog.
    patterns : tuple[str, tuple[str, ...]]
        A (label, (pattern, ...)) tuple for tk filetypes, e.g.
        ("CSV Files", ("*.csv",)).

    Returns
    -------
    str
        Selected absolute path, or an empty string if cancelled.
    """
    return askopenfilename(title=title, filetypes=[patterns]) or ""

def pick_dir(title: str = "Select folder") -> str:
    """
    Open a directory-selection dialog.

    Parameters
    ----------
    title : str, optional
        Dialog title, by default "Select folder".

    Returns
    -------
    str
        Selected directory path, or an empty string if cancelled.
    """
    return filedialog.askdirectory(title=title) or ""

def sibling_mzml(path_like: str) -> Optional[str]:
    """
    Infer and validate a sibling .mzML path for a given spectra path.

    Parameters
    ----------
    path_like : str
        Path pointing to .ms2/.txt/_formatted.* file.

    Returns
    -------
    str or None
        Existing .mzML sibling path if found; otherwise None.
    """
    repls = [("_formatted.ms2", ".mzML"), ("_formatted.txt", ".mzML"),
             (".txt", ".mzML"), (".ms2", ".mzML")]
    for old, new in repls:
        if old in path_like:
            candidate = path_like.replace(old, new)
            return candidate if os.path.isfile(candidate) else None
    return None

def require_truthy(items: Iterable[Tuple[str, str, str]]) -> Optional[str]:
    """
    Validate that a set of string values are non-empty.

    Parameters
    ----------
    items : Iterable[tuple[str, str, str]]
        Triplets of (label, value, error_message).

    Returns
    -------
    str or None
        The first error message for an empty value, or None if all pass.
    """
    for label, val, msg in items:
        if not str(val).strip():
            return msg
    return None

# ----------------------------
# Paths & Window
# ----------------------------

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path("./assets")

def relative_to_assets(path: str) -> Path:
    """
    Resolve an asset path relative to the application assets directory.

    Parameters
    ----------
    path : str
        Relative asset path (filename or subpath).

    Returns
    -------
    pathlib.Path
        Resolved absolute path inside ./assets.
    """
    return ASSETS_PATH / Path(path)

# ----------------------------
# TK Variables
# ----------------------------

input_path_MS2 = StringVar()
input_path_format_MS2 = StringVar()
mz_range_min = StringVar()
mz_range_max = StringVar()
min_intensity = StringVar()
max_precursor_z = StringVar()
max_fragment_z = StringVar()
database_csv_path = StringVar()
target_peptide_list_path = StringVar()
fasta_path = StringVar()
precursor_err = StringVar()
fragment_err = StringVar()
max_mods_pep = StringVar()
min_motif_len = StringVar(); min_motif_len.set('3')
amid_status_var = IntVar()
acet_term_status_var = IntVar()
acet_k_status_var = IntVar()
acet_s_status_var = IntVar()
acet_t_status_var = IntVar()
biotin_k_status_var = IntVar()
carb_term_status_var = IntVar()
carb_c_status_var = IntVar()
carb_d_status_var = IntVar()
carb_e_status_var = IntVar()
carb_h_status_var = IntVar()
carb_k_status_var = IntVar()
carboxy_d_status_var = IntVar()
carboxy_e_status_var = IntVar()
carboxy_k_status_var = IntVar()
carboxy_w_status_var = IntVar()
deamid_n_status_var = IntVar()
deamid_q_status_var = IntVar()
dehyd_d_status_var = IntVar()
dehyd_s_status_var = IntVar()
dehyd_t_status_var = IntVar()
dehyd_y_status_var = IntVar()
meth_k_status_var = IntVar()
meth_r_status_var = IntVar()
ox_m_status_var = IntVar()
sodi_d_status_var = IntVar()
sodi_e_status_var = IntVar()
plex12_nterm_status_var = IntVar()
plex12_k_status_var = IntVar()
pg_E_status_var = IntVar()
pg_q_status_var = IntVar()
sulfo_y_status_var = IntVar()
phospho_s_status_var = IntVar()
phospho_t_status_var = IntVar()
phospho_y_status_var = IntVar()
mddileu_1101_nterm_status_var = IntVar()
mddileu_1101_k_status_var = IntVar()
mddileu_0400_nterm_status_var = IntVar()
mddileu_0400_k_status_var = IntVar()

amid_status_var.set(1)
ox_m_status_var.set(1)
pg_E_status_var.set(1)
pg_q_status_var.set(1)
sodi_d_status_var.set(1)
#sodi_e_status_var.set(1)
# plex12_nterm_status_var.set(1)
# plex12_k_status_var.set(1)

motif_db_path = StringVar()
confident_coverage_threshold = StringVar()
standard_err = StringVar(); standard_err.set('0.1')
max_adjacent_swapped_AAs = StringVar(); max_adjacent_swapped_AAs.set('2')
FDR_threshold = StringVar()
max_swapped_AA = StringVar(); max_swapped_AA.set('1')
output_dir_path = StringVar()
eg_threshold = StringVar()

# ------------- User defaults -------------
input_path_MS2.set(r"D:\Manuscripts\2025_SodiumAdducts\QE_raw_data\30sal_TG_TR3.ms2")
fasta_path.set(r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\decoy_duplicate_removed_crustacean_database_validated_formatted20220725.fasta")
motif_db_path.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_input\fasta_raw\motif_db_20230621.csv")
output_dir_path.set(r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\sodium_adduct_addition\30sal_TG_TR3")
#database_csv_path.set(r"D:\Manuscripts\2024_Multiplexed_Feeding\LumosRun_20250830\EG_search_out\Shuffle_database.csv")
#target_peptide_list_path.set(r"D:\Manuscripts\2024_Multiplexed_Feeding\LumosRun_20250830\EG_search_out\target_list.csv")
#FDR_threshold.set('0.05')

mz_range_min.set('50')
mz_range_max.set('3000')
min_intensity.set('1000')
max_precursor_z.set('8')
max_fragment_z.set('4')
precursor_err.set('30')
fragment_err.set('0.05')
max_mods_pep.set('3')
confident_coverage_threshold.set('50')
eg_threshold.set('1000')

# ----------------------------
# File pickers
# ----------------------------
def raw_MS2_path():
    input_path_MS2.set(pick_file("Select .MS2 file", (".MS2 Files", ("*.MS2","*.ms2"))))

def formatted_MS2_path():
    input_path_format_MS2.set(pick_file("Select formatted MS2 (.txt)", ("Text Files", ("*.txt",))))

def prebuilt_db_path():
    database_csv_path.set(pick_file("Select database CSV", ("CSV Files", ("*.csv",))))

def target_list_path_get():
    target_peptide_list_path.set(pick_file("Select target list CSV", ("CSV Files", ("*.csv",))))

def fasta_db_get():
    fasta_path.set(pick_file("Select FASTA", ("FASTA Files", ("*.fasta",))))

def motif_db_get():
    motif_db_path.set(pick_file("Select motif DB CSV", ("CSV Files", ("*.csv",))))

def output_path_get():
    output_dir_path.set(pick_dir("Select output folder"))

# ----------------------------
# Make Target List 
# ----------------------------
def make_target_list():
    """
    Build and save a target peptide list from a FASTA file.

    Returns
    -------
    str
        Path to the saved `target_list.csv`.

    Notes
    -----
    Uses `utils.fasta_to_df` to produce a sequence list and writes a single
    'Sequence' column CSV in the selected output folder.
    """
    fasta_path_get = fasta_path.get()
    output_folder = output_dir_path.get()
    target_fasta = fasta_to_df(fasta_path_get)

    target_list = pd.DataFrame({'Sequence': target_fasta})
    file_path = os.path.join(output_folder, 'target_list.csv')
    save_csv(target_list, file_path)
    return file_path

# ----------------------------
# Make DB
# ----------------------------
def make_db():
    """
    Construct a database CSV from FASTA and selected modifications.

    Returns
    -------
    str
        Path to the generated database CSV.

    Notes
    -----
    Collects PTM flags from global Tk IntVars into `variable_mod_dict`,
    then calls `db_generator.make_a_DB`. The maximum modifications per
    peptide is read from the `max_mods_pep` entry.
    """
    fasta_path_get = fasta_path.get()
    output_folder = output_dir_path.get()

    MOD_RULES: List[Tuple[IntVar, Callable[[Dict], None]]] = [
        (amid_status_var,               lambda d: d.__setitem__('-Amidated', True)),
        (acet_term_status_var,          lambda d: d.__setitem__('Acetylation-', True)),
        (acet_k_status_var,             lambda d: d.setdefault('(Acetylation)', []).append('K')),
        (acet_s_status_var,             lambda d: d.setdefault('(Acetylation)', []).append('S')),
        (acet_t_status_var,             lambda d: d.setdefault('(Acetylation)', []).append('T')),
        (biotin_k_status_var,           lambda d: d.setdefault('(Biotinylation)', []).append('K')),
        (carb_term_status_var,          lambda d: d.__setitem__('Carbamidomethylation-', True)),
        (carb_c_status_var,             lambda d: d.setdefault('(Carbamidomethylation)', []).append('C')),
        (carb_d_status_var,             lambda d: d.setdefault('(Carbamidomethylation)', []).append('D')),
        (carb_e_status_var,             lambda d: d.setdefault('(Carbamidomethylation)', []).append('E')),
        (carb_h_status_var,             lambda d: d.setdefault('(Carbamidomethylation)', []).append('H')),
        (carb_k_status_var,             lambda d: d.setdefault('(Carbamidomethylation)', []).append('K')),
        (carboxy_d_status_var,          lambda d: d.setdefault('(Carboxylation)', []).append('D')),
        (carboxy_e_status_var,          lambda d: d.setdefault('(Carboxylation)', []).append('E')),
        (carboxy_k_status_var,          lambda d: d.setdefault('(Carboxylation)', []).append('K')),
        (carboxy_w_status_var,          lambda d: d.setdefault('(Carboxylation)', []).append('W')),
        (deamid_n_status_var,           lambda d: d.setdefault('(Deamidation)', []).append('N')),
        (deamid_q_status_var,           lambda d: d.setdefault('(Deamidation)', []).append('Q')),
        (dehyd_d_status_var,            lambda d: d.setdefault('(Dehydration)', []).append('D')),
        (dehyd_s_status_var,            lambda d: d.setdefault('(Dehydration)', []).append('S')),
        (dehyd_t_status_var,            lambda d: d.setdefault('(Dehydration)', []).append('T')),
        (dehyd_y_status_var,            lambda d: d.setdefault('(Dehydration)', []).append('Y')),
        (meth_k_status_var,             lambda d: d.setdefault('(Methylation)', []).append('K')),
        (meth_r_status_var,             lambda d: d.setdefault('(Methylation)', []).append('R')),
        (ox_m_status_var,               lambda d: d.setdefault('(Oxidation)', []).append('M')),
        (sodi_d_status_var,             lambda d: d.setdefault('(SodiumAdduct)', []).append('D')),
        (sodi_e_status_var,             lambda d: d.setdefault('(SodiumAdduct)', []).append('E')),
        (plex12_nterm_status_var,       lambda d: d.__setitem__('12PlexDiLeu-', True)),
        (plex12_k_status_var,           lambda d: d.setdefault('(12PlexDiLeu)', []).append('K')),
        (mddileu_1101_nterm_status_var, lambda d: d.__setitem__('mdDiLeu1101-', True)),
        (mddileu_1101_k_status_var,     lambda d: d.setdefault('(mdDiLeu1101)', []).append('K')),
        (mddileu_0400_nterm_status_var, lambda d: d.__setitem__('mdDiLeu0400-', True)),
        (mddileu_0400_k_status_var,     lambda d: d.setdefault('(mdDiLeu0400)', []).append('K')),
        (pg_E_status_var,               lambda d: d.__setitem__('(Glu->pyro-Glu)', ['ntermE'])),
        (pg_q_status_var,               lambda d: d.__setitem__('(Gln->pyro-Glu)', ['ntermQ'])),
        (sulfo_y_status_var,            lambda d: d.setdefault('(Sulfation)', []).append('Y')),
        (phospho_s_status_var,          lambda d: d.setdefault('(Phosphorylation)', []).append('S')),
        (phospho_t_status_var,          lambda d: d.setdefault('(Phosphorylation)', []).append('T')),
        (phospho_y_status_var,          lambda d: d.setdefault('(Phosphorylation)', []).append('Y'))]

    variable_mod_dict: Dict[str, object] = {}
    for var, apply_rule in MOD_RULES:
        if var.get() == 1:
            apply_rule(variable_mod_dict)

    from db_generator import make_a_DB
    max_mods_number = int(max_mods_pep.get())
    database_generated = make_a_DB(variable_mod_dict, fasta_path_get, output_folder, max_mods_number)
    return database_generated
   

# ----------------------------
# Search pipeline
# ----------------------------
def checked_clear_begin_search(predefined_db_path,output_parent_directory,raw_file_formatted_path,target_path):
    """
    Run the EndoGenius first-pass search pipeline.

    Parameters
    ----------
    predefined_db_path : str
        Path to the database CSV to search.
    output_parent_directory : str
        Parent directory for all outputs.
    raw_file_formatted_path : str
        Path to formatted spectra TXT.
    target_path : str
        Path to target peptide list CSV.

    Workflow
    --------
    Executes, in order:
      - DB search (pt1)
      - PSM assignment
      - Motif search
      - Metric extraction
      - Metric handling (score application)
      - FDR or EndoGenius score filtering (based on user choice)

    Side Effects
    ------------
    Writes multiple intermediate and final CSVs/figures into the sample
    output directory; shows status messages in the console.
    """
    from database_search import raw_file_detail_extraction, launch_db_search_pt1
    from PSM_assignment import PSM_assignment_execute
    from motif_search import start_motif_search
    from results_metric_extract import results_metric_extract_start
    from metric_handling import metric_handling_apply
    from target_decoy_assess import target_decoy_apply
    from EndoGenius_Score_Apply import endogenius_apply
    
    details = raw_file_detail_extraction(raw_file_formatted_path,output_parent_directory)
    sample_name = details[0]
    sample_output_directory = details[1]

    precursor_error_cutoff = float(precursor_err.get())
    fragment_error_cutoff = float(fragment_err.get())
    min_mz = float(mz_range_min.get())
    min_intensity_pass = int(min_intensity.get())
    standard_err_percent = float(standard_err.get())

    # Gather mod status flags
    mods = [
        amid_status_var, acet_term_status_var, acet_k_status_var, acet_s_status_var, acet_t_status_var,
        biotin_k_status_var, carb_term_status_var, carb_c_status_var, carb_d_status_var, carb_e_status_var,
        carb_h_status_var, carb_k_status_var, carboxy_d_status_var, carboxy_e_status_var, carboxy_k_status_var,
        carboxy_w_status_var, deamid_n_status_var, deamid_q_status_var, dehyd_d_status_var, dehyd_s_status_var,
        dehyd_t_status_var, dehyd_y_status_var, meth_k_status_var, meth_r_status_var, ox_m_status_var,
        sodi_d_status_var, sodi_e_status_var, plex12_nterm_status_var, plex12_k_status_var,
        mddileu_1101_nterm_status_var, mddileu_1101_k_status_var, mddileu_0400_nterm_status_var, mddileu_0400_k_status_var,
        pg_E_status_var, pg_q_status_var, sulfo_y_status_var, phospho_s_status_var, phospho_t_status_var, phospho_y_status_var]

    mod_values = [m.get() for m in mods]

    max_modifications = int(max_mods_pep.get())
    choose_mzml_directory = input_path_format_MS2.get() if len(input_path_format_MS2.get()) > 0 else input_path_MS2.get()

    first_pass_db = launch_db_search_pt1(
        predefined_db_path, output_parent_directory, choose_mzml_directory, raw_file_formatted_path,
        precursor_error_cutoff, fragment_error_cutoff, min_mz, min_intensity_pass, standard_err_percent,
        *mod_values, max_modifications, sample_output_directory)
    print('First pass DB complete')

    confident_seq_cov = float(confident_coverage_threshold.get())
    max_adjacent_swapped_AA_get = int(max_adjacent_swapped_AAs.get())
    min_motif_len_get = int(min_motif_len.get())
    num_sub_AAs = int(max_swapped_AA.get())
    motif_path = motif_db_path.get()

    first_pass_PSM_assign = PSM_assignment_execute(
        standard_err_percent, confident_seq_cov, max_adjacent_swapped_AA_get, min_motif_len_get,
        fragment_error_cutoff, num_sub_AAs, output_parent_directory, target_path, motif_path, sample_output_directory)
    print('First pass PSM assign complete')
    first_pass_motif_search = start_motif_search(output_parent_directory, motif_path, sample_output_directory)
    print('First pass motif search complete')
    first_pass_weighting_extract = results_metric_extract_start(output_parent_directory, output_parent_directory, sample_output_directory)
    print('First pass weighting extract complete')
    first_pass_metric_apply = metric_handling_apply(first_pass_weighting_extract, output_parent_directory, sample_output_directory)
    print('First pass metric apply complete')

    fdr_cutoff = str(FDR_threshold.get())
    eg_cutoff = str(eg_threshold.get())
    if len(fdr_cutoff)>0 and len(eg_cutoff)>0:
        messagebox.showerror('Input Error', 'Must input either an FDR cutoff \n or EndoGenius Score Cutoff\n not both')

    if len(fdr_cutoff)>0:
        fdr_cutoff_float = float(FDR_threshold.get())
        target_decoy_apply(first_pass_metric_apply, target_path, output_parent_directory, fdr_cutoff_float, sample_output_directory, raw_file_formatted_path)
    elif len(eg_cutoff)>0:
        eg_cutoff_float = float(eg_threshold.get())
        endogenius_apply(first_pass_metric_apply, target_path, output_parent_directory, eg_cutoff_float, sample_output_directory, raw_file_formatted_path)
    else:
        messagebox.showerror('Input Error', 'Must input either an FDR cutoff \n or EndoGenius Score Cutoff')
    print('First pass target-decoy complete')


# ----------------------------
# Begin search
# ----------------------------
def begin_search_confirmed():
    """
    Begin the search using current UI values (assumes validation done).

    Notes
    -----
    If a prebuilt DB is provided, requires a target list. Otherwise builds a
    database and target list from FASTA. Formats raw MS2 if needed. Calls
    `checked_clear_begin_search` and shows completion alerts.
    """
    predefined_db_path = database_csv_path.get()
    if len(database_csv_path.get()) > 0:
        if len(target_peptide_list_path.get()) == 0:
            messagebox.showerror('Input Error', 'Target list not specified')
        if len(target_peptide_list_path.get()) > 0:
            if len(input_path_format_MS2.get()) > 0:
                output_parent_directory = output_dir_path.get()
                raw_file_formatted_path = input_path_format_MS2.get()
                target_path = target_peptide_list_path.get()
                checked_clear_begin_search(predefined_db_path,output_parent_directory,raw_file_formatted_path,target_path)
                pymsgbox.alert('Your database search is complete','Status Update')
            elif len(input_path_format_MS2.get()) == 0:
                if len(input_path_MS2.get()) == 0:
                    messagebox.showerror('Input Error', 'No spectra file indicated')
                if len(input_path_MS2.get()) > 0:
                    from format_MS2_file_RT_IIT import format_raw_MS2
                    unformatted_spectra_path = input_path_MS2.get()
                    output_parent_directory = output_dir_path.get()
                    raw_file_formatted_path = format_raw_MS2(unformatted_spectra_path,output_parent_directory)
                    target_path = target_peptide_list_path.get()
                    checked_clear_begin_search(predefined_db_path,output_parent_directory,raw_file_formatted_path,target_path)
                    pymsgbox.alert('Your database search is complete','Status Update')
    if len(database_csv_path.get()) == 0:
        new_db_path = make_db()
        if len(input_path_format_MS2.get()) > 0:
            output_parent_directory = output_dir_path.get()
            raw_file_formatted_path = input_path_format_MS2.get()
            target_path = make_target_list()
            checked_clear_begin_search(new_db_path,output_parent_directory,raw_file_formatted_path,target_path)
            pymsgbox.alert('Your database search is complete','Status Update')
        if len(input_path_format_MS2.get()) == 0:
            from format_MS2_file_RT_IIT import format_raw_MS2
            unformatted_spectra_path = input_path_MS2.get()
            output_parent_directory = output_dir_path.get()
            raw_file_formatted_path = format_raw_MS2(unformatted_spectra_path,output_parent_directory)
            target_path = make_target_list()
            checked_clear_begin_search(new_db_path,output_parent_directory,raw_file_formatted_path,target_path)
            pymsgbox.alert('Your database search is complete','Status Update')

# ----------------------------
# Begin search
# ----------------------------
def begin_search():
    """
    Validate user inputs and kick off the search pipeline.

    Validation
    ----------
    - m/z bounds, intensity, charge limits, and error thresholds present
    - Exactly one of raw MS2 or formatted MS2 provided
    - Exactly one of FASTA or database CSV provided
    - Sibling .mzML exists for the chosen spectra path
    - Exactly one of FDR or EndoGenius score threshold provided
    - Target list required if using a prebuilt DB

    Side Effects
    ------------
    Shows error dialogs on invalid input; otherwise calls
    `begin_search_confirmed()`.
    """
    fields = {
        "mz_min": mz_range_min.get(),
        "mz_max": mz_range_max.get(),
        "min_int": min_intensity.get(),
        "max_prec_z": max_precursor_z.get(),
        "max_frag_z": max_fragment_z.get(),
        "prec_err": precursor_err.get(),
        "frag_err": fragment_err.get(),
        "max_mods": max_mods_pep.get(),
        "motif_db": motif_db_path.get(),
        "cov_thr": confident_coverage_threshold.get(),
        "out_dir": output_dir_path.get(),
        "raw_ms2": input_path_MS2.get(),
        "fmt_ms2": input_path_format_MS2.get(),
        "fasta": fasta_path.get(),
        "db_csv": database_csv_path.get(),
        "fdr": FDR_threshold.get(),
        "eg": eg_threshold.get(),
        "target_list": target_peptide_list_path.get()}

    err = require_truthy([
        ("mz_min", fields["mz_min"], 'Input minimum m/z value'),
        ("mz_max", fields["mz_max"], 'Input maximum m/z value'),
        ("min_int", fields["min_int"], 'Input minimum intensity value'),
        ("max_prec_z", fields["max_prec_z"], 'Input maximum precursor charge value'),
        ("max_frag_z", fields["max_frag_z"], 'Input maximum fragment charge value'),
        ("prec_err", fields["prec_err"], 'Input maximum precursor error value'),
        ("frag_err", fields["frag_err"], 'Input maximum fragment error value'),
        ("max_mods", fields["max_mods"], 'Input maximum # modifications per peptide\n\nIf no modifications selected, input 0'),
        ("motif_db", fields["motif_db"], 'Select motif database'),
        ("cov_thr", fields["cov_thr"], 'Input confident coverage threshold'),
        ("out_dir", fields["out_dir"], 'Select output folder')])
    if err:
        return messagebox.showerror('Input Error', err)

    have_raw = bool(fields["raw_ms2"])
    have_fmt = bool(fields["fmt_ms2"])
    if not (have_raw ^ have_fmt):
        return messagebox.showerror('Input Error', 'Input spectral file' if not (have_raw or have_fmt)
                                    else 'Input either raw or formatted spectral file, not both')

    have_fasta = bool(fields["fasta"])
    have_csvdb = bool(fields["db_csv"])
    if not (have_fasta ^ have_csvdb):
        return messagebox.showerror('Input Error', 'Input database file' if not (have_fasta or have_csvdb)
                                    else 'Input either raw .fasta file or formatted .csv database file, not both')

    ms2_path = fields["raw_ms2"] or fields["fmt_ms2"]
    if sibling_mzml(ms2_path) is None:
        return messagebox.showerror('Input Error', 'Corresponding .mzML file does not exist')

    have_fdr = bool(fields["fdr"])
    have_eg = bool(fields["eg"])
    if not (have_fdr ^ have_eg):
        return messagebox.showerror('Input Error', 'Input either FDR threshold or EndoGenius score threshold'
                                    if not (have_fdr or have_eg)
                                    else 'Input either FDR cutoff or EndoGenius score cutoff, not both')

    if have_csvdb and not fields["target_list"]:
        return messagebox.showerror('Input Error', 'Target list not specified')

    begin_search_confirmed()
    
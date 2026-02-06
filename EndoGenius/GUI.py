"""EndoGenius GUI

Tkinter-based application for running the EndoGenius proteomics pipeline.
It provides file pickers, modification selection, reporter ion extraction,
motif database tools, spectral-library building, and end-to-end search/
scoring (including FDR or EndoGenius score cutoffs).

Requires: tkinter, pandas, pymsgbox, and project-local utils.
"""

from pathlib import Path
from tkinter.filedialog import askopenfilename
from tkinter import filedialog
# import pymsgbox
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

# ----------------------------
# General Utilities
# ----------------------------

pd.options.mode.chained_assignment = None  # default='warn'
_icon_img = None


    
    
class ScrollableFrame(ttk.Frame):
    """A vertical scrollable container.

    Usage:
        sf = ScrollableFrame(root)
        sf.pack(fill="both", expand=True)
        parent = sf.scrollable_frame  # put your widgets in here
    """
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)

        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0)
        self.vbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vbar.set)

        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame_id = self.canvas.create_window(
            (0, 0), window=self.scrollable_frame, anchor="nw"
        )

        # Keep scrollregion sized to content
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        # Make the inner frame match the canvas width when the window resizes
        def _resize_inner(event):
            self.canvas.itemconfig(self.scrollable_frame_id, width=event.width)
        self.canvas.bind("<Configure>", _resize_inner)

        # Mouse wheel scrolling (Windows/macOS)
        self.canvas.bind_all("<MouseWheel>", lambda e: self.canvas.yview_scroll(int(-1*(e.delta/120)), "units"))
        # Linux (Comment out if it interferes on Windows/macOS)
        self.canvas.bind_all("<Button-4>",  lambda e: self.canvas.yview_scroll(-3, "units"))
        self.canvas.bind_all("<Button-5>",  lambda e: self.canvas.yview_scroll( 3, "units"))

        self.canvas.pack(side="left", fill="both", expand=True)
        self.vbar.pack(side="right", fill="y")

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

def set_app_icon(win):
    ico_path = relative_to_assets("endogenius.ico")
    png_path = relative_to_assets("endogenius_icon.png")
    try:
        if ico_path.exists():
            win.iconbitmap(str(ico_path))
            return
    except Exception:
        pass
    global _icon_img
    try:
        if png_path.exists():
            if _icon_img is None:
                _icon_img = PhotoImage(file=str(png_path))
            win.iconphoto(True, _icon_img)
    except Exception:
        pass

window = Tk()
# window.geometry("778x870")

window.geometry("778x600")
window.minsize(778, 600)
window.resizable(False, True) 
window.configure(bg = "#423C56")
window.title('EndoGenius v2.0.2')
set_app_icon(window)

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
sodi_e_status_var.set(1)
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

casanovo_input_mgf = StringVar()                           
casanovo_fdr = StringVar()
casanovo_model = StringVar()
casanovo_yaml = StringVar()
run_casanovo_var = tk.BooleanVar(value=False)
casanovo_db = StringVar()

# ------------- User defaults -------------
input_path_MS2.set("/Users/jenifer/Undergrad/Research/crab_database_EG/Fed_Brain_TR1.ms2")
# input_path_format_MS2.set(r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\baseline_analysis_noAdducts_noLoss_v02\EGsearch_baseline_v02\30sal_PO_TR1\round2\30sal_PO_TR1_formatted.txt")
#fasta_path.set(r"D:\Manuscripts\2025_SodiumAdducts\EG_analysis\missing_values_test\missing_values_DB_v03.fasta")
motif_db_path.set("/Users/jenifer/Undergrad/Research/crab_database_EG/motif_db_20230621.csv")
output_dir_path.set("/Users/jenifer/Undergrad/Research/crab_database_EG/test_hybrid")
database_csv_path.set("/Users/jenifer/Undergrad/Research/Shuffle_database.csv")
target_peptide_list_path.set("/Users/jenifer/Undergrad/Research/target_list.csv")
#FDR_threshold.set('0.05')

mz_range_min.set('50')
mz_range_max.set('3000')
min_intensity.set('1000')
max_precursor_z.set('8')
max_fragment_z.set('4')
precursor_err.set('20')
fragment_err.set('0.02')
max_mods_pep.set('3')
confident_coverage_threshold.set('50')
eg_threshold.set('1000')
casanovo_fdr.set("0.05")
casanovo_input_mgf.set("/Users/jenifer/Undergrad/Research/crab_database_EG/Fed_Brain_TR1.mgf")
casanovo_yaml.set("/Users/jenifer/Downloads/casanovo.yaml")
casanovo_model.set("/Users/jenifer/Downloads/casanovo_v4_2_0.ckpt")
casanovo_db.set("/Users/jenifer/Downloads/short_fasta.fasta")

# ----------------------------
# Modification selector UI
# ----------------------------

def open_mod_select_window():
    """
    Open the modification selector dialog.

    Notes
    -----
    Presents checkboxes for supported PTMs and adducts. Selections are stored
    in global Tk `IntVar`s used later by the database builder and search.
    Blocks until the user clicks Submit (destroys the window).
    """
    mod_window = Toplevel(window)
    set_app_icon(mod_window)
    mod_window.geometry("500x580")
    mod_window.configure(bg = "#D9D9D9")
    mod_window.title('Modification selector')
    
    LEFT = [('Amidation (C-termini)', amid_status_var),
        ('Acetylation (N-termini)', acet_term_status_var),
        ('Acetylation (K)', acet_k_status_var),
        ('Acetylation (S)', acet_s_status_var),
        ('Acetylation (T)', acet_t_status_var),
        ('Biotinylation (K)', biotin_k_status_var),
        ('Carbamidomethylation (N-termini)', carb_term_status_var),
        ('Carbamidomethylation (C)', carb_c_status_var),
        ('Carbamidomethylation (D)', carb_d_status_var),
        ('Carbamidomethylation (E)', carb_e_status_var),
        ('Carbamidomethylation (H)', carb_h_status_var),
        ('Carbamidomethylation (K)', carb_k_status_var),
        ('Carboxylation (D)', carboxy_d_status_var),
        ('Carboxylation (E)', carboxy_e_status_var),
        ('Carboxylation (K)', carboxy_k_status_var),
        ('Carboxylation (W)', carboxy_w_status_var),
        ('Deamidation (N)', deamid_n_status_var),
        ('Deamidation (Q)', deamid_q_status_var),
        ('Dehydration (D)', dehyd_d_status_var)]
    RIGHT = [('Dehydration (S)', dehyd_s_status_var),
        ('Dehydration (T)', dehyd_t_status_var),
        ('Dehydration (Y)', dehyd_y_status_var),
        ('Methylation (K)', meth_k_status_var),
        ('Methylation (R)', meth_r_status_var),
        ('Oxidation (M)', ox_m_status_var),
        ('Sodium adduct (D)', sodi_d_status_var),
        ('Sodium adduct (E)', sodi_e_status_var),
        ('12-plex DiLeu (N-terminal)', plex12_nterm_status_var),
        ('12-plex DiLeu (K)', plex12_k_status_var),
        ('Pyro-Glu (E)', pg_E_status_var),
        ('Pyro-Glu (Q)', pg_q_status_var),
        ('Sulfation (Y)', sulfo_y_status_var),
        ('Phosphorylation (S)', phospho_s_status_var),
        ('Phosphorylation (T)', phospho_t_status_var),
        ('Phosphorylation (Y)', phospho_y_status_var),
        ('mdDiLeu 1101 (N-term)', mddileu_1101_nterm_status_var),
        ('mdDiLeu 1101 (K)', mddileu_1101_k_status_var),
        ('mdDiLeu 0400 (N-term)', mddileu_0400_nterm_status_var),
        ('mdDiLeu 0400 (K)', mddileu_0400_k_status_var)]

    def add_checks(col_items: List[Tuple[str, IntVar]], x: int, y0: int = 5, dy: int = 25) -> None:
        for i, (label, var) in enumerate(col_items):
            Checkbutton(
                mod_window, text=label, variable=var, bg='#D9D9D9').place(x=x, y=y0 + i * dy)

    add_checks(LEFT, x=5)
    add_checks(RIGHT, x=300)

    Button(mod_window, text='Submit', borderwidth=0, highlightthickness=0,
           relief="flat", command=mod_window.destroy).place(x=200, y=500, width=78, height=30)

    mod_window.mainloop()


# ----------------------------
# Reporter Ion Extraction UI
# ----------------------------

def reporter_ion_extraction_begin():
    """
    Launch the Reporter Ion Extraction tool.

    Workflow
    --------
    1. User selects EndoGenius results CSV and spectra TXT.
    2. User confirms output folder and error tolerance (ppm).
    3. For peptides carrying 12-plex DiLeu tags, the tool finds reporter ions
       per scan, picking the highest-intensity match within the ppm window.
    4. Outputs a wide CSV of intensities per tag.

    Side Effects
    ------------
    Shows a Toplevel window and writes `reporter_ions_extracted.csv`
    to the chosen output folder.
    """
    extract_window = Toplevel(window)
    set_app_icon(extract_window)
    extract_window.geometry("750x670")
    extract_window.configure(bg = "#423C56")
    extract_window.attributes("-topmost", True)
    extract_window.title('Reporter Ion Extraction')
    
    exportdirectorypath = StringVar()
    importdirectorypath = StringVar()
    importspectrapath = StringVar()
    errorthreshold = StringVar()

    # 12plex fields
    name_vars = [StringVar() for _ in range(12)]
    mass_vars = [StringVar() for _ in range(12)]
    defaults_names = ['115a','115b','116a','116b','116c','117a','117b','117c','118a','118b','118c','118d']
    defaults_masses = ['115.12476','115.13108','116.12812','116.13444','116.14028',
                       '117.13147','117.13731','117.14363','118.13483','118.14067','118.14699','118.15283']
    for v, val in zip(name_vars, defaults_names): v.set(val)
    for v, val in zip(mass_vars, defaults_masses): v.set(val)
    
    canvas = Canvas(extract_window,bg = "#423C56",height = 800,width = 729,bd = 0,highlightthickness = 0,relief = "ridge")
    canvas.place(x = 0, y = 0)
    canvas.create_text(19.0,9.0,anchor="nw",text="Reporter Ion Extraction",fill="#FFFFFF",font=("Inter", 64 * -1))
    x, y, width, height = 19, 90, 700, 550
    canvas.create_rectangle(x, y, x+width, y+height,fill="#D9D9D9",outline="")

    def eg_results_path():
        importdirectorypath.set(pick_file("Select EndoGenius Results CSV", ("CSV Files",("*.csv",)) ))

    def eg_results_spectra_path():
        importspectrapath.set(pick_file("Select spectra TXT", ("Text Files",("*.txt",)) ))

    def export_folder():
        exportdirectorypath.set(pick_dir())

    def begin_reporter_ion_extraction():
        output_path = exportdirectorypath.get()
        results_path = importdirectorypath.get()
        spectra_path = importspectrapath.get()
        fragment_error_threshold = float(errorthreshold.get())
        mass_h = 1.00784

        tag_names = [v.get() for v in name_vars]
        tag_masses = [float(v.get()) for v in mass_vars]

        spectra = pd.read_csv(
            spectra_path, sep=",", skiprows=[0],
            names=['fragment_mz','fragment_z','fragment_resolution','precursor_mz','ms2_scan',
                   'precursor_z','precursor_RT','IonInjectTime','ms1_scan','precursor_intensity','null'])
        spectra['fragment_z'] = spectra['fragment_z'].replace(0, 1)
        spectra['fragment_monoisotopic_mass'] = (spectra['fragment_mz'] * spectra['fragment_z']) - (mass_h * spectra['fragment_z'])

        # intensity column must exist:
        if 'fragment_intensity' not in spectra.columns:
            raise ValueError("Expected 'fragment_intensity' column in spectra file.")

        results = pd.read_csv(results_path)

        # Precompute error columns
        for m in tag_masses:
            mono = m - mass_h  # z=1
            spectra[f'{m} error'] = (spectra['fragment_monoisotopic_mass'].sub(mono).abs() / mono) * 1e6

        rows = []
        mask_12plex = results['Peptide'].astype(str).str.contains(r'\(12PlexDiLeu\)', na=False)
        results_filtered = results.loc[mask_12plex, ['Peptide', 'Scan']]

        for pep, scan in results_filtered.itertuples(index=False):
            sf = spectra.loc[spectra['ms2_scan'] == scan]
            for tag_name, tag_mass in zip(tag_names, tag_masses):
                err_col = f'{tag_mass} error'
                if not sf.empty and sf[err_col].min() <= fragment_error_threshold:
                    sf_tag = sf.loc[sf[err_col] <= fragment_error_threshold]
                    if not sf_tag.empty:
                        idx = sf_tag['fragment_intensity'].idxmax()
                        rows.append((pep, scan, tag_name, sf_tag.at[idx, 'fragment_intensity'], sf_tag.at[idx, 'fragment_mz']))
                        continue
                rows.append((pep, scan, tag_name, 0, 0))

        out = pd.DataFrame(rows, columns=['Peptide','Scan','Tag','Intensity','Fragment_mz'])
        pivot_df = out.pivot_table(index=['Peptide', 'Scan'], columns='Tag', values='Intensity').reset_index()
        pivot_df.columns.name = None
        pivot_df.columns = ['Peptide', 'Scan'] + [f'Intensity_{c}' for c in pivot_df.columns[2:]]

        save_csv(pivot_df, os.path.join(output_path, 'reporter_ions_extracted.csv'))
        messagebox.showinfo("Process Complete", "Reporter ions have been extracted.")

    # --- UI wiring ---
    canvas.create_text(26.0,100.0,anchor="nw",text="EndoGenius Results Path: ",fill="#000000",font=("Inter", 16 * -1))
    eg_results_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=importdirectorypath)
    eg_results_entry.place(x=274.0,y=100.0,width=340.0,height=28.0)
    Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,relief="flat",command=eg_results_path).place(x=625.0,y=100.0,width=77.2,height=30.0)
    
    canvas.create_text(110.0,150.0,anchor="nw",text="Spectra Path: ",fill="#000000",font=("Inter", 16 * -1))
    eg_spectra_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=importspectrapath)
    eg_spectra_entry.place(x=274.0,y=150.0,width=340.0,height=28.0)
    Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,relief="flat",command=eg_results_spectra_path).place(x=625.0,y=150.0,width=77.2,height=30.0)
    
    canvas.create_text(115.0,205.0,anchor="nw",text="Export Directory: ",fill="#000000",font=("Inter", 16 * -1))
    eg_export_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=exportdirectorypath)
    eg_export_entry.place(x=274.0,y=200.0,width=340.0,height=28.0)
    Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,relief="flat",command=export_folder).place(x=625.0,y=200.0,width=77.2,height=30.0)
    
    canvas.create_text(110.0,255.0,anchor="nw",text="Error tolerance (ppm): ",fill="#000000",font=("Inter", 16 * -1))
    Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=errorthreshold).place(x=274.0,y=250.0,width=100.0,height=28.0)

    canvas.create_text(75.0,290.0,anchor="nw",text="Name",fill="#000000",font=("Inter", 16 * -1))
    canvas.create_text(215.0,290.0,anchor="nw",text="m/z",fill="#000000",font=("Inter", 16 * -1))
    canvas.create_text(400.0,290.0,anchor="nw",text="Name",fill="#000000",font=("Inter", 16 * -1))
    canvas.create_text(540.0,290.0,anchor="nw",text="m/z",fill="#000000",font=("Inter", 16 * -1))

    # place the 12 entries
    left_coords = [(50,310),(50,350),(50,390),(50,430),(50,470),(50,510)]
    right_coords = [(375,310),(375,350),(375,390),(375,430),(375,470),(375,510)]
    for i, (nx, ny) in enumerate(left_coords):
        Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=name_vars[i]).place(x=nx,y=ny,width=100,height=28)
        Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=mass_vars[i]).place(x=nx+125,y=ny,width=100,height=28)
        canvas.create_text(30.0, ny+4, anchor="nw", text=f"{i+1}:", fill="#000000", font=("Inter", 16 * -1))

    for j, (nx, ny) in enumerate(right_coords, start=6):
        Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=name_vars[j]).place(x=nx,y=ny,width=100,height=28)
        Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=mass_vars[j]).place(x=nx+125,y=ny,width=100,height=28)
        label = f"{j+1}:"
        xoff = 347.0 if j+1 >= 10 else 355.0
        canvas.create_text(xoff, ny+4, anchor="nw", text=label, fill="#000000", font=("Inter", 16 * -1))

    Button(canvas,text='Begin Analysis',borderwidth=0,highlightthickness=0,relief="flat",
           command=begin_reporter_ion_extraction).place(x=300.0,y=575.0,width=90,height=30.0)
    
    extract_window.resizable(True, True)
    extract_window.mainloop()

# ----------------------------
# MotifQuest GUI
# ----------------------------
def launch_motif_build_gui():
    """
    Launch the MotifQuest builder GUI.

    Notes
    -----
    Gathers FASTA, T-value, minimum motif length/instances, and Clustal Omega
    path, then calls `start_building_a_motif_db(...)`. Displays a completion
    message on success.
    """
    OUTPUT_PATH = Path(__file__).parent
    ASSETS_PATH = OUTPUT_PATH / Path("./assets")

    def relative_to_assets(path: str) -> Path:
        return ASSETS_PATH / Path(path)

    motif_window = Toplevel(window)
    set_app_icon(motif_window)
    motif_window.geometry("770x540")
    motif_window.configure(bg = "#423C56")
    motif_window.attributes("-topmost", True)
    motif_window.title('MotifQuest')

    min_motif_len_in = StringVar()
    min_num_motif_inst = StringVar()
    part_motif_flank = StringVar()
    fasta_path_input = StringVar()
    t_val = StringVar()
    output_dir = StringVar()
    clustal_path = StringVar()

    def build_motif_db():
        from MotifQuest_code_v02 import start_building_a_motif_db
        in_file = fasta_path_input.get()
        output_folder_path = output_dir.get()
        min_motif_len_get = int(min_motif_len_in.get())
        t_value_get = int(t_val.get())
        t_value_format = [t_value_get]
        min_motif_inst_get = int(min_num_motif_inst.get())
        clustal_omega_path = clustal_path.get()
        start_building_a_motif_db(in_file, output_folder_path, t_value_format, min_motif_len_get, min_motif_inst_get, clustal_omega_path)
        messagebox.showinfo("Finished", "MotifQuest has finished")

    def browse_files_fasta():
        fasta_path_input.set(pick_file("Select FASTA", ("FASTA Files", ("*.fasta",))))

    def browse_files_clustalo():
        clustal_path.set(pick_file("Select Clustal Omega EXE", ("Executable Files", ("*.exe",))))

    def browse_files():
        output_dir.set(pick_dir())

    motif_canvas = Canvas(motif_window, bg="#423C56", height=750, width=778, bd=0, highlightthickness=0, relief="ridge")
    motif_canvas.place(x=0, y=0)
    motif_canvas.create_rectangle(11.0, 101.0, 748.0, 517.0, fill="#D9D9D9", outline="")
    motif_canvas.create_text(17.0, 123.0, anchor="nw", text="Import FASTA Database", fill="#000000", font=("Inter", 16 * -1))
    motif_canvas.create_text(5.0, 0.0, anchor="nw", text="MotifQuest", fill="#FFFFFF", font=("Inter", 64 * -1))
    Entry(motif_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=fasta_path_input).place(x=213.0, y=123.0, width=340.0, height=28.0)
    Entry(motif_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=min_motif_len_in).place(x=213.0, y=209.0, width=340.0, height=28.0)
    Entry(motif_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=min_num_motif_inst).place(x=213.0, y=259.0, width=200.0, height=28.0)
    Entry(motif_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=part_motif_flank).place(x=213.0, y=308.0, width=340.0, height=28.0)
    Entry(motif_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=clustal_path).place(x=213.0, y=400.0, width=340.0, height=28.0)
    Button(motif_window, text='Browse', borderwidth=0, highlightthickness=0, command=browse_files_fasta, relief="flat").place(x=575.0, y=123.0, width=78.2, height=30.0)
    motif_canvas.create_text(37.0, 360.0, anchor="nw", text="Export Directory:", fill="#000000", font=("Inter", 16 * -1))
    Entry(motif_canvas, bd=0, highlightthickness=0, textvariable=output_dir).place(x=213.0, y=355.0, width=340.0, height=28.0)
    Button(motif_window, text='Browse', borderwidth=0, highlightthickness=0, command=browse_files, relief="flat").place(x=575.0, y=355.0, width=78.2, height=30.0)
    motif_canvas.create_text(33.0, 173.0, anchor="nw", text="T-Value range", fill="#000000", font=("Inter", 16 * -1))
    motif_canvas.create_text(21.0, 218.0, anchor="nw", text="Minimum Motif Length", fill="#000000", font=("Inter", 16 * -1))
    motif_canvas.create_text(15.0, 260.0, anchor="nw", text="Minimum Number of\nMotif Instances", fill="#000000", font=("Inter", 16 * -1))
    Entry(motif_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=t_val).place(x=213.0, y=166.0, width=72.0, height=28.0)
    motif_canvas.create_text(18.0, 315.0, anchor="nw", text="Partial Motif Flank Size", fill="#000000", font=("Inter", 16 * -1))
    motif_canvas.create_text(18.0, 411.0, anchor="nw", text="Clustal Omega Path", fill="#000000", font=("Inter", 16 * -1))
    Button(motif_window, text='Begin Analysis', borderwidth=0, highlightthickness=0, command=build_motif_db, relief="flat").place(x=250.0, y=450.0, width=141.0, height=30.0)
    Button(motif_window, text='Browse', borderwidth=0, highlightthickness=0, command=browse_files_clustalo, relief="flat").place(x=575.0, y=400.0, width=78.2, height=30.0)
    motif_window.resizable(False, False)
    motif_window.mainloop()

# ----------------------------
# Library Builder GUI
# ----------------------------
def launch_lib_build_gui():
    """
    Launch the Spectral Library Builder GUI.

    Notes
    -----
    Asks for an EndoGenius results directory, output directory, and fragment
    error threshold, then calls `build_a_SL(...)`. Displays a completion
    message on success.
    """
    OUTPUT_PATH = Path(__file__).parent
    ASSETS_PATH = OUTPUT_PATH / Path("./assets")
        
    def relative_to_assets(path: str) -> Path:
        return ASSETS_PATH / Path(path)
    
    lib_window = Toplevel(window)
    set_app_icon(lib_window)
    lib_window.geometry("729x295")
    lib_window.configure(bg = "#423C56")
    lib_window.attributes("-topmost", True)

    eg_results_dir = StringVar()
    out_dir = StringVar()
    error = StringVar()

    def in_dir_path_get():
        eg_results_dir.set(pick_dir("Select EndoGenius Results Directory"))
        
    def out_dir_path_get():
        out_dir.set(pick_dir("Select Output Directory"))
    
    def launch_library_process():
        eg_results_dir_get = eg_results_dir.get()
        out_dir_get = out_dir.get()
        error_get = float(error.get())
        from sl_builder_driver import build_a_SL
        build_a_SL(eg_results_dir_get,out_dir_get,error_get)
        messagebox.showinfo("Process Complete", "The library has been built.")

    canvas = Canvas(lib_window,bg = "#423C56",height = 295,width = 729,bd = 0,highlightthickness = 0,relief = "ridge")
    canvas.place(x = 0, y = 0)
    canvas.create_text(19.0,9.0,anchor="nw",text="Library Builder",fill="#FFFFFF",font=("Inter", 64 * -1))
    canvas.create_rectangle(19.0,90.0,709.0,271.0,fill="#D9D9D9",outline="")                       

    canvas.create_text(26.0,106.0,anchor="nw",text="EndoGenius Results Directory",fill="#000000",font=("Inter", 16 * -1))
    Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=eg_results_dir).place(x=274.0,y=106.0,width=340.0,height=28.0)
    Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,command=in_dir_path_get,relief="flat").place(x=625.0,y=106.0,width=77.2,height=30.0)
    
    canvas.create_text(26.0,182.0,anchor="nw",text="Output Directory",fill="#000000",font=("Inter", 16 * -1))
    Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=out_dir).place(x=274.0,y=182.0,width=340.0,height=28.0)
    
    canvas.create_text(26.0,144.0,anchor="nw",text="Fragment error threshold (Da)",fill="#000000",font=("Inter", 16 * -1))
    Entry(lib_window,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=error).place(x=274.0,y=144.0,width=69.0,height=28.0)

    Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,command=out_dir_path_get,relief="flat").place(x=625.0,y=182.0,width=77.2,height=30.0)
    Button(canvas,borderwidth=0,highlightthickness=0,text='Build Library',command=launch_library_process,relief="flat").place(x=303.0,y=228.0,width=141.0,height=30.0)

    lib_window.resizable(True, True)
    lib_window.mainloop()

def launch_diann():
    """
    Attempt to start the DIA-NN GUI.

    Side Effects
    ------------
    Spawns the DIA-NN executable as a subprocess and shows a status dialog.

    Raises
    ------
    Shows an error dialog if the executable cannot be launched.
    """
    try:
        subprocess.Popen([r"C:\DIA-NN\2.2.0\DIA-NN.exe"])
        messagebox.showinfo("Launch", "DIA-NN GUI is launching...")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to launch DIA-NN: {e}")

# ----------------------------
# Quant GUI 
# ----------------------------
def launch_quant():
    """
    Open the Quantification app window.

    Workflow
    --------
    - User adds N CSV result files and provides N matching sample names.
    - For each file, the peptide with the highest intensity per sequence is kept.
    - Merges all samples on 'Peptide' and writes `merged_intensities.csv`.

    Side Effects
    ------------
    Opens a Toplevel window and saves an output CSV to the chosen directory.
    """
    quant_window = Toplevel(window)
    set_app_icon(quant_window)
    quant_window.geometry("760x530")
    quant_window.title('Quantification App')

    quant_menubar = Menu(quant_window)
    quant_filemenu = Menu(quant_menubar, tearoff=0)
    quant_filemenu.add_command(label="Exit", command=quant_window.quit)
    quant_menubar.add_cascade(label="File", menu=quant_filemenu)
    quant_window.configure(bg = "#423C56",menu=quant_menubar)

    quant_output_dir = StringVar()
    quant_input_file_path = StringVar()
    quant_input_sample_name = StringVar()

    def quantify_summary_begin():
        quant_output_directory = quant_output_dir.get()
        quant_file_paths = list(quant_entry_4.get(0,END))
        quant_file_names = list(quant_entry_5.get(0,END))

        quant_merge_df = pd.DataFrame()

        if len(quant_file_paths) != len(quant_file_names):
            raise Exception("A file name must be selected for each file input")
        else:
            for a in range(len(quant_file_paths)):
                quant_single_file_path = quant_file_paths[a]
                quant_single_name = quant_file_names[a]
                quant_single_file = pd.read_csv(quant_single_file_path)

                quant_single_file_filtered = pd.DataFrame({
                    'Peptide': quant_single_file['Peptide'],
                    f'{quant_single_name} Intensity': quant_single_file['precursor_intensity']})

                quant_single_file_filtered = quant_single_file_filtered.sort_values(by=f'{quant_single_name} Intensity', ascending=False)
                quant_single_file_filtered = quant_single_file_filtered.drop_duplicates(subset='Peptide')

                quant_merge_df = quant_single_file_filtered if quant_merge_df.empty else pd.merge(quant_merge_df, 
                                                                                                  quant_single_file_filtered, on='Peptide', how='outer')

        out_path = os.path.join(quant_output_directory, 'merged_intensities.csv')
        save_csv(quant_merge_df, out_path)
        messagebox.showinfo('Your quantification report is complete','Status Update')

    def quant_output_path_get():
        quant_output_dir.set(pick_dir("Select export directory"))

    def quant_input_file_path_get():
        quant_input_file_path.set(pick_file("Select CSV", ("CSV Files",("*.csv",)) ))

    def quant_add_file():
        quant_entry_4.insert('end', quant_input_file_path.get())
        quant_entry_1.delete(0,'end')

    def quant_delete_file():
        if quant_entry_4.curselection():
            quant_entry_4.delete(quant_entry_4.curselection())

    def quant_add_name():
        quant_entry_5.insert('end', quant_input_sample_name.get())
        quant_entry_3.delete(0,'end')

    def quant_delete_sample_name():
        if quant_entry_5.curselection():
            quant_entry_5.delete(quant_entry_5.curselection())

    quant_canvas = Canvas(quant_window, bg="#423C56", height=883, width=778, bd=0, highlightthickness=0, relief="ridge")
    quant_canvas.place(x=0, y=0)
    quant_canvas.create_rectangle(10.0, 109.0, 747.0, 525.0, fill="#D9D9D9", outline="")
    quant_canvas.create_text(15.0, 128.0, anchor="nw", text="Import Results:", fill="#000000", font=("Inter", 16 * -1))
    quant_canvas.create_text(5.0, 0.0, anchor="nw", text="Quantify Results", fill="#FFFFFF", font=("Inter", 64 * -1))

    quant_entry_1 = Entry(quant_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=quant_input_file_path)
    quant_entry_1.place(x=131.0, y=123.0, width=340.0, height=28.0)
    Button(quant_canvas, text='Browse', borderwidth=0, highlightthickness=0, command=quant_input_file_path_get, relief="flat").place(x=480.0, y=123.0, width=78.2, height=30.0)
    
    quant_canvas.create_text(11.0, 452.0, anchor="nw", text="Export Directory:", fill="#000000", font=("Inter", 16 * -1))
    Entry(quant_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=quant_output_dir).place(x=153.0, y=447.0, width=340.0, height=28.0)
    Button(quant_canvas, text='Browse', borderwidth=0, highlightthickness=0, command=quant_output_path_get, relief="flat").place(x=502.0, y=447.0, width=78.2, height=30.0)
    
    quant_canvas.create_text(15.0, 168.0, anchor="nw", text="Sample Name:", fill="#000000", font=("Inter", 16 * -1))
    quant_canvas.create_text(14.0, 200.0, anchor="nw", text="File List:", fill="#000000", font=("Inter", 16 * -1))
    quant_entry_3 = Entry(quant_canvas, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=quant_input_sample_name)
    quant_entry_3.place(x=131.0, y=163.0, width=340.0, height=28.0)
    quant_entry_4 = Listbox(quant_canvas, width=50); quant_entry_4.place(x=17.0, y=225.0, width=340.0, height=144.0)
    quant_canvas.create_text(352.0, 200.0, anchor="nw", text="Name List:", fill="#000000", font=("Inter", 16 * -1))
    quant_entry_5 = Listbox(quant_canvas, width=50); quant_entry_5.place(x=372.0, y=225.0, width=340.0, height=144.0)
    
    Button(quant_canvas, text='Add Name', borderwidth=0, highlightthickness=0, command=quant_add_name, relief="flat").place(x=480.0, y=162.0, width=109.4, height=31.0)
    Button(quant_canvas, text='Add File', borderwidth=0, highlightthickness=0, command=quant_add_file, relief="flat").place(x=558.0, y=123.0, width=77.2, height=30.0)
    Button(quant_canvas, text='Delete File', borderwidth=0, highlightthickness=0, command=quant_delete_file, relief="flat").place(x=131.0, y=381.0, width=100.0, height=30.0)
    Button(quant_canvas, text='Delete Name', borderwidth=0, highlightthickness=0, command=quant_delete_sample_name, relief="flat").place(x=503.0, y=381.0, width=100.0, height=30.0)
    Button(quant_canvas, text='Begin', borderwidth=0, highlightthickness=0, command=quantify_summary_begin, relief="flat").place(x=319.0, y=487.0, width=141.0, height=30.0)
    quant_canvas.create_rectangle(638.0, 828.0, 738.0, 928.0, fill="#000000", outline="")
    
    quant_window.resizable(False, False)
    quant_window.mainloop()


# ----------------------------
# Simple URL helpers
# ----------------------------
def openweb_liweb():
    """
    Open the Li Lab resources webpage in the default browser.
    """
    webbrowser.open("https://www.lilabs.org/resources", new=1)

def openweb_git():
    """
    Open the Li Lab GitHub organization page in the default browser.
    """
    webbrowser.open("https://github.com/lingjunli-research", new=1)

def openweb_user_manual():
    """
    Open the EndoGenius user manual in the default browser.
    """
    webbrowser.open("https://docs.google.com/document/d/e/2PACX-1vRwKSjIl6wu88MTObZ7G0QYl9wzg7Rm065o4AxM1zzAMspEfHChLMcHMmWFWD8BjLIKSsvsqONeHknb/pub", new=1)


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

def casanovo_input_mgf_path_get():
    casanovo_input_mgf.set(pick_file("Select target mgf file", ("MGF Files", ("*.mgf"))))

def casanovo_yaml_get():
    casanovo_yaml.set(pick_file("Select yaml config file", ("YAML Files", ("*.yaml"))))

def casanovo_model_get():
    casanovo_model.set(pick_file("Select the model", ("CKPT Files", ("*.ckpt"))))

def casanovo_db_get():
    casanovo_db.set(pick_file("Select the database FASTA", ("FASTA Files", ("*.fasta"))))

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
                if run_casanovo_var.get():
                    print("Running Casanovo pipeline...")
                    run_casanovo_pipeline()
                messagebox.showinfo('Your database search is complete','Status Update')
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
                    if run_casanovo_var.get():
                        print("Running Casanovo pipeline...")
                        run_casanovo_pipeline()
                    messagebox.showinfo('Your database search is complete','Status Update')
    if len(database_csv_path.get()) == 0:
        new_db_path = make_db()
        if len(input_path_format_MS2.get()) > 0:
            output_parent_directory = output_dir_path.get()
            raw_file_formatted_path = input_path_format_MS2.get()
            target_path = make_target_list()
            checked_clear_begin_search(new_db_path,output_parent_directory,raw_file_formatted_path,target_path)
            if run_casanovo_var.get():
                print("Running Casanovo pipeline...")
                run_casanovo_pipeline()
            messagebox.showinfo('Your database search is complete','Status Update')
        if len(input_path_format_MS2.get()) == 0:
            from format_MS2_file_RT_IIT import format_raw_MS2
            unformatted_spectra_path = input_path_MS2.get()
            output_parent_directory = output_dir_path.get()
            raw_file_formatted_path = format_raw_MS2(unformatted_spectra_path,output_parent_directory)
            target_path = make_target_list()
            checked_clear_begin_search(new_db_path,output_parent_directory,raw_file_formatted_path,target_path)
            if run_casanovo_var.get():
                print("Running Casanovo pipeline...")
                run_casanovo_pipeline()
            messagebox.showinfo('Your database search is complete','Status Update')

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

# ----------------------------
# casanovo & fasta generation
# ----------------------------
def run_casanovo(input_mgf, output_mztab, casanovo_yaml, casanovo_model):
    cmd = [
        "casanovo", "sequence",
        "-o", output_mztab,
        "-c", casanovo_yaml,
        "-m", casanovo_model,
        input_mgf
    ]
    subprocess.run(cmd, check=True)  
    return(output_mztab)                                 

def run_decoy(input_mgf):
    import generate_decoy  
    decoy_mgf = os.path.join(output_dir_path.get(),"decoy_spectra.mgf")
    generate_decoy.generate_decoy_mgf(input_mgf, decoy_mgf)
    return(decoy_mgf)                                              

def run_analysis(target_file, decoy_file, db_file, fdr=0.05):
    import filter_casanovo
    filter_casanovo.automated_data_analysis(
        target_casanovo_result=target_file,
        decoy_casanovo_result=decoy_file,
        engine_score_column="search_engine_score[1]",
        required_fdr=fdr,
        output_dir=output_dir_path.get(),
        db = db_file
    )

def run_casanovo_pipeline():
    try:
        mgf_file = casanovo_input_mgf.get()
        db_file = casanovo_db.get()
        fdr_score = float(casanovo_fdr.get())
        model_path = casanovo_model.get()
        yaml_path = casanovo_yaml.get()
        output_mztab_target = os.path.join(output_dir_path.get(),"casanovo_result_target.mztab")
        output_mztab_decoy = os.path.join(output_dir_path.get(),"casanovo_result_decoy.mztab")
        # Step 1: Run Casanovo
        output_mztab = run_casanovo(mgf_file, output_mztab_target, yaml_path, model_path)

        # Step 2: Generate Decoy
        decoy_mgf = run_decoy(mgf_file)
        
        # Step 3: Run Decoy Casanovo
        decoy_file = run_casanovo(decoy_mgf, output_mztab_decoy, yaml_path, model_path)

        # Step 4: Run Analysis
        run_analysis(output_mztab, decoy_file, db_file, fdr = fdr_score)

        messagebox.showinfo("Pipeline Finished")
    except Exception as e:
        messagebox.showerror("Error", str(e))
                                                                                    
# ----------------------------
# Menus & Main Canvaxs/UI
# ----------------------------
menubar = Menu(window)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Exit", command=window.quit)
menubar.add_cascade(label="File", menu=filemenu)

helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="Li Lab Website", command=openweb_liweb)
helpmenu.add_command(label="Li Lab GitHub", command=openweb_git)
helpmenu.add_command(label="User manual", command=openweb_user_manual)
menubar.add_cascade(label="Help", menu=helpmenu)

toolmenu = Menu(menubar, tearoff=0)
toolmenu.add_command(label="Quantiation Report", command = launch_quant)
toolmenu.add_command(label="Build Spectral Library", command = launch_lib_build_gui)
toolmenu.add_command(label="Launch DIA-NN GUI", command=launch_diann)
toolmenu.add_command(label="Launch MotifQuest", command=launch_motif_build_gui)
toolmenu.add_command(label="Extract Reporter Ions", command=reporter_ion_extraction_begin)
menubar.add_cascade(label="Tools", menu=toolmenu)

window.config(menu=menubar)

main_sf = ScrollableFrame(window)
main_sf.pack(fill="both", expand=True)
PARENT = main_sf.scrollable_frame 
# ----------------------------
# Main canvas widgets
# ----------------------------

canvas = Canvas(PARENT, bg="#423C56", height=1200, width=778, bd=0, highlightthickness=0, relief="ridge")
canvas.place(x=0, y=0)

canvas.create_rectangle(10.0,200.0,773.0,257.0,fill="#D9D9D9",outline="")
canvas.create_text(20.0,216.0,anchor="nw",text="m/z range",fill="#000000",font=("Inter", 16 * -1))
canvas.create_text(171.0,219.0,anchor="nw",text="-",fill="#000000",font=("Inter", 16 * -1))
canvas.create_text(270.0,210.0,anchor="nw",text="minimum\nintensity",fill="#000000",font=("Inter", 16 * -1),justify='center')
canvas.create_text(420.0,210.0,anchor="nw",text="max precursor\ncharge",fill="#000000",font=("Inter", 16 * -1),justify='center')
canvas.create_text(600.0,210.0,anchor="nw",text="max fragment\ncharge",fill="#000000",font=("Inter", 16 * -1),justify='center')
canvas.create_text(5.0, 166.0, anchor="nw", text="2. Spectral processing", fill="#FFFFFF", font=("Inter", 24 * -1))
canvas.create_text(5.0, 76.0, anchor="nw", text="1. Spectral input", fill="#FFFFFF", font=("Inter", 24 * -1))
canvas.create_rectangle(10.0, 109.0, 369.0, 166.0, fill="#D9D9D9", outline="")
canvas.create_text(32.0, 128.0, anchor="nw", text="Raw .MS2", fill="#000000", font=("Inter", 16 * -1))
canvas.create_rectangle(411.0, 109.0, 770.0, 166.0, fill="#D9D9D9", outline="")
entry_1 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=input_path_format_MS2)
entry_1.place(x=527.0, y=123.0, width=146.0, height=28.0)
canvas.create_text(440.0, 119.0, anchor="nw", text="Formatted\nRaw .MS2", fill="#000000", font=("Inter", 16 * -1), justify='center')
canvas.create_text(332.0, 125.0, anchor="nw", text="or", fill="#FFFFFF", font=("Inter", 16 * -1))
canvas.create_text(5.0, 582.0, anchor="nw", text="5. PSM assignment", fill="#FFFFFF", font=("Inter", 24 * -1))
canvas.create_rectangle(10.0, 613.0, 770.0, 715.0, fill="#D9D9D9", outline="")
canvas.create_text(14.0, 626.0, anchor="nw", text="Motif database", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(425.0, 626.0, anchor="nw", text="Confident coverage threshold (%)", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(5.0, 720.0, anchor="nw", text="6. De novo sequencing (casanovo)", fill="#FFFFFF", font=("Inter", 24 * -1))
canvas.create_text(5.0, 847.0, anchor="nw", text="7. Export results", fill="#FFFFFF", font=("Inter", 24 * -1))
canvas.create_rectangle(10.0, 876.0, 411.0, 935.0, fill="#D9D9D9", outline="")
canvas.create_text(29.0, 899.0, anchor="nw", text="Output directory", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(5.0, 265.0, anchor="nw", text="3. Database definition", fill="#FFFFFF", font=("Inter", 24 * -1))
canvas.create_rectangle(10.0, 315.0, 369.0, 436.0, fill="#D9D9D9", outline="")
canvas.create_rectangle(10.0, 749.0, 770.0, 843.0, fill="#D9D9D9", outline="")
canvas.create_text(20.0,768,anchor="nw",text="MGF file",fill="#000000",font=("Inter", 18 * -1))
canvas.create_text(436.0,768,anchor="nw",text="FDR rate",fill="#000000",font=("Inter", 18 * -1))
canvas.create_text(32.0, 334.0, anchor="nw", text="Database", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(21.0, 388.0, anchor="nw", text="Target\npeptide list", fill="#000000", font=("Inter", 16 * -1), justify='center')
canvas.create_text(338.0, 335.0, anchor="nw", text="or", fill="#FFFFFF", font=("Inter", 16 * -1))
canvas.create_text(7.0, 293.0, anchor="nw", text="Pre-built database", fill="#FFFFFF", font=("Inter", 16 * -1))
canvas.create_text(420.0, 293.0, anchor="nw", text="Generate from .fasta", fill="#FFFFFF", font=("Inter", 16 * -1))
canvas.create_rectangle(414.0, 316.0, 773.0, 374.0, fill="#D9D9D9", outline="")
canvas.create_text(436.0, 335.0, anchor="nw", text="Database", fill="#000000", font=("Inter", 16 * -1))
canvas.create_rectangle(10.0, 466.0, 773.0, 579.0, fill="#D9D9D9", outline="")
canvas.create_text(25.0, 489.0, anchor="nw", text="Precursor error (ppm)", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(14.0, 540.0, anchor="nw", text="Modifications", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(280.0, 489.0, anchor="nw", text="Fragment error (Da)", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(25.0, 805.0, anchor="nw", text="Model", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(278.0, 805.0, anchor="nw", text="Yaml config", fill="#000000", font=("Inter", 16 * -1))
canvas.create_text(548.0, 805.0, anchor="nw", text="Filter db", fill="#000000", font=("Inter", 16 * -1))


button_modification = Button(PARENT, text='Select Modifications', borderwidth=0, highlightthickness=0, command=open_mod_select_window, relief="flat")
button_modification.place(x=125.0, y=535.0, width=150, height=30.0)
canvas.create_text(550.0, 489.0, anchor="nw", text="Max mods/peptide", fill="#000000", font=("Inter", 16 * -1), justify='center')
canvas.create_text(5.0, 436.0, anchor="nw", text="4. Database search", fill="#FFFFFF", font=("Inter", 24 * -1))
button_image_1 = PhotoImage(file=relative_to_assets("button_1.png"))
button_1 = Button(PARENT, image=button_image_1, borderwidth=0, highlightthickness=0, command=begin_search, relief="flat")
button_1.place(x=9.0, y=944.0, width=401.0, height=59.0)
image_image_0 = PhotoImage(file=relative_to_assets("EndoGenius_Logo_12.png"))
image_0 = canvas.create_image(400.0, 60.0, image=image_image_0)
image_image_1 = PhotoImage(file=relative_to_assets("image_1.png"))
image_1 = canvas.create_image(617.0, 944.0, image=image_image_1)
canvas.create_text(390.0, 668.0, anchor="nw", text="FDR\nThreshold", fill="#000000", font=("Inter", 16 * -1), justify='center')
canvas.create_text(550.0, 668.0, anchor="nw", text="EndoGenius Score\nThreshold", fill="#000000", font=("Inter", 16 * -1), justify='center')
button_image_2 = PhotoImage(file=relative_to_assets("button_2.png"))
button_2 = Button(PARENT, image=button_image_2, borderwidth=0, highlightthickness=0, command=formatted_MS2_path, relief="flat")
button_2.place(x=678.0, y=123.0, width=77.21710205078125, height=30.0)
entry_2 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=input_path_MS2)
entry_2.place(x=131.0, y=123.0, width=146.0, height=28.0)
button_image_3 = PhotoImage(file=relative_to_assets("button_3.png"))
button_3 = Button(PARENT, image=button_image_3, borderwidth=0, highlightthickness=0, command=raw_MS2_path, relief="flat")
button_3.place(x=282.0, y=123.0, width=77.21710205078125, height=30.0)
entry_3 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=output_dir_path)
entry_3.place(x=168.0, y=897.0, width=146.0, height=28.0)
button_image_4 = PhotoImage(file=relative_to_assets("button_4.png"))
button_4 = Button(PARENT, image=button_image_4, borderwidth=0, highlightthickness=0, command=output_path_get, relief="flat")
button_4.place(x=319.0, y=897.0, width=77.21710205078125, height=30.0)
entry_4 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=fasta_path)
entry_4.place(x=531.0, y=329.0, width=146.0, height=28.0)
button_image_5 = PhotoImage(file=relative_to_assets("button_5.png"))
button_5 = Button(PARENT, image=button_image_5, borderwidth=0, highlightthickness=0, command=fasta_db_get, relief="flat")
button_5.place(x=682.0, y=329.0, width=77.21710205078125, height=30.0)
entry_5 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=database_csv_path)
entry_5.place(x=121.0, y=329.0, width=146.0, height=28.0)
button_image_6 = PhotoImage(file=relative_to_assets("button_6.png"))
button_6 = Button(PARENT, image=button_image_6, borderwidth=0, highlightthickness=0, command=prebuilt_db_path, relief="flat")
button_6.place(x=272.0, y=329.0, width=77.21710205078125, height=30.0)
button_image_7 = PhotoImage(file=relative_to_assets("button_7.png"))
button_7 = Button(PARENT, image=button_image_7, borderwidth=0, highlightthickness=0, command=lambda: print("button_7 clicked"), relief="flat")
button_7.place(x=272.0, y=388.0, width=77.21710205078125, height=30.0)
entry_7 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=target_peptide_list_path)
entry_7.place(x=121.0, y=388.0, width=146.0, height=28.0)
button_image_8 = PhotoImage(file=relative_to_assets("button_8.png"))
button_8 = Button(PARENT, image=button_image_8, borderwidth=0, highlightthickness=0, command=target_list_path_get, relief="flat")
button_8.place(x=272.0, y=388.0, width=77.21710205078125, height=30.0)
button_image_10 = PhotoImage(file=relative_to_assets("button_10.png"))
button_10 = Button(PARENT, image=button_image_10, borderwidth=0, highlightthickness=0, command=casanovo_input_mgf_path_get, relief="flat")
button_10.place(x=300.0, y=765.0, width=77.21710205078125, height=30.0)
button_image_11 = PhotoImage(file=relative_to_assets("button_11.png"))
button_11 = Button(PARENT, image=button_image_10, borderwidth=0, highlightthickness=0, command=casanovo_db_get, relief="flat")
button_11.place(x=694.0, y=803.0, width=77.21710205078125, height=30.0)
button_image_12 = PhotoImage(file=relative_to_assets("button_12.png"))
button_12 = Button(PARENT, image=button_image_10, borderwidth=0, highlightthickness=0, command=casanovo_yaml_get, relief="flat")
button_12.place(x=453.0, y=803.0, width=77.21710205078125, height=30.0)
button_image_13 = PhotoImage(file=relative_to_assets("button_13.png"))
button_13 = Button(PARENT, image=button_image_10, borderwidth=0, highlightthickness=0, command=casanovo_model_get, relief="flat")
button_13.place(x=185.0, y=803.0, width=77.21710205078125, height=30.0)
entry_8 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=motif_db_path)
entry_8.place(x=141.0, y=624.0, width=146.0, height=28.0)
entry_10 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=max_fragment_z)
entry_10.place(x=712.0, y=214.0, width=53.0, height=28.0)
entry_11 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=max_precursor_z)
entry_11.place(x=537.0, y=213.0, width=53.0, height=28.0)
entry_12 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=min_intensity)
entry_12.place(x=352.0, y=213.0, width=53.0, height=28.0)
entry_13 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=mz_range_max)
entry_13.place(x=196.0, y=213.0, width=53.0, height=28.0)
entry_14 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=mz_range_min)
entry_14.place(x=106.0, y=213.0, width=53.0, height=28.0)
entry_15 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=max_mods_pep)
entry_15.place(x=685.0, y=483.0, width=53.0, height=28.0)
entry_16 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=fragment_err)
entry_16.place(x=439.0, y=483.0, width=53.0, height=28.0)
entry_17 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=precursor_err)
entry_17.place(x=189.0, y=483.0, width=53.0, height=28.0)
entry_18 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=confident_coverage_threshold)
entry_18.place(x=685.0, y=621.0, width=53.0, height=28.0)
entry_20 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=FDR_threshold)
entry_20.place(x=472.0, y=669.0, width=53.0, height=28.0)
entry_21 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=eg_threshold)
entry_21.place(x=685.0, y=669.0, width=53.0, height=28.0)
entry_22 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=casanovo_input_mgf)
entry_22.place(x=90.0, y=765.0, width=200.0, height=28.0)
entry_23 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=casanovo_fdr)
entry_23.place(x=540.0, y=765.0, width=100.0, height=28.0)
button_image_9 = PhotoImage(file=relative_to_assets("button_9.png"))
button_9 = Button(PARENT, image=button_image_9, borderwidth=0, highlightthickness=0, command=motif_db_get, relief="flat")
button_9.place(x=292.0, y=624.0, width=77.21710205078125, height=30.0)
casanovo_checkbox = tk.Checkbutton(PARENT, text="Run Casanovo", variable=run_casanovo_var, font=("Arial", 10))
casanovo_checkbox.place(x=650.0, y=765.0, width=110, height=30.0)
entry_24 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=casanovo_db)
entry_24.place(x=613.0, y=803.0, width=73.0, height=28.0)
entry_25 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=casanovo_yaml)
entry_25.place(x=370.0, y=803.0, width=78.0, height=28.0)
entry_26 = Entry(PARENT, bd=0, bg="#FFFFFF", highlightthickness=0, textvariable=casanovo_model)
entry_26.place(x=78.0, y=803.0, width=103.0, height=28.0)
canvas.create_rectangle(138.0, 544.0, 156.0, 562.0, fill="#D9D9D9", outline="")
canvas.create_rectangle(295.0, 544.0, 313.0, 562.0, fill="#D9D9D9", outline="")
canvas.create_rectangle(417.0, 544.0, 435.0, 562.0, fill="#D9D9D9", outline="")
canvas.create_rectangle(652.0, 544.0, 670.0, 562.0, fill="#D9D9D9", outline="")
canvas.create_rectangle(542.0, 544.0, 560.0, 562.0, fill="#D9D9D9", outline="")

CONTENT_W, CONTENT_H = 778, 1200
PARENT.configure(width=CONTENT_W, height=CONTENT_H)
main_sf.canvas.configure(scrollregion=(0, 0, CONTENT_W, CONTENT_H))

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    window.resizable(False, False)
    window.mainloop()
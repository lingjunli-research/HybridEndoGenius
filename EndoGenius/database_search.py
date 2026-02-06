# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:03:08 2023

@author: lawashburn
"""

from __future__ import annotations

import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
from pyopenms import *
pd.options.mode.chained_assignment = None

print('Database Search')

# ------------------------- #
# Path & output folder util #
# ------------------------- #
def raw_file_detail_extraction(raw_file_formatted_path: str, output_parent_directory: str):
    """
    Returns (sample_name, output_folder) and creates output folder if needed.

    """
    p = Path(raw_file_formatted_path)
    sample_name = p.stem.replace("_formatted", "")
    out_dir = Path(output_parent_directory) / sample_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return sample_name, str(out_dir)

# -------------------------------- #
# Physical constants & mass tables #
# -------------------------------- #
proton_mass = 1.00727646688
H = 1.0078250352
O = 15.99491463
C = 12.0000000
N = 14.003074
P = 30.973762
S = 31.9720707
Na = 22.989769
C13 = 13.003355
O18 = 17.9991610
N15 = 15.0001
H2 = 2.014101777844
h_mass = 1.00784  # used for mono-mass conversions in your IO
charge = 1 #fragment charge

# Base AA monoisotopic masses (residue formula masses)
AA_MASSES = {
    "G": C * 2 + H * 3 + N + O,
    "A": C * 3 + H * 5 + N + O,
    "S": C * 3 + H * 5 + N + O * 2,
    "P": C * 5 + H * 7 + N + O,
    "V": C * 5 + H * 9 + N + O,
    "T": C * 4 + H * 7 + N + O * 2,
    "C": C * 3 + H * 5 + N + O + S,
    "L": C * 6 + H * 11 + N + O,
    "I": C * 6 + H * 11 + N + O,
    "N": C * 4 + H * 6 + N * 2 + O * 2,
    "D": C * 4 + H * 5 + N + O * 3,
    "Q": C * 5 + H * 8 + N * 2 + O * 2,
    "K": C * 6 + H * 12 + N * 2 + O,
    "E": C * 5 + H * 7 + N + O * 3,
    "M": C * 5 + H * 9 + N + O + S,
    "H": C * 6 + H * 7 + N * 3 + O,
    "F": C * 9 + H * 9 + N + O,
    "R": C * 6 + H * 12 + N * 4 + O,
    "Y": C * 9 + H * 9 + N + O * 2,
    "W": C * 11 + H * 10 + N * 2 + O,
    "O": C * 5 + H * 12 + N * 2 + O * 2,  # if used
    # PTM-bearing residue identities (used when listing residues with labels)
    "K(Acetylation)": C * 8 + H * 14 + N * 2 + O * 2,
    "S(Acetylation)": C * 5 + H * 7 + N + O * 3,
    "T(Acetylation)": C * 6 + H * 9 + N + O * 3,
    "K(Biotinylation)": C * 16 + H * 26 + N * 2 + O * 3 + S,
    "C(Carbamidomethylation)": C * 5 + H * 8 + N * 2 + O * 2 + S,
    "D(Carbamidomethylation)": C * 6 + H * 8 + N * 2 + O * 4,
    "E(Carbamidomethylation)": C * 7 + H * 10 + N * 2 + O * 4,
    "H(Carbamidomethylation)": C * 8 + H * 10 + N * 4 + O * 2,
    "K(Carbamidomethylation)": C * 8 + H * 15 + N * 3 + O * 2,
    "D(Carboxylation)": C * 5 + H * 5 + N + O * 5,
    "E(Carboxylation)": C * 6 + H * 7 + N + O * 5,
    "K(Carboxylation)": C * 7 + H * 12 + N * 2 + O * 3,
    "W(Carboxylation)": C * 12 + H * 10 + N * 2 + O * 3,
    "N(Deamidation)": C * 4 + H * 5 + N + O,
    "Q(Deamidation)": C * 5 + H * 7 + N + O,
    "D(Dehydration)": C * 4 + H * 3 + N + O * 2,
    "S(Dehydration)": C * 3 + H * 3 + N + O,
    "T(Dehydration)": C * 4 + H * 5 + N + O,
    "Y(Dehydration)": C * 9 + H * 7 + N + O,
    "K(Methylation)": C * 7 + H * 14 + N * 2 + O,
    "R(Methylation)": C * 7 + H * 14 + N * 4 + O,
    "M(Oxidation)": C * 5 + H * 9 + N + O * 2 + S,
    "D(SodiumAdduct)": C * 4 + H * 4 + N + O * 3 + Na,
    "E(SodiumAdduct)": C * 5 + H * 6 + N + O * 3 + Na,
    "K(12PlexDiLeu)": H * 15 + C * 7 + C13 * 1 + N15 * 1 + O18 * 1,
    "K(mdDiLeu1101)": H * 15 + C * 7 + C13 * 1 + N15 * 1 + O18 * 1,
    "K(mdDiLeu0400)": H * 11 + C * 8 + N * 1 + O * 1 + H2 * 4,
    "E(Glu->pyro-Glu)": C * 5 + H * 4 + O * 3,
    "Q(Gln->pyro-Glu)": C * 5 + H * 5 + N + O * 2,
    "Y(Sulfation)": C * 9 + H * 9 + N + O * 5 + S,
    "S(Phosphorylation)": C * 3 + H * 6 + N + O * 5 + P,
    "T(Phosphorylation)": C * 4 + H * 8 + N + O * 5 + P,
    "Y(Phosphorylation)": C * 9 + H * 10 + N + O * 5 + P,
}

TERMINI = {
    "Standard": H * 2 + O,
    "(Amidated)": N + H * 3,
}
N_TERMINI = {
    "(12PlexDiLeu)": H * 15 + C * 7 + C13 * 1 + N15 * 1 + O18 * 1,
    "(mdDiLeu1101)": H * 15 + C * 7 + C13 * 1 + N15 * 1 + O18 * 1,
    "(mdDiLeu0400)": H * 11 + C * 8 + N * 1 + O * 1 + H2 * 4,
    "(Acetylation)": H * 2 + C * 2 + O,
    "(Carbamidomethylation)": H * 3 + C * 2 + N + O,
}
PTMs = {
    "K(Acetylation)": C * 2 + H * 2 + O,
    "S(Acetylation)": C * 2 + H * 2 + O,
    "T(Acetylation)": C * 2 + H * 2 + O,
    "K(Biotinylation)": H * 14 + C * 10 + N * 2 + O * 2 + S,
    "C(Carbamidomethylation)": C * 2 + H * 3 + N + O,
    "D(Carbamidomethylation)": C * 2 + H * 3 + N + O,
    "E(Carbamidomethylation)": C * 2 + H * 3 + N + O,
    "H(Carbamidomethylation)": C * 2 + H * 3 + N + O,
    "K(Carbamidomethylation)": C * 2 + H * 3 + N + O,
    "D(Carboxylation)": C + O * 2,
    "E(Carboxylation)": C + O * 2,
    "K(Carboxylation)": C + O * 2,
    "W(Carboxylation)": C + O * 2,
    "N(Deamidation)": H * -1 + N * -1 + O,
    "Q(Deamidation)": H * -1 + N * -1 + O,
    "D(Dehydration)": H * -2 + O * -1,
    "S(Dehydration)": H * -2 + O * -1,
    "T(Dehydration)": H * -2 + O * -1,
    "Y(Dehydration)": H * -2 + O * -1,
    "K(Methylation)": C + H * 2,
    "R(Methylation)": C + H * 2,
    "M(Oxidation)": O,
    "D(SodiumAdduct)": H * -1 + Na,
    "E(SodiumAdduct)": H * -1 + Na,
    "K(12PlexDiLeu)": H * 15 + C * 7 + C13 * 1 + N15 * 1 + O18 * 1,
    "K(mdDiLeu1101)": H * 15 + C * 7 + C13 * 1 + N15 * 1 + O18 * 1,
    "K(mdDiLeu0400)": H * 11 + C * 8 + N * 1 + O * 1 + H2 * 4,
    "E(Glu->pyro-Glu)": H * -3 + N * -1,
    "Q(Gln->pyro-Glu)": H * -3 + N * -1,
    "Y(Sulfation)": O * 3 + S,
    "S(Phosphorylation)": P + O * 3 + H,
    "T(Phosphorylation)": P + O * 3 + H,
    "Y(Phosphorylation)": P + O * 3 + H,
}
ADDUCTS = {"H2O": (H * 2) + O, 
           "NH3": N + (H * 3),
           "CO": C + O,
           "H5N2": (H*6) + (N*2),
           "H3N+H2O": (H*5) + N + O,
           "2*H2O": (H*4) + (O*2),
           "HPO3": H + P + (O*3),
           "HPO3+H2O": (O*4) + P + (H*2)}

# def _terminus_mass(token: str) -> float:
#     return TERMINI.get(token, TERMINI["Standard"])


# def _nterm_mass(token: str) -> float:
#     return N_TERMINI.get(token, 0.0)


# def _ptm_mass(token: str) -> float:
#     return PTMs.get(token, 0.0)
def token_mass(tok: str) -> float:
    """
    Return the neutral mass of a single token:
      - N-term label like '(12PlexDiLeu)' comes from N_TERMINI
      - PTM-bearing residue like 'Y(Sulfation)' tries AA_MASSES first,
        then falls back to base + PTMs delta
      - Plain residue like 'Y' from AA_MASSES
    """
    if tok.startswith('('):  # N-term label token, e.g. '(12PlexDiLeu)'
        return N_TERMINI.get(tok, 0.0)

    if '(' in tok:           # PTM-bearing residue, e.g. 'Y(Sulfation)'
        base = tok[0]
        mod  = tok[tok.index('('):]   # '(Sulfation)'
        # try direct lookup first
        if tok in AA_MASSES:
            return AA_MASSES[tok]
        # else compose base + delta
        return AA_MASSES[base] + PTMs.get(f"{base}{mod}", 0.0)
    # plain residue
    return AA_MASSES[tok]

def check_termini_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return TERMINI['Standard']

def check_PTM_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return 0

def check_termini(pot_mod_peptide):
    if '(' in pot_mod_peptide:
        term_start = (pot_mod_peptide.rindex('('))
        termini_ID = pot_mod_peptide[term_start:]
        termini_mass_change = check_termini_Key(TERMINI, termini_ID)
        return termini_mass_change
    else:
        return TERMINI['Standard']

def check_PTM(pot_mod_peptide):
    number_of_mods = pot_mod_peptide.count('(')
    if number_of_mods > 0:
        current_peptide = []
        mass_change_collection = []
        current_peptide.append(pot_mod_peptide)
        for a in range(0,number_of_mods):
            peptide_less_mod = current_peptide[-1]
            ptm_start = (peptide_less_mod.index('('))-1
            ptm_end = (peptide_less_mod.index(')'))+1
            ptm_ID = peptide_less_mod[ptm_start:ptm_end]
            ptm_mass_change = check_PTM_Key(PTMs, ptm_ID)
            mass_change_collection.append(ptm_mass_change)
            peptide_less_mod2 = peptide_less_mod[:ptm_start] + peptide_less_mod[ptm_end:]
            current_peptide.append(peptide_less_mod2)
            
        ptm_mass_change = sum(mass_change_collection)
        return ptm_mass_change
    else:
        ptm_mass_change = 0
        return ptm_mass_change

# --------------------------- #
# Residue & mass computations #
# --------------------------- #
# def list_of_residues(peptide: str) -> list[str]:
#     """
#     Returns a list of residues where PTM-bearing residues are kept together, e.g. "K(Acetylation)".
#     Handles possible N-term tokens like "(Acetylation)PEPTIDE".
#     """
#     res = []
#     s = peptide

#     # handle leading N-terminal token "(...)" if present
#     if s.startswith("("):
#         end = s.index(")")
#         res.append(s[: end + 1])  # e.g., "(Acetylation)"
#         s = s[end + 1 :]  # remainder after N-term token

#     i = 0
#     while i < len(s):
#         if i + 1 < len(s) and s[i + 1] == "(":
#             j = s.index(")", i + 1)
#             res.append(s[i : j + 1])  # e.g., "K(Acetylation)"
#             i = j + 1
#         else:
#             res.append(s[i])
#             i += 1
#     return res

def list_of_residues(pot_mod_peptide):
     list_of_res = []
     pep_update = []
     pep_update.append(pot_mod_peptide)
     no_mods = pot_mod_peptide.count('(')
     if no_mods > 1:
         for c in range(0,no_mods+1):
             pep_of_interest = pep_update[-1]
             if '(' in pep_of_interest:
                 if pep_of_interest.startswith('('):
                     c_term_last_res = peptide.index(')')
                     c_term_mod = peptide[:c_term_last_res+1]
                     list_of_res.append(c_term_mod)
                     
                     remaining_pep = peptide[c_term_last_res+1:]
                     
                     
                     first_ptm_start = remaining_pep.index('(')
                     first_ptm_end = remaining_pep.index(')')
                     first_residues = remaining_pep[:(first_ptm_start-1)]
                     for a in first_residues:
                         list_of_res.append(a)
                     ptm_residue = remaining_pep[(first_ptm_start-1):(first_ptm_end+1)]
                     list_of_res.append(ptm_residue)
                     remaining_pep = remaining_pep[(first_ptm_end+1):]
                     pep_update.append(remaining_pep)  
                 else:
                     first_ptm_start = pep_of_interest.index('(')
                     first_ptm_end = pep_of_interest.index(')')
                     first_residues = pep_of_interest[:(first_ptm_start-1)]
                     for a in first_residues:
                         list_of_res.append(a)
                     ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                     list_of_res.append(ptm_residue)
                     remaining_pep = pep_of_interest[(first_ptm_end+1):]
                     pep_update.append(remaining_pep)  
             else:
                 for d in pep_of_interest:
                     list_of_res.append(d)
     elif no_mods == 1:
         for c in range(0,1):
             pep_of_interest = pep_update[-1]
             if '(' in pep_of_interest:
                 first_ptm_start = pep_of_interest.index('(')
                 first_ptm_end = pep_of_interest.index(')')
                 if first_ptm_start == 1:
                     ptm_residue =  pep_of_interest[0] + (pep_of_interest[(first_ptm_start):(first_ptm_end+1)])
                     list_of_res.append(ptm_residue)
                     remaining_pep = pep_of_interest[(first_ptm_end+1):]
                     for d in remaining_pep:
                         list_of_res.append(d)
                 if first_ptm_start != 1:
                     first_residues = pep_of_interest[:(first_ptm_start-1)]
                     for a in first_residues:
                         list_of_res.append(a)
                     ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                     list_of_res.append(ptm_residue)
                     remaining_pep = pep_of_interest[(first_ptm_end+1):]
                     for d in remaining_pep:
                         list_of_res.append(d)              
             else:
                 for d in pep_of_interest:
                     list_of_res.append(d) 
     elif no_mods == 0:
         for c in pot_mod_peptide:
             list_of_res.append(c)
     return list_of_res

# def check_termini(peptide: str) -> float:
#     """
#     Returns terminal mass contribution. If peptide has leading "(...)", interpret as an N-term modification.
#     Otherwise returns standard termini mass.
#     """
#     if peptide.startswith("("):
#         end = peptide.index(")")
#         token = peptide[: end + 1]
#         return _nterm_mass(token)
#     return TERMINI["Standard"]

# def check_PTM(peptide: str) -> float:
#     """
#     Sum of PTM mass changes on residues, e.g., "K(Acetylation)".
#     """
#     total = 0.0
#     for token in re.findall(r"[A-Z]\(.*?\)", peptide):
#         total += _ptm_mass(token)
#     return total

# ----------------------------- #
# Theoretical spectra generator #
# ----------------------------- #
# def theoretical_spectra_generator(peptide: str, charge: int = 1) -> pd.DataFrame:
#     """
#     Build b/y ions (and -H2O / -NH3 adducts) with monoisotopic masses.
#     Returns DataFrame with columns: ['ion', 'Fragment theoretical monoisotopic mass']
#     """
#     # plain sequence for counting
#     plain = re.sub(r"[\(\[].*?[\)\]]", "", peptide)
#     res_list = list_of_residues(peptide)

#     # total mass (MH)
#     mass_of_residues = [AA_MASSES[r] for r in plain]
#     peptide_mass = sum(mass_of_residues) + check_termini(peptide) + check_PTM(peptide)
#     mz_MH = (peptide_mass + (proton_mass * charge)) / charge

#     num_ions = max(len(plain) - 1, 0)

#     # b-ions
#     b_masses, b_names = [], []
#     running = 0.0
#     # account for possible first token being N-term "(...)"
#     idx = 0
#     if res_list and res_list[0].startswith("("):
#         running += _nterm_mass(res_list[0])
#         idx = 1
#     for k in range(idx, idx + num_ions):
#         rk = res_list[k]
#         running += AA_MASSES[rk] if rk in AA_MASSES else AA_MASSES.get(rk, 0.0)
#         b_masses.append(running + proton_mass)
#         b_names.append(f"b{k - idx + 1}")

#     # y-ions (using running subtraction from full MH)
#     y_masses, y_names = [], []
#     running_y = mz_MH
#     for k in range(len(res_list) - 1, len(res_list) - 1 - num_ions, -1):
#         rk = res_list[k]
#         running_y -= AA_MASSES[rk] if rk in AA_MASSES else AA_MASSES.get(rk, 0.0)
#         y_masses.append(running_y)
#         y_names.append(f"y{len(y_masses)}")
#     # fix y names to count from C-term: yN,...,y1
#     if y_names:
#         total = len(y_names)
#         y_names = [f"y{total - i}" for i in range(total)]

#     # assemble reports
#     b_df = pd.DataFrame({"ion": b_names, "mass": b_masses})
#     y_df = pd.DataFrame({"ion": y_names, "mass": y_masses})

#     b_h2o = pd.DataFrame({"ion": b_df["ion"] + "-H2O", "mass": b_df["mass"] - ADDUCTS["H2O"]})
#     b_nh3 = pd.DataFrame({"ion": b_df["ion"] + "-NH3", "mass": b_df["mass"] - ADDUCTS["NH3"]})

#     y_h2o = pd.DataFrame({"ion": y_df["ion"] + "-H2O", "mass": y_df["mass"] - ADDUCTS["H2O"]})
#     y_nh3 = pd.DataFrame({"ion": y_df["ion"] + "-NH3", "mass": y_df["mass"] - ADDUCTS["NH3"]})

#     ion_report = pd.concat([b_df, y_df, b_h2o, b_nh3, y_nh3, y_h2o], ignore_index=True).drop_duplicates()
#     ion_report = ion_report.rename(columns={"mass": "Fragment theoretical monoisotopic mass"})
#     return ion_report
def tokenize_residues(pep: str) -> list[str]:
    """
    Convert a peptide string into tokens:
      ['(12PlexDiLeu)', 'Y', 'A', 'F', 'G', 'L(Amidated)']  (example)
    Uses your existing list_of_residues(), but preserves an N-term label
    as a standalone token at the front.
    """
    if pep.startswith('('):
        end = pep.index(')') + 1
        nterm = pep[:end]
        rest  = pep[end:]
        toks = list_of_residues(rest)
        return [nterm] + toks
    else:
        return list_of_residues(pep)

def cterm_mass_for_y(pep_with_mods: str) -> float:
    """
    y-ions include the C-terminus. Use your TERMINI mapping:
      - 'Standard' means +H2O
      - '(Amidated)' means +NH3
    We detect C-term choice from the peptide tail, reusing your check_termini logic.
    """
    # check_termini returns the *peptide termini* mass; for y-ions we only need
    # the C-term contribution. Your TERMINI dict already encodes 'H2O' as Standard
    # and 'NH3' as (Amidated), which is exactly what y-ions need.
    return check_termini(pep_with_mods)

def theoretical_spectra_generator(peptide_to_check: str, charge: int = 1) -> pd.DataFrame:
    """
    Build b/y (and -H2O / -NH3) using proper neutral-mass accumulation.
    """
    # tokens containing N-term labels and PTM-bearing residues
    toks = tokenize_residues(peptide_to_check)

    # figure out where residues start if there is an N-term label
    start_idx = 1 if toks and toks[0].startswith('(') else 0

    # residue-only slice (no N-term label)
    residues = toks[start_idx:]

    # full peptide neutral mass (for info; not required to build y-ions anymore)
    full_neutral = sum(token_mass(t) for t in toks) + (TERMINI['Standard'] - cterm_mass_for_y(peptide_to_check))
    # ^ explanation:
    # token_mass already includes any N-term label in toks.
    # We need to ensure we only add *one* C-term choice. list_of_residues returns the
    # last residue possibly as 'X(Amidated)'; the actual C-term contribution for the whole
    # peptide is handled by check_termini(). So normalize by adding Standard and then
    # replacing with the actual cterm (difference term). This keeps full_neutral coherent.

    # number of internal cleavages
    n = len(residues) - 1
    if n <= 0:
        return pd.DataFrame(columns=['ion', 'Fragment theoretical monoisotopic mass'])

    # --- b-ions ---
    b_neutral = []
    running = 0.0

    # include N-term label mass once at the beginning of b-series if present
    if start_idx == 1:
        running += token_mass(toks[0])

    for i in range(n):
        running += token_mass(residues[i])  # add residue i neutral mass
        b_neutral.append(running)

    b_mz = [ (m + proton_mass*charge) / charge for m in b_neutral ]
    b_names = [ f"b{i+1}" for i in range(n) ]

    # --- y-ions ---
    # y_n is the C-terminal n-residue fragment + C-terminus group
    cterm = cterm_mass_for_y(peptide_to_check)

    y_neutral = []
    running = cterm
    for k in range(n):
        # walk from C-terminus backwards
        idx = len(residues) - 1 - k
        running += token_mass(residues[idx])
        y_neutral.append(running)

    y_mz = [ (m + proton_mass*charge) / charge for m in y_neutral ]
    y_names = [ f"y{(k+1)}" for k in range(n) ]  # y1 is the most C-terminal

    # package
    b_df = pd.DataFrame({'ion': b_names, 'Fragment theoretical monoisotopic mass': b_mz})
    y_df = pd.DataFrame({'ion': y_names, 'Fragment theoretical monoisotopic mass': y_mz})

    # neutral losses (apply to m/z values)
    def add_loss(df, label, loss_neutral):
        out = df.copy()
        out['ion'] = out['ion'] + label
        out['Fragment theoretical monoisotopic mass'] = out['Fragment theoretical monoisotopic mass'] - (loss_neutral / charge)
        return out

    b_h2o = add_loss(b_df, '-H2O', ADDUCTS['H2O'])
    b_nh3 = add_loss(b_df, '-NH3', ADDUCTS['NH3'])
    y_h2o = add_loss(y_df, '-H2O', ADDUCTS['H2O'])
    y_nh3 = add_loss(y_df, '-NH3', ADDUCTS['NH3'])
    # y_co = add_loss(y_df, '-CO', ADDUCTS['CO'])
    # b_co = add_loss(b_df, '-CO', ADDUCTS['CO'])
    # y_h5n2 = add_loss(y_df, '-H5N2', ADDUCTS['H5N2'])
    # b_h5n2 = add_loss(b_df, '-H5N2', ADDUCTS['H5N2'])
    # y_h3n_h2o = add_loss(y_df, '-H3N+H2O', ADDUCTS['H3N+H2O'])
    # b_h3n_h2o = add_loss(b_df, '-H3N+H2O', ADDUCTS['H3N+H2O'])
    # y_2_h2o = add_loss(y_df, '-2*H2O', ADDUCTS['2*H2O'])
    # b_2_h2o = add_loss(b_df, '-2*H2O', ADDUCTS['2*H2O'])
    # y_hpo3 = add_loss(y_df, '-HPO3', ADDUCTS['HPO3'])
    # b_hpo3 = add_loss(b_df, '-HPO3', ADDUCTS['HPO3'])
    # y_hpo3_h2o = add_loss(y_df, '-HPO3+H2O', ADDUCTS['HPO3+H2O'])
    # b_hpo3_h2o = add_loss(b_df, '-HPO3+H2O', ADDUCTS['HPO3+H2O'])
    
    ion_report = pd.concat([b_df, y_df, 
                            b_h2o, y_h2o, 
                            b_nh3, y_nh3,
                            # b_co,y_co,
                            # b_h5n2,y_h5n2,
                            # b_h3n_h2o,y_h3n_h2o,
                            # b_2_h2o,y_2_h2o,
                            # b_hpo3,y_hpo3,
                            # b_hpo3_h2o,y_hpo3_h2o,
                            ], ignore_index=True).drop_duplicates()
    return ion_report


# ---------------- #
# Helper: IO parse #
# ---------------- #
def _load_raw_converted_csv(raw_converter_path: str) -> pd.DataFrame:
    """
    Loads the formatted raw CSV (from your converter) and returns a DataFrame with standardized columns.
    """
    df = pd.read_csv(
        raw_converter_path,
        sep=",",
        skiprows=[0],
        names=["fragment_mz",
            "fragment_intensity",
            "fragment_z",
            "fragment_resolution",
            "precursor_mz",
            "ms2_scan",
            "precursor_z",
            "precursor_RT",
            "IonInjectTime",
            "ms1_scan",
            "precursor_intensity",
            "null"])

    df = df.rename(columns={
            "fragment_mz": "Fragment actual m/z",
            "fragment_z": "Fragment actual charge",
            "fragment_intensity": "Fragment actual intensity",
            "precursor_mz": "Precursor actual m/z",
            "precursor_z": "Precursor actual charge",
            "ms2_scan": "Scan"})
    df["Fragment actual charge"] = df["Fragment actual charge"].replace(0, 1)
    df = df.copy()
    df["Precursor actual monoisotopic mass"] = (
        df["Precursor actual m/z"] * df["Precursor actual charge"] - h_mass * df["Precursor actual charge"])
    return df

# ------------------------- #
# Sequence coverage utility #
# ------------------------- #
def _seq_coverage_calc(frag_df: pd.DataFrame, ion_report: pd.DataFrame, scan: int, peptide: str, frag_da_cutoff: float):
    matched = frag_df.loc[frag_df["Fragment error (Da)"] <= frag_da_cutoff].copy()
    if matched.empty:
        return pd.DataFrame({"Sequence": [peptide], "Scan": [scan], "Sequence coverage": [0.0]})

    matched["ion_format"] = matched["ion"].str.extract(r"(\d+)", expand=False)
    ion_tmp = ion_report[["ion"]].copy()
    ion_tmp["ion_format"] = ion_tmp["ion"].str.extract(r"(\d+)", expand=False)

    seq_cov = (
        matched.drop_duplicates(subset="ion_format").shape[0]
        / max(ion_tmp.drop_duplicates(subset="ion_format").shape[0], 1)
        * 100.0)
    return pd.DataFrame({"Sequence": [peptide], "Scan": [scan], "Sequence coverage": [seq_cov]})

# ------------------------------------- #
# Hyperscore (OpenMS) with token map    #
# ------------------------------------- #
def _build_peptide_for_openms(peptide: str) -> str:
    """
    Converts tokens like 'K(Acetylation)' and leading '(Acetylation)' to OpenMS AASequence format with bracketed masses.
    """
    token_to_bracket = {"E(Glu->pyro-Glu)": f"E[{PTMs['E(Glu->pyro-Glu)']}]",
        "K(12PlexDiLeu)": f"K[{PTMs['K(12PlexDiLeu)']}]",
        "K(mdDiLeu1101)": f"K[{PTMs['K(mdDiLeu1101)']}]",
        "K(mdDiLeu0400)": f"K[{PTMs['K(mdDiLeu0400)']}]",
        "Q(Gln->pyro-Glu)": f"Q[{PTMs['Q(Gln->pyro-Glu)']}]",
        "K(Acetylation)": f"K[{PTMs['K(Acetylation)']}]",
        "S(Acetylation)": f"S[{PTMs['S(Acetylation)']}]",
        "T(Acetylation)": f"T[{PTMs['T(Acetylation)']}]",
        "K(Biotinylation)": f"K[{PTMs['K(Biotinylation)']}]",
        "C(Carbamidomethylation)": f"C[{PTMs['C(Carbamidomethylation)']}]",
        "D(Carbamidomethylation)": f"D[{PTMs['D(Carbamidomethylation)']}]",
        "E(Carbamidomethylation)": f"E[{PTMs['E(Carbamidomethylation)']}]",
        "H(Carbamidomethylation)": f"H[{PTMs['H(Carbamidomethylation)']}]",
        "K(Carbamidomethylation)": f"K[{PTMs['K(Carbamidomethylation)']}]",
        "D(Carboxylation)": f"D[{PTMs['D(Carboxylation)']}]",
        "E(Carboxylation)": f"E[{PTMs['E(Carboxylation)']}]",
        "K(Carboxylation)": f"K[{PTMs['K(Carboxylation)']}]",
        "W(Carboxylation)": f"W[{PTMs['W(Carboxylation)']}]",
        "N(Deamidation)": f"N[{PTMs['N(Deamidation)']}]",
        "Q(Deamidation)": f"Q[{PTMs['Q(Deamidation)']}]",
        "D(Dehydration)": f"D[{PTMs['D(Dehydration)']}]",
        "S(Dehydration)": f"S[{PTMs['S(Dehydration)']}]",
        "T(Dehydration)": f"T[{PTMs['T(Dehydration)']}]",
        "Y(Dehydration)": f"Y[{PTMs['Y(Dehydration)']}]",
        "K(Methylation)": f"K[{PTMs['K(Methylation)']}]",
        "R(Methylation)": f"R[{PTMs['R(Methylation)']}]",
        "M(Oxidation)": f"M[{PTMs['M(Oxidation)']}]",
        "D(SodiumAdduct)": f"D[{PTMs['D(SodiumAdduct)']}]",
        "E(SodiumAdduct)": f"E[{PTMs['E(SodiumAdduct)']}]",
        "Y(Sulfation)": f"Y[{PTMs['Y(Sulfation)']}]",
        "Y(Phosphorylation)": f"Y[{PTMs['Y(Phosphorylation)']}]",
        "S(Phosphorylation)": f"S[{PTMs['S(Phosphorylation)']}]",
        "T(Phosphorylation)": f"T[{PTMs['T(Phosphorylation)']}]"}
    nterm_to_bracket = {
        "(Carbamidomethylation)": f"[{N_TERMINI['(Carbamidomethylation)']}]",
        "(12PlexDiLeu)": f"[{N_TERMINI['(12PlexDiLeu)']}]",
        "(mdDiLeu1101)": f"[{N_TERMINI['(mdDiLeu1101)']}]",
        "(mdDiLeu0400)": f"[{N_TERMINI['(mdDiLeu0400)']}]",
        "(Acetylation)": f"[{N_TERMINI['(Acetylation)']}]",
    }

    s = peptide
    for old, new in token_to_bracket.items():
        s = s.replace(old, new)
    for old, new in nterm_to_bracket.items():
        s = s.replace(old, new)
    return s


def _hyperscore_calc(openms_exp: MSExperiment, scan_number: int, peptide: str, min_intensity: float, min_mz: float,
                     precursor_ppm_tolerance: float) -> float:
    """
    Compute OpenMS HyperScore for a (scan, peptide).
    """
    tsg = TheoreticalSpectrumGenerator()
    thspec = MSSpectrum()
    p = Param()
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)

    pep_str = _build_peptide_for_openms(peptide)
    pep = AASequence.fromString(pep_str)
    tsg.getSpectrum(thspec, pep, 1, 1)

    spectrum_of_interest = openms_exp[scan_number - 1]
    mz, inten = spectrum_of_interest.get_peaks()
    # you filtered peaks later; HyperScore uses the full spectrum, tol set separately
    hscore = HyperScore()
    result = hscore.compute(precursor_ppm_tolerance, True, spectrum_of_interest, thspec)
    return float(result)
# ----------------------------------------------- #
# Main entry point
# ----------------------------------------------- #

def launch_db_search_pt1(predefined_db_path,output_parent_directory,choose_mzml_directory,raw_file_formatted_path,precursor_error_cutoff,
                         fragment_error_cutoff,min_mz,min_intensity_pass,standard_err_percent,amid_status_var_status,acet_term_status_var_status,
                         acet_k_status_var_status,acet_s_status_var_status,acet_t_status_var_status,biotin_k_status_var_status,carb_term_status_var_status,
                         carb_c_status_var_status,carb_d_status_var_status,carb_e_status_var_status,carb_h_status_var_status,carb_k_status_var_status,
                         carboxy_d_status_var_status,carboxy_e_status_var_status,carboxy_k_status_var_status,carboxy_w_status_var_status,deamid_n_status_var_status,
                         deamid_q_status_var_status,dehyd_d_status_var_status,dehyd_s_status_var_status,dehyd_t_status_var_status,dehyd_y_status_var_status,
                         meth_k_status_var_status,meth_r_status_var_status,ox_m_status_var_status,sodi_d_status_var_status,sodi_e_status_var_status,
                         plex12_nterm_status_var_status,plex12_k_status_var_status,pg_E_status_var_status,pg_q_status_var_status,sulfo_y_status_var_status,
                         phospho_s_status_var_status,phospho_t_status_var_status,phospho_y_status_var_status,mddileu_1101_nterm_status_var_status,
                         mddileu_1101_k_status_var_status,mddileu_0400_nterm_status_var_status,mddileu_0400_k_status_var_status,max_modifications,
                         sample_output_directory, fragment_store="parquet"):

    if "_formatted.ms2" in choose_mzml_directory:
        mzml_path = choose_mzml_directory.replace("_formatted.ms2", ".mzML")
    elif "_formatted.txt" in choose_mzml_directory:
        mzml_path = choose_mzml_directory.replace("_formatted.txt", ".mzML")
    elif choose_mzml_directory.endswith(".txt"):
        mzml_path = choose_mzml_directory.replace(".txt", ".mzML")
    elif choose_mzml_directory.lower().endswith(".ms2"):
        mzml_path = choose_mzml_directory[:-4] + ".mzML"
    else:
        mzml_path = choose_mzml_directory  # last resort

    raw_conv_mzml_storage = [[mzml_path,raw_file_formatted_path]]

 
    fasta_w_mass = pd.read_csv(predefined_db_path)
    if fasta_w_mass.empty:
        raise ValueError("Database file is empty")
    
    # Precursor loose tolerance for asof in Daltons (rough ppm→Da for windowing)
    precursor_temp_cutoff_da = precursor_error_cutoff * 3

    all_seq_cov_rows = []

    for mzml_path_input, raw_converter_path in raw_conv_mzml_storage:
       
        # Output location details
        sample_name, peptide_report_output = raw_file_detail_extraction(raw_converter_path, output_parent_directory)

        # Read experimental spectra
        exp_precursor = _load_raw_converted_csv(raw_converter_path)
        exp_precursor = exp_precursor.drop(
            columns=["fragment_resolution", "precursor_RT", "IonInjectTime", "ms1_scan", "precursor_intensity", "null"]
        )

        exp_mass_only = exp_precursor.drop(
            columns=["Fragment actual m/z", "Fragment actual charge", "Fragment actual intensity", "Precursor actual m/z", "Precursor actual charge"]
        )

        db_mass_only = fasta_w_mass.copy()
    
        # Sorted unique by mass (merge_asof requirement)
        exp_sorted = (
            exp_mass_only.drop_duplicates(subset=["Precursor actual monoisotopic mass", "Scan"])
            .sort_values(by="Precursor actual monoisotopic mass")
        )
        db_sorted = db_mass_only.drop_duplicates(subset="Precursor theoretical monoisotopic mass").sort_values(
            by="Precursor theoretical monoisotopic mass"
        )

        exp_sorted_unique_mass = exp_sorted.drop_duplicates(subset="Precursor actual monoisotopic mass")

        mm_forward = pd.merge_asof(
            exp_sorted_unique_mass,
            db_sorted,
            left_on="Precursor actual monoisotopic mass",
            right_on="Precursor theoretical monoisotopic mass",
            tolerance=precursor_temp_cutoff_da,
            allow_exact_matches=True,
            direction="forward",
        )
        mm_backward = pd.merge_asof(
            exp_sorted_unique_mass,
            db_sorted,
            left_on="Precursor actual monoisotopic mass",
            right_on="Precursor theoretical monoisotopic mass",
            tolerance=precursor_temp_cutoff_da,
            allow_exact_matches=True,
            direction="backward",
        )
        merge_match = pd.concat([mm_forward, mm_backward], ignore_index=True)
        merge_match = merge_match[merge_match["Sequence"].notna()].drop_duplicates().sort_values(by="Scan")

        merge_again = pd.merge(merge_match, db_mass_only, on="Precursor theoretical monoisotopic mass", how="inner")
        merge_again2 = pd.merge(merge_again, exp_sorted, on="Precursor actual monoisotopic mass", how="inner")
        merge_match_count = merge_again2.drop_duplicates(subset=["Scan_y", "Sequence_y"]).copy()
        merge_match_count = merge_match_count.drop(columns=["Sequence_x", "Scan_x"]).rename(
            columns={"Sequence_y": "Sequence", "Scan_y": "Scan"}
        )

        # ppm error
        merge_match_count["Precursor error (ppm)"] = (
            merge_match_count["Precursor theoretical monoisotopic mass"]
            .sub(merge_match_count["Precursor actual monoisotopic mass"])
            .abs()
            .div(merge_match_count["Precursor theoretical monoisotopic mass"])
            * 1e6
        )

        precursor_amm_results = merge_match_count[merge_match_count["Precursor error (ppm)"] <= precursor_error_cutoff].copy()
        
        # Save precursor AMM results
        Path(sample_output_directory).mkdir(parents=True, exist_ok=True)
        (Path(sample_output_directory) / "fragment_matches").mkdir(parents=True, exist_ok=True)
        precursor_csv = Path(sample_output_directory) / "precursor_AMM_results.csv"
        precursor_amm_results.to_csv(precursor_csv, index=False)
                
        # Prepare for fragment matching
        all_seqs = (
            precursor_amm_results.drop(columns=["Scan", "Precursor theoretical monoisotopic mass", "Precursor actual monoisotopic mass"])
            .drop_duplicates(subset="Sequence", ignore_index=True)
        )

        # Fragment matching + sequence coverage per candidate (scan, peptide)
        unfiltered_psms = []
        fragment_buffers = []
        for idx_row in range(len(precursor_amm_results)):
            row = precursor_amm_results.iloc[[idx_row]]
            candidate_seq = row["Sequence"].values[0]
            candidate_scan = row["Scan"].values[0]

            # Theoretical ion series for the peptide
            ion_report = theoretical_spectra_generator(candidate_seq)

            # Join to get the specific (scan, candidate_seq)-filtered experimental fragments
            scan_extract = pd.merge(all_seqs[all_seqs["Sequence"] == candidate_seq], precursor_amm_results, on="Sequence")
            fragments_extract = pd.merge(scan_extract, exp_precursor, on=["Scan"])

            # Compute fragment monoisotopic mass and filter by current scan
            fragments_extract["Fragment actual monoisotopic mass"] = (
                fragments_extract["Fragment actual m/z"] * fragments_extract["Fragment actual charge"]
                - h_mass * fragments_extract["Fragment actual charge"]
            )
            fragments_extract = fragments_extract[fragments_extract["Scan"] == candidate_scan].copy()

            # Vectorized asof for fragment matching (two directions)
            frag_tol_da = fragment_error_cutoff * 1_000_000  # original line used *1e6 as "Da"; keep identical
            ion_sorted = ion_report.sort_values("Fragment theoretical monoisotopic mass")
            frag_sorted = fragments_extract.sort_values("Fragment actual monoisotopic mass")

            frag_fwd = pd.merge_asof(
                ion_sorted,
                frag_sorted,
                left_on="Fragment theoretical monoisotopic mass",
                right_on="Fragment actual monoisotopic mass",
                tolerance=frag_tol_da,
                allow_exact_matches=True,
                direction="forward",
            )
            frag_bwd = pd.merge_asof(
                ion_sorted,
                frag_sorted,
                left_on="Fragment theoretical monoisotopic mass",
                right_on="Fragment actual monoisotopic mass",
                tolerance=frag_tol_da,
                allow_exact_matches=True,
                direction="backward",
            )
            frag_merge_match = pd.concat([frag_fwd, frag_bwd], ignore_index=True)
            frag_merge_match = frag_merge_match[frag_merge_match["Sequence"].notna()].copy()
            frag_merge_match["Fragment error (Da)"] = (
                frag_merge_match["Fragment actual monoisotopic mass"]
                - frag_merge_match["Fragment theoretical monoisotopic mass"]
            ).abs()

            # Sequence coverage row
            seq_cov_df = _seq_coverage_calc(
                frag_merge_match, ion_report, int(candidate_scan), candidate_seq, fragment_error_cutoff
            )
            unfiltered_psms.append(seq_cov_df)

            # Save per-peptide/scan fragment report (buffered or legacy-per-file)
            scan_to_report = candidate_scan
            peptide = candidate_seq

            frag_df = frag_merge_match.copy()
            frag_df["Scan"] = scan_to_report
            frag_df["Sequence"] = peptide

            if fragment_store == "none":
                # no disk writes
                pass
            elif fragment_store in ("parquet", "csv_gz"):
                # keep useful columns to reduce memory/disk
                keep_cols = [
                    "Scan", "Sequence", "ion",
                    "Fragment theoretical monoisotopic mass",
                    "Fragment actual monoisotopic mass",
                    "Fragment error (Da)",
                    "Fragment actual intensity",
                ]
                cols = [c for c in keep_cols if c in frag_df.columns]
                fragment_buffers.append(frag_df[cols])
            else:
                # legacy per-file CSV write
                peptide_formal = (
                    peptide.replace("Q(Gln->pyro-Glu)", "Q(pyroGlu)")
                           .replace("E(Glu->pyro-Glu)", "E(pyroGlu)")
                )
                out_dir_frag = Path(sample_output_directory) / "fragment_matches"
                out_dir_frag.mkdir(parents=True, exist_ok=True)
                output_path_rep = out_dir_frag / f"{peptide_formal}_{scan_to_report}_fragment_report.csv"
                frag_df.to_csv(output_path_rep, index=False)

        # ------------------------- #
        # HyperScore via OpenMS     #
        # ------------------------- #
        e = MSExperiment()
        MzMLFile().load(mzml_path_input, e)

        unfiltered_psms_df = pd.concat(unfiltered_psms, ignore_index=True) if unfiltered_psms else pd.DataFrame(
            columns=["Sequence", "Scan", "Sequence coverage"]
        )
        unique_scans = unfiltered_psms_df.drop_duplicates(subset="Scan")

        corr_rows = []
        for i in range(len(unique_scans)):
            scan_val = int(unique_scans.iloc[i]["Scan"])
            per_scan = unfiltered_psms_df[unfiltered_psms_df["Scan"] == scan_val]

            for j in range(len(per_scan)):
                peptide_val = per_scan.iloc[j]["Sequence"]
                corr_val = _hyperscore_calc(
                    e,
                    scan_val,
                    peptide_val,
                    min_intensity_pass,
                    min_mz,
                    precursor_error_cutoff,  # HyperScore uses this as ppm tol
                )
                row = per_scan.iloc[[j]].copy()
                row["Correlation value"] = corr_val
                corr_rows.append(row)

        all_correlation_results = (
            pd.concat(corr_rows, ignore_index=True) if corr_rows else pd.DataFrame(columns=["Sequence", "Scan", "Sequence coverage", "Correlation value"])
        )
        out_corr = Path(sample_output_directory) / "all_correlation_results.csv"
        all_correlation_results.to_csv(out_corr, index=False)
        # ---- Consolidated fragment write ----
        if fragment_store in ("parquet", "csv_gz") and fragment_buffers:
            all_frags = pd.concat(fragment_buffers, ignore_index=True)
            out_dir_frag = Path(sample_output_directory)
            if fragment_store == "parquet":
                try:
                    all_frags.to_parquet(out_dir_frag / "fragment_matches.parquet", index=False)
                except Exception:
                    # Fallback if parquet engine not available
                    all_frags.to_csv(out_dir_frag / "fragment_matches.csv.gz", index=False, compression="gzip")
            elif fragment_store == "csv_gz":
                all_frags.to_csv(out_dir_frag / "fragment_matches.csv.gz", index=False, compression="gzip")

        return all_correlation_results

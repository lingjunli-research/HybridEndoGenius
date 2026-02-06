# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:01:10 2023

@author: lawashburn
"""

from Bio.SeqIO.FastaIO import SimpleFastaParser
from pyteomics import parser
import os
import re
import random
import pandas as pd
from utils import *


def make_a_DB(variable_mod_dict,fasta_path,output_folder,max_mods_number):
    """
    Build a target+decoy peptide database with variable modifications applied.
    Writes <Shuffle>_database.csv in output_folder and returns its path.

    Notes:
      - Decoy algorithm is 'Shuffle' (matching your original).
      - m/z is computed for charge 1 (i.e., neutral mass + proton).
      - Legacy tokens (e.g., 'Acetylation-' or '(Acetylation)K') are normalized.
    """
    print("Variable mod dict:", variable_mod_dict)

    decoy_algorithm = 'Shuffle'

    # Atomic masses
    H  = 1.0078250352
    O  = 15.99491463
    C  = 12.0000000
    N  = 14.003074
    P  = 30.973762
    S  = 31.9720707
    C13 = 13.003355
    O18 = 17.9991610
    N15 = 15.0001
    H2  = 2.014101777844
    Na  = 22.989769
    
    proton_mass = 1.00727646688
    charge = 1  # fragment/precursor charge used for reported m/z

    # Base AA (neutral) masses (monoisotopic, residues in peptide)
    aa_masses = {
        "G": C*2  + H*3  + N   + O,
        "A": C*3  + H*5  + N   + O,
        "S": C*3  + H*5  + N   + O*2,
        "P": C*5  + H*7  + N   + O,
        "V": C*5  + H*9  + N   + O,
        "T": C*4  + H*7  + N   + O*2,
        "C": C*3  + H*5  + N   + O   + S,
        "L": C*6  + H*11 + N   + O,
        "I": C*6  + H*11 + N   + O,
        "N": C*4  + H*6  + N*2 + O*2,
        "D": C*4  + H*5  + N   + O*3,
        "Q": C*5  + H*8  + N*2 + O*2,
        "K": C*6  + H*12 + N*2 + O,
        "E": C*5  + H*7  + N   + O*3,
        "M": C*5  + H*9  + N   + O   + S,
        "H": C*6  + H*7  + N*3 + O,
        "F": C*9  + H*9  + N   + O,
        "R": C*6  + H*12 + N*4 + O,
        "Y": C*9  + H*9  + N   + O*2,
        "W": C*11 + H*10 + N*2 + O,
        "O": C*5  + H*12 + N*2 + O*2,  # Pyrrolysine
        "X": 0.0,                       # unknown
    }
    
    # Terminal modifications (mass deltas)
    n_term_mods = {
        "Standard": 0.0,
        "(12PlexDiLeu)": H*15 + C*7 + C13 + N15 + O18,
        "(mdDiLeu1101)": H*15 + C*7 + C13 + N15 + O18,
        "(mdDiLeu0040)": H*11 + C*8 + N*1 + O*1 + H2*4,
        "(Acetylation)": H*2 + C*2 + O,
        "(Carbamidomethylation)": H*3 + C*2 + N + O,
    }

    c_term_mods = {
        "Standard": O + H*2,  # –COOH
        "(Amidated)": N + H*3,  # –CONH2
    }
    
    # Residue-attached PTMs (mass deltas)
    PTMs = {
        "K(Acetylation)"        : C*2 + H*2 + O,
        "S(Acetylation)"        : C*2 + H*2 + O,
        "T(Acetylation)"        : C*2 + H*2 + O,
        "K(Biotinylation)"      : H*14 + C*10 + N*2 + O*2 + S,
        "C(Carbamidomethylation)": C*2 + H*3 + N + O,
        "D(Carbamidomethylation)": C*2 + H*3 + N + O,
        "E(Carbamidomethylation)": C*2 + H*3 + N + O,
        "H(Carbamidomethylation)": C*2 + H*3 + N + O,
        "K(Carbamidomethylation)": C*2 + H*3 + N + O,
        "D(Carboxylation)"      : C + O*2,
        "E(Carboxylation)"      : C + O*2,
        "K(Carboxylation)"      : C + O*2,
        "W(Carboxylation)"      : C + O*2,
        "N(Deamidation)"        : -1*H + -1*N + O,
        "Q(Deamidation)"        : -1*H + -1*N + O,
        "D(Dehydration)"        : -2*H + -1*O,
        "S(Dehydration)"        : -2*H + -1*O,
        "T(Dehydration)"        : -2*H + -1*O,
        "Y(Dehydration)"        : -2*H + -1*O,
        "K(Methylation)"        : C + H*2,
        "R(Methylation)"        : C + H*2,
        "M(Oxidation)"          : O,
        "D(SodiumAdduct)"       : -1*H + Na,
        "E(SodiumAdduct)"       : -1*H + Na,
        "K(12PlexDiLeu)"        : H*15 + C*7 + C13*1 + N15*1 + O18*1,
        "K(mdDiLeu1101)"        : H*15 + C*7 + C13*1 + N15*1 + O18*1,
        "K(mdDiLeu0400)"        : H*11 + C*8 + N*1 + O*1 + H2*4,
        "E(Glu->pyro-Glu)"      : -3*H + -1*N,
        "Q(Gln->pyro-Glu)"      : -3*H + -1*N,
        "Y(Sulfation)"          : O*3 + S,
        "S(Phosphorylation)"    : P + O*3 + H,
        "T(Phosphorylation)"    : P + O*3 + H,
        "Y(Phosphorylation)"    : P + O*3 + H,
    }
    
    replace_map = {
        "-Amidated": "(Amidated)",
        "12PlexDiLeu-": "(12PlexDiLeu)",
        "mdDiLeu1101-": "(mdDiLeu1101)",
        "mdDiLeu0400-": "(mdDiLeu0400)",
        "Acetylation-": "(Acetylation)",
        "Carbamidomethylation-": "(Carbamidomethylation)",
        "(Acetylation)K": "K(Acetylation)",
        "(Acetylation)S": "S(Acetylation)",
        "(Acetylation)T": "T(Acetylation)",
        "(Biotinylation)K": "K(Biotinylation)",
        "(Carbamidomethylation)C": "C(Carbamidomethylation)",
        "(Carbamidomethylation)D": "D(Carbamidomethylation)",
        "(Carbamidomethylation)E": "E(Carbamidomethylation)",
        "(Carbamidomethylation)H": "H(Carbamidomethylation)",
        "(Carbamidomethylation)K": "K(Carbamidomethylation)",
        "(Carboxylation)D": "D(Carboxylation)",
        "(Carboxylation)E": "E(Carboxylation)",
        "(Carboxylation)K": "K(Carboxylation)",
        "(Carboxylation)W": "W(Carboxylation)",
        "(Deamidation)N": "N(Deamidation)",
        "(Deamidation)Q": "Q(Deamidation)",
        "(Dehydration)D": "D(Dehydration)",
        "(Dehydration)S": "S(Dehydration)",
        "(Dehydration)T": "T(Dehydration)",
        "(Dehydration)Y": "Y(Dehydration)",
        "(Methylation)K": "K(Methylation)",
        "(Methylation)R": "R(Methylation)",
        "(Oxidation)M": "M(Oxidation)",
        "(SodiumAdduct)D": "D(SodiumAdduct)",
        "(SodiumAdduct)E": "E(SodiumAdduct)",
        "(Gln->pyro-Glu)Q": "Q(Gln->pyro-Glu)",
        "(Glu->pyro-Glu)E": "E(Glu->pyro-Glu)",
        "(Sulfation)Y": "Y(Sulfation)",
        "(12PlexDiLeu)K": "K(12PlexDiLeu)",
        "(mdDiLeu1101)K": "K(mdDiLeu1101)",
        "(mdDiLeu0400)K": "K(mdDiLeu0400)",
        "(Phosphorylation)S": "S(Phosphorylation)",
        "(Phosphorylation)T": "T(Phosphorylation)",
        "(Phosphorylation)Y": "Y(Phosphorylation)",
    }
    
    repl_pat = re.compile("|".join(map(re.escape, replace_map.keys())))
    def _normalize_tokens(seq: str) -> str:
        return repl_pat.sub(lambda m: replace_map[m.group(0)], seq)
    
    def n_term_delta(peptide: str) -> float:
        # N-term mod appears as "(Mod)..." at the start
        if peptide.startswith("("):
            term = peptide[: peptide.index(")") + 1]
            return n_term_mods.get(term, n_term_mods["Standard"])
        return n_term_mods["Standard"]

    def c_term_delta(peptide: str) -> float:
        # C-term mod appears as "...(Mod)" at the very end (and must be a known c-term mod)
        if peptide.endswith(")"):
            term = peptide[peptide.rfind("(") :]
            return c_term_mods.get(term, c_term_mods["Standard"])
        return c_term_mods["Standard"]

    _ptm_regex = re.compile(r"([A-Z])\(([^)]+)\)")

    
    def ptm_delta(peptide: str) -> float:
        # Sum residue-attached PTMs
        total = 0.0
        for res, mod in _ptm_regex.findall(peptide):
            total += PTMs.get(f"{res}({mod})", 0.0)
        return total

    def calc_mz_z1(peptide: str) -> float:
        # Neutral mass of plain residues
        plain = re.sub(r"[\(\[][^\)\]]*[\)\]]", "", peptide)  # strip (...) tokens
        try:
            base = sum(aa_masses[a] for a in plain)
        except KeyError as e:
            raise ValueError(f"Unknown residue in sequence '{peptide}': {e}")

        neutral = base + n_term_delta(peptide) + c_term_delta(peptide) + ptm_delta(peptide)
        # Report m/z for charge 1 (neutral + H+)
        return (neutral + proton_mass * charge) / charge
    
    # ---------------- Load FASTA ----------------
    
    seqs = fasta_to_df(fasta_path)
 
    # ---------------- Targets: variable-mod isoforms ----------------
    target_isoforms = []
    for seq in seqs:
        target_isoforms.extend(list(parser.isoforms(seq, variable_mods=variable_mod_dict, max_mods=max_mods_number)))

    # ---------------- Decoys ----------------
    if decoy_algorithm == "Shuffle":
        decoys = []
        for seq in seqs:
            s = "".join(random.sample(seq, len(seq)))
            # avoid trivial duplicates by re-shuffling once if needed
            if s in decoys:
                s = "".join(random.sample(seq, len(seq)))
            decoys.append(s)
    elif decoy_algorithm == "Reverse":
        decoys = [s[::-1] for s in seqs]
    else:
        raise ValueError(f"Unsupported decoy algorithm: {decoy_algorithm}")

    decoy_isoforms = []
    for seq in decoys:
        decoy_isoforms.extend(list(parser.isoforms(seq, variable_mods=variable_mod_dict, max_mods=max_mods_number)))

    all_peps = list(dict.fromkeys(target_isoforms + decoy_isoforms))

    norm_peps = [_normalize_tokens(p) for p in all_peps]

    # Compute m/z (z=1)
    masses = [calc_mz_z1(p) for p in norm_peps]
    
    df = pd.DataFrame({"Sequence": norm_peps,
        "Precursor theoretical monoisotopic mass": masses})

    df = df[~df["Sequence"].str.contains(r"\)\(")]

    # ---------------- Write & return ----------------
    os.makedirs(output_folder, exist_ok=True)
    out_path = os.path.join(output_folder, f"{decoy_algorithm}_database.csv")
    df.to_csv(out_path, index=False)
    return out_path
    
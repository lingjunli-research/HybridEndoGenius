# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 09:07:19 2025

@author: lafields2
"""

from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd

def fasta_to_df(fasta_path_get):
    """
    Read a FASTA file and collect raw sequences into a list.

    Parameters
    ----------
    fasta_path_get : str
        Path to the input FASTA file.

    Returns
    -------
    list[str]
        A list of sequences (strings) in the order they appear in the file.
        Headers are discarded.

    Notes
    -----
    Uses Biopython's ``SimpleFastaParser`` to iterate over records as
    ``(title, sequence)`` tuples and stores only the sequence string.
    This function does not validate characters, trim whitespace within
    sequences, or deduplicate entries.
    """
    fasta_to_df = []
    with open(fasta_path_get) as fasta_file:
        for title, sequence in SimpleFastaParser(fasta_file):
            fasta_to_df.append(sequence)
    return fasta_to_df
            
def raw_MS2_extraction(raw_file_formatted_path):
    raw_converter = pd.read_csv(
        raw_file_formatted_path,
        sep=",",
        skiprows=[0],
        names=["fragment_mz","fragment_resolution","fragment_z","fragment_intensity",
            "precursor_mz","ms2_scan","precursor_z","precursor_RT","IonInjectTime",
            "ms1_scan","precursor_intensity","null"]).drop_duplicates(subset="ms2_scan")
    return raw_converter
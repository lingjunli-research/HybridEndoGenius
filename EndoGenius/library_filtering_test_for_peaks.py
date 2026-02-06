# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:37:10 2024

@author: lafields2
"""

import pandas as pd
import csv

def library_filtering(lib_w_dups_path,out_path):
    
    def endogenius_score_filter(library,out):
        lib_sorted = library.sort_values(by='EndoGenius score', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\endogenius_score_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def norm_endogenius_score_filter(library,out):
        lib_sorted = library.sort_values(by='Norm. EndoGenius score', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\normalized_endogenius_score_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def hyperscore_filter(library,out):
        lib_sorted = library.sort_values(by='hyperscore', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\hyperscore_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
        return out_path_code
    
    def norm_hyperscore_filter(library,out):
        lib_sorted = library.sort_values(by='Norm. hyperscore', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\normalized_hyperscore_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def seq_coverage_filter(library,out):
        lib_sorted = library.sort_values(by='% sequence coverage', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\percent_sequence_coverage_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def tic_filter(library,out):
        lib_sorted = library.sort_values(by='precursor intensity', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\precursor_intensity_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def norm_tic_filter(library,out):
        lib_sorted = library.sort_values(by='Norm. precursor intensity', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\normalized_precursor_intensity_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def peaks_filter(library,out):
        lib_sorted = library.sort_values(by='fragment peak count', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\number_fragment_peaks_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def motif_filter(library,out):
        lib_sorted = library.sort_values(by='motif score', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\motif_score_peaks_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def avg_annot(library,out):
        lib_sorted = library.sort_values(by='avg # annotations per fragment', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\avg_annotations_per_fragment_ion_peaks_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)
    
    def avg_ions(library,out):
        lib_sorted = library.sort_values(by='avg # fragment ions per AA', ascending = False)
        lib_no_dups = lib_sorted.drop_duplicates(subset = ['Sequence'])
        out_path_code = out + '\\avg_fragments_per_AA_peaks_filter.csv'
        lib_no_dups.to_csv(out_path_code, index=False)    
      
    
    lib = pd.read_csv(lib_w_dups_path)
    # endogenius_score_filter(lib,out_path)
    # norm_endogenius_score_filter(lib,out_path)
    hyperscore_filter(lib,out_path)
    # norm_hyperscore_filter(lib,out_path)
    # seq_coverage_filter(lib,out_path)
    # tic_filter(lib,out_path)
    # norm_tic_filter(lib,out_path)
    # peaks_filter(lib,out_path)
    # motif_filter(lib,out_path)
    # avg_annot(lib,out_path)
    # avg_ions(lib,out_path)
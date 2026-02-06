# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:57:19 2024

@author: lafields2
"""
import re
import pandas as pd

def format_lib(working_directory):
    
    
    
    library_input = working_directory + '\\hyperscore_filter.csv'
    library_output = working_directory + '\\hyperscore_filter.tsv'

    
    library = pd.read_csv(library_input)
    
    seq_storage = []
    precursor_mz_storage = []
    precursor_z_storage = []
    rt_storage = []
    frag_loss_type_storage = []
    frag_type_storage = []
    fragment_series_number_storage = []
    fragment_mz_storage = []
    fragment_z_storage = []
    fragment_intensity_storage = []
    
    for row in range(0,len(library)):
        precursor_mz = library['precursor m/z'].iloc[row]
        precursor_z = library['precursor z'].iloc[row]
        rt = library['precursor RT'].iloc[row]
        seq = library['Sequence'].iloc[row]

        lib_entries = library['fragment peak list'].iloc[row]
        
        lib_entries = lib_entries.replace('[','')
        lib_entries = lib_entries.replace(']','')
        lib_entries = lib_entries.replace("'",'')
        
        lib_entries_list = lib_entries.split(',')
        
        for a in lib_entries_list:
            lib_details = a.split(':')
            
            fragment_type = lib_details[0]
            
            if 'NH3' in fragment_type:
                frag_loss_type = 'NH3'
            elif 'H2O' in fragment_type:
                frag_loss_type = 'H2O'
            else:
                frag_loss_type = 'None'
            
            
            if 'b' in fragment_type:
                fragment_type_by = 'b'
            elif 'y' in fragment_type:
                fragment_type_by = 'y'
            
            fragment_series_number = fragment_type.replace(fragment_type_by,'')
            fragment_series_number = fragment_series_number.replace(frag_loss_type,'')
            fragment_series_number = fragment_series_number.replace('-','')
            
            fragment_mz = lib_details[1]
            fragment_z = float(lib_details[2])
            fragment_int = lib_details[3]
            
            seq = seq.replace('(Amidated)','(UniMod:2)')
            seq = seq.replace('(Oxidation)','(UniMod:35)')
            seq = seq.replace('(Gln->pyro-Glu)','(UniMod:28)')
            seq = seq.replace('(Glu->pyro-Glu)','(UniMod:27)')
            seq = seq.replace('(Sulfation)','(UniMod:40)')
            seq = seq.replace('(Acetylation)','(UniMod:1)')
            seq = seq.replace('(Carbamidomethylation)','(UniMod:4)')
            seq = seq.replace('(Biotinylation)','(UniMod:3)')
            seq = seq.replace('(Carboxylation)','(UniMod:299)')
            seq = seq.replace('(Deamidation)','(UniMod:7)')
            seq = seq.replace('(Dehydration)','(UniMod:23)')
            seq = seq.replace('(Methylation)','(UniMod:34)')
            seq = seq.replace('(SodiumAdduct)','(UniMod:30)')
            seq = seq.replace('(Phosphorylation)','(UniMod:21)')
            seq = seq.replace('(12PlexDiLeu)','(UniMod:1321)')

            seq_storage.append(seq)
            precursor_mz_storage.append(precursor_mz)
            precursor_z_storage.append(precursor_z)
            rt_storage.append(rt)
            frag_loss_type_storage.append(frag_loss_type)
            frag_type_storage.append(fragment_type_by)
            fragment_series_number_storage.append(fragment_series_number)
            fragment_mz_storage.append(fragment_mz)
            fragment_z_storage.append(fragment_z)
            fragment_intensity_storage.append(fragment_int)
        
    library_formatted = pd.DataFrame({'PrecursorMz':precursor_mz_storage,
                                      'ProductMz':fragment_mz_storage,
                                      'Tr_recalibrated':rt_storage,
                                      'LibraryIntensity':fragment_intensity_storage,
                                      'ModifiedPeptide':seq_storage,
                                      'PrecursorCharge':precursor_z_storage,
                                      'FragmentSeriesNumber':fragment_series_number_storage,
                                      'FragmentType':frag_type_storage,
                                      'FragmentLossType':frag_loss_type_storage,
                                      'FragmentCharge':fragment_z_storage})
    
    def move_12plexdileu_to_first_res(peptide):
        return re.sub(r'^\(UniMod:1321\)([A-Z])', r'\1(UniMod:1321)', peptide)
    library_formatted["ModifiedPeptide"] = library_formatted["ModifiedPeptide"].apply(move_12plexdileu_to_first_res)
    library_formatted.to_csv(library_output, index=False)
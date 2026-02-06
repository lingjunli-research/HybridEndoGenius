# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:16:41 2024

@author: lafields2
"""

import pandas as pd
import os
import csv



def lib_w_dubs_generate(results_directory,output_directory):

    def get_file_names_with_strings_list(str_list, folder): #definition for finding a file containing a string in filename in specified directory
        full_list = os.listdir(folder)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return final_list
    
    
    csv_list = (get_file_names_with_strings_list(['.csv'],results_directory)) #search for the file based on query
    
    combined_library = pd.DataFrame()
    
    for csv in csv_list:
        intermediate_library_path = results_directory + '\\' + csv
        intermediate_library = pd.read_csv(intermediate_library_path)
        intermediate_library['Origin library'] = csv
        
        if len(combined_library) == 0:
            combined_library = intermediate_library
        else:
            combined_library = pd.concat([combined_library,intermediate_library])
    
    file_path = output_directory + '\\library_with_duplicates.csv'
    combined_library.to_csv(file_path, index=False)
    
    return file_path
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:56:44 2024

@author: lafields2
"""

from spec_lib_build_app_eg_score_v03 import generate_mini_libraries
from library_with_dups_generation import lib_w_dubs_generate
from library_filtering_test_for_peaks import library_filtering
from lib_formatting_for_DIANN import format_lib
import os

def build_a_SL(parent_results_directory,output_directory,fragment_error):

    def make_directory(directory):
        path = os.path.join(output_directory, directory)  
        if not os.path.exists(path):
            os.makedirs(path)
    
        else:
            pass
        return path
    
    first_lib_out_dir = str(make_directory("individual_libraries"))
    first_lib = generate_mini_libraries(parent_results_directory,first_lib_out_dir,fragment_error)
    
    dups_lib_out_dir = str(make_directory("intermediate_libraries"))
    lib_w_dups = lib_w_dubs_generate(first_lib_out_dir,dups_lib_out_dir)
    
    filtered_lib = library_filtering(lib_w_dups,dups_lib_out_dir)
    
    formatted_lib = format_lib(dups_lib_out_dir)

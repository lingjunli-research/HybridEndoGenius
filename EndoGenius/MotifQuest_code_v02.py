# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 20:37:53 2024

@author: 15156
"""
from Bio.Align.Applications import ClustalOmegaCommandline
import os
import subprocess
from Bio import AlignIO
import pandas as pd
import numpy as amp
from Bio import pairwise2
from Bio import SeqIO
import re
import math
import csv 
from scipy.cluster.hierarchy import linkage, fcluster, inconsistent
from scipy.spatial.distance import squareform

def start_building_a_motif_db(in_file,output_folder_path,t_val,min_motif_len,motif_instances,clustalo_path):

    flank_size = 0
    out_file = output_folder_path + "\\np.fasta"
    distance_matrix_file = "distance_matrixNPs.txt" 
    
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
    
    command_str_alignment = str(clustalomega_cline)
    
    
    working_dir = output_folder_path
    
    
    print("Starting alignment...")
    subprocess.call(command_str_alignment, cwd=working_dir, shell=True)
    
    print("Alignment finished.")
    unaligned_sequences = []
    for seq_record in SeqIO.parse(in_file, "fasta"):
        unaligned_sequences.append({
            "ID": seq_record.id,
            "Sequence": str(seq_record.seq)
        })
    
    
    unaligned_sequences_df = pd.DataFrame(unaligned_sequences)
    command_str_distance_matrix = f'{clustalo_path} -i "{in_file}" -o "{out_file}" --distmat-out="{distance_matrix_file}" --force --full -v'
    result = subprocess.run(command_str_distance_matrix, cwd=working_dir, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print("Distance matrix generated successfully.")
    test = (os.path.join(working_dir,in_file))
    
    print(test)
    database_length = sum(1 for _ in SeqIO.parse(in_file, "fasta"))

    if os.path.exists(os.path.join(working_dir, distance_matrix_file)):
        distance_matrix = amp.loadtxt(os.path.join(working_dir, distance_matrix_file), skiprows=1, usecols=range(1, database_length + 1))
        print("Distance matrix loaded successfully.")

    aligned_file = out_file
    
    # Parse the aligned sequences
    aligned_sequences = []
    for seq_record in SeqIO.parse(aligned_file, "fasta"):
        aligned_sequences.append({
            "ID": seq_record.id,
            "Sequence": str(seq_record.seq)
        })
    
    aligned_sequences_df = pd.DataFrame(aligned_sequences)

    sequence_labels = unaligned_sequences_df['Sequence'].tolist()
    
    distance_df = pd.DataFrame(distance_matrix, index=sequence_labels, columns=sequence_labels)
    
    condensed_dist_matrix = squareform(distance_df.values, checks=False)
    
    #hierarchical clustering
    Z = linkage(condensed_dist_matrix, 'average')
    
    for t in t_val:
        t_name = []
        num_clusters = []
        t_value_data = {}
        clusters_count = {} 
        
        
        def create_consensus(sequences, threshold):
            if not sequences: 
                return ''
            consensus = ['X'] * len(sequences[0])
            for i in range(len(sequences[0])):  
                position = [seq[i] for seq in sequences if i < len(seq)] 
                aa_count = {}
                for aa in position:
                    aa_count[aa] = aa_count.get(aa, 0) + 1
                most_common_aa, max_count = max(aa_count.items(), key=lambda x: x[1])
                freq = max_count / len(position)
                if freq >= threshold:
                    consensus[i] = most_common_aa
            return ''.join(consensus)
    
    
    
        clusters = fcluster(Z, t, criterion='distance')
        num_clusters = len(amp.unique(clusters))
        clusters_count[t] = num_clusters
        sequences_df = pd.DataFrame(sequence_labels, columns=['Sequence'])
        sequences_df['Cluster'] = clusters
        
        consensus_sequences = {}
        iamput_sequences_data = {}
        
        for cluster_id, group in sequences_df.groupby('Cluster'):
            cluster_sequences = group['Sequence'].tolist()
            if len(cluster_sequences) >= motif_instances: 
                consensus_sequence = create_consensus(cluster_sequences, threshold=1)
                consensus_sequences[cluster_id] = consensus_sequence
                iamput_sequences_data[cluster_id] = cluster_sequences
        
        t_value_data[t] = {
            "t_value": t,
            "number_of_clusters": num_clusters,
            "consensus_sequences": consensus_sequences,
            "iamput_sequences": iamput_sequences_data  
        }

        full_motif_storage = []
        partial_motif_storage = []
        
        
        for aa, aaa in t_value_data.items():
            if 'consensus_sequences' in aaa and isinstance(aaa['consensus_sequences'], dict):
                consensus_sequences = aaa['consensus_sequences']
                for seq_key, seq_value in consensus_sequences.items():
                    if isinstance(seq_value, str):
                        sequence_split = seq_value.split('-')
                        sequence_split = list(filter(None, sequence_split))
                    else:
                        pass
        
        
                    for seq in sequence_split:
                        if 'X' in seq:
                            x_ind = [n for (n, e) in enumerate(seq) if e == 'X'] #returns each index of X within the continuous sequence
                            num_index = len(x_ind)
                            for x in range(0,num_index-1): #purposely skipping last entry. Will check separately
                                if x_ind[x] > 0: #not doing a partial motif if the first value in the sequence is X
                                    current_value = x_ind[x]
                                    previous_value = x_ind[x-1]
                                    next_value = x_ind[x+1]
                
                                    pre_gap_size = current_value - previous_value
                                    post_gap_size = next_value - current_value
                                    if (pre_gap_size >= flank_size) and (post_gap_size >= flank_size): #makes sure flanking AAs check out
                                        partial_motif = seq[(previous_value+1):(next_value)] #extract partial motif
                                        partial_motif_storage.append(partial_motif)
                                        partial_motif_storage = list(set(partial_motif_storage))
                                    else:
                                        pass
                                else:
                                    pass
                            seq_len = len(seq)
                            last_ind = x_ind[-1]
                            if (seq_len - last_ind) >= flank_size: #check if a partial motif applies to the last instance of X in the sequence
                                 if len(x_ind) > 1:
                                     previous_index = x_ind[-2]
                                     last_flank = last_ind - previous_index
                                     if last_flank >= flank_size:
                                         last_partial = seq[(previous_index+1):(seq_len)]
                                         partial_motif_storage.append(last_partial)
                                     else:
                                         pass
                                 else:
                                     partial_motif_storage.append(seq)
                            elif seq[-1] == 'X':
                                seq_shorten = seq[0:-1]
                                if 'X' not in seq_shorten:
                                    if len(seq_shorten) >= min_motif_len:
                                        full_motif_storage.append(seq_shorten)
                                else:
                                    seq_split = seq_shorten.split('X')
                                    for seq_splitted in seq_split:
                                        if len(seq_splitted) >= min_motif_len:
                                            full_motif_storage.append(seq_splitted)
                                        else:
                                            pass
                            elif seq[0] == 'X':
                                seq_shorten = seq[1:]
                                if 'X' not in seq_shorten:
                                    if len(seq_shorten) >= min_motif_len:
                                        full_motif_storage.append(seq_shorten)
                                else:
                                    seq_split = seq_shorten.split('X')
                                    for seq_splitted in seq_split:
                                        if len(seq_splitted) >= min_motif_len:
                                            full_motif_storage.append(seq_splitted)
                                        else:
                                            pass
                            else:
                                 pass
                
                        
                        elif len(seq) >= min_motif_len: #store full motifs (without X)
                            full_motif_storage.append(seq)
                        else:
                            pass
        
        
        number_of_motifs = len(set(full_motif_storage + partial_motif_storage))
 
        full_motif_storage = list(set(full_motif_storage))
        partial_motif_storage = list(set(partial_motif_storage)) #shouldnt have any but good measure 
        
    
        full_motif_matches = {}
        
        for motif in full_motif_storage:
            matching_sequences_set = set()
        
            for sequence in unaligned_sequences_df['Sequence']:
                if motif in sequence:
                    matching_sequences_set.add(sequence)
        
            full_motif_matches[motif] = list(matching_sequences_set)
           
        
        def motif_to_regex(par_mot):
            return par_mot.replace('X', '.')
        
        partial_matching_motifs = {}
        
        for par_mot in partial_motif_storage:
            regex_pattern = motif_to_regex(par_mot)
            
            partial_matching_sequences = []
            
            for sequence in unaligned_sequences_df['Sequence']:
                if re.search(regex_pattern, sequence):
                    partial_matching_sequences.append(sequence)
            
            partial_matching_motifs[par_mot] = partial_matching_sequences
        print("matching partial motifs finished.")

        def EG_motif_dictionaries(partial_matching_motifs, full_motif_matches):
            results = {
                "partial_motif_EG_score": {},
                "full_motif_EG_score": {}
            }
        
            def EG_calcuation(d, i):
                for motif, sequence in d.items():
                    EG_results = []
                    for x in sequence:
                        EG_calc = (len(motif) / len(x)) * math.sqrt(len(x))
                        EG_results.append((EG_calc, x))
                    results[i][motif] = EG_results
        
            EG_calcuation(partial_matching_motifs, "partial_motif_EG_score")
            EG_calcuation(full_motif_matches, "full_motif_EG_score")
        
            return results
        
        EG_score_storage = EG_motif_dictionaries(partial_matching_motifs, full_motif_matches)
        print("EG score finished.")
        ##########################################
        #final motif score
        
        def find_t_value_for_motif(motif, t_value_data):
            for t_data in t_value_data.values():
                consensus_sequences = t_data.get("consensus_sequences", {})
                if motif in consensus_sequences:
                    return t_data["t_value"]
                for seq_key, seq_value in consensus_sequences.items():
                    if motif in seq_value:
                        return t_data["t_value"]
            return None
        
        def final_motif_score(EG_score_storage, database_length, t_value_data):
            final_scores = {
                "partial_motif_final_score": {},
                "full_motif_final_score": {}
            }
        
            def avg_first_elements(tuples_list):
                try:
                    avg = sum(x[0] for x in tuples_list) / len(tuples_list)
                except ZeroDivisionError:
                    avg = 0
                return avg
        
            for score_type, motifs in EG_score_storage.items():
                final_score_type = "partial_motif_final_score" if score_type == "partial_motif_EG_score" else "full_motif_final_score"
                total_size = sum(len(v) for v in motifs.values())
                factor = total_size / math.log2(database_length)
        
                for key, values in motifs.items():
                    avg_floats = avg_first_elements(values)
                    final_score = avg_floats * factor
                    t_value = find_t_value_for_motif(key, t_value_data)
                    matching_sequences = [x[1] for x in values]
                    final_scores[final_score_type][key] = {
                        "score": final_score,
                        "sequences": matching_sequences,
                        "t_value": t_value
                    }
        
            return final_scores
        print("final score finished.")
        final_scores = final_motif_score(EG_score_storage, database_length, t_value_data)
        
        full_motif_df = pd.DataFrame.from_dict(final_scores['full_motif_final_score'], orient='index')
        full_motif_df.index.name = 'motif'  
        
        partial_motif_df = pd.DataFrame.from_dict(final_scores['partial_motif_final_score'], orient='index')
        partial_motif_df.index.name = 'motif'  
        
        full_motif_df.to_csv(output_folder_path + '\\full_motif_report_' + str(t) +'.csv')
        partial_motif_df.to_csv(output_folder_path + '\\partial_motif_report_' + str(t) + '.csv')

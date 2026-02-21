# import libraries
from __future__ import print_function


import os
import numpy as np
import pandas as pd
from multiprocessing import Pool


from rapidfuzz import fuzz, process
from Bio import SeqIO
from Bio.Seq import Seq


# functions
def extract_psm_to_csv(mztab_file, output_csv):
    with open(mztab_file, 'r') as f:
        lines = f.readlines()

    # Extract the PSH header line
    header_line = next(line for line in lines if line.startswith("PSH"))
    headers = header_line.strip().split('\t')[1:]  # drop the 'PSH' prefix

    # Extract PSM lines
    psm_lines = [line.strip().split('\t')[1:] for line in lines if line.startswith("PSM")]

    # Convert to pandas DataFrame
    import pandas as pd
    df = pd.DataFrame(psm_lines, columns=headers)

    # Save to CSV
    df.to_csv(output_csv, index=False)

def read_denovo(denovo_csv, selected_features=None):
    denovo_psm = pd.read_csv(denovo_csv, keep_default_na=False)
    denovo_psm['Scan'] = denovo_psm['spectra_ref'].apply(lambda x: x.split('scan=')[-1])
    denovo_psm['feature_id'] = denovo_psm.apply(
            lambda row: row['spectra_ref'].split(':scan=')[0] + '||' + str(row['Scan']), axis=1
        )
    if selected_features:
        denovo_psm['selected_features'] = denovo_psm['feature_id'].isin(selected_features)
        denovo_psm = denovo_psm[denovo_psm['selected_features']]
    return denovo_psm

def calculate_FDR(target_csv, decoy_csv, engine_score, fdr_list, selected_features=None):
    
    print("target_csv =", target_csv)
    print("decoy_csv =", decoy_csv)
    target_psm = read_denovo(target_csv, selected_features)
    decoy_psm = read_denovo(decoy_csv, selected_features)
    print("len(target_psm) =", len(target_psm)); print("len(decoy_psm) =", len(decoy_psm))
    dfs = pd.concat([target_psm, decoy_psm], keys=['target', 'decoy']).reset_index().rename(columns={'level_0': 'spectrum'})
    dfs['is_target'] = dfs.apply(lambda row: row['spectrum']=='target', axis=1)

    # target-decoy competition
    dfs.sort_values(by=[engine_score, 'is_target'], ascending=[False, False], inplace=True)
    # 1-1 competition on each scan id
#     dfs_fdr = dfs.drop_duplicates(subset=['feature_id'])
    # competition on whole dataset
    dfs_fdr = dfs
    dfs_fdr['feature_id'] = dfs_fdr.apply(lambda x: x['feature_id'] if x['is_target'] else x['feature_id']+'||decoy', axis=1)
    print("len(dfs) =", len(dfs))
    print("len(dfs_fdr) =", len(dfs_fdr))
    print("sum(dfs_fdr['is_target']) =", sum(dfs_fdr['is_target']))
    
    # fdr estimation
    cumsum = range(1, len(dfs_fdr) + 1)
    cumsum_target = np.cumsum(np.array(dfs_fdr['is_target'].astype(int)))
    cumsum_decoy = cumsum - cumsum_target
    estimated_fdr = cumsum_decoy / cumsum_target
    dfs_fdr['estimated_fdr'] = estimated_fdr

    score_list = []
    count_list = []
    for fdr in fdr_list:
        fdr_index = np.flatnonzero(estimated_fdr <= fdr)
        fdr_index = fdr_index[-1] if len(fdr_index) > 0 else 0
        score_list.append(dfs_fdr.iloc[fdr_index][engine_score])
        count_list.append(fdr_index + 1)

    return dfs, dfs_fdr, score_list, count_list

def casanovo_filter(original_csv, output_csv, fdr_value):
    df = pd.read_csv(original_csv)
    filtered_df = df[df['search_engine_score[1]'] > fdr_value].copy()
    filtered_df['sequence_new'] = filtered_df['sequence'].str.replace(r"[^A-Z]","",regex=True)
    filtered_df.to_csv(output_csv, index=False)

def database_compare(filtered_csv, database_fasta, output_csv):
    df_original = pd.read_csv(filtered_csv)
    denovo_peptides = df_original['sequence_new'].tolist()

    ##### 1. Search known neuropeptides
    fasta_peptides = []
    for record in SeqIO.parse(database_fasta, "fasta"):
        seq = str(record.seq).strip()
        if seq:
            fasta_peptides.append(seq)
        
    results = []
    for query in denovo_peptides:
        match, score, _ = process.extractOne(query, fasta_peptides, scorer=fuzz.ratio)
        results.append((query, match, score))

    # === Save results to CSV ===   ####################################################################################
    df = pd.DataFrame(results, columns=["DeNovo_Peptide", "Matched_Database_Seq", "Similarity_Score"])
    df.to_csv(output_csv, index=False)

def match_update(db_csv,target_csv,output_dir):
    fasta_output = os.path.join(output_dir, "novel_peptides.fasta")
    dfs = []
    if db_csv:
        all_db = pd.read_csv(db_csv)
    else:
        all_db = pd.DataFrame()
    df_filter = pd.read_csv(target_csv)
    df_seq = df_filter[['sequence_new']]
    to_remove = set(all_db.loc[all_db["Similarity_Score"] == 100, "DeNovo_Peptide"])
    df_seq = df_seq[~df_seq['sequence_new'].isin(to_remove)]
    df_seq = df_seq.drop_duplicates(subset=['sequence_new'], keep='first')
    df_seq = df_seq.reset_index(drop=True)
    df_seq.index += 1
    with open(fasta_output, "w") as f:
        for idx, row in df_seq.iterrows():
            f.write(f">ID_{idx}\n{row['sequence_new']}\n")
    return fasta_output

def automated_data_analysis(target_casanovo_result, decoy_casanovo_result, engine_score_column, required_fdr, output_dir, db):
    # file name organize
    target_base, _ = os.path.splitext(target_casanovo_result)
    decoy_base, _ = os.path.splitext(decoy_casanovo_result)
    target_csv_path = target_base + ".csv"
    decoy_csv_path = decoy_base + ".csv"
    filtered_target_path = target_base + "_filtered.csv"
    extract_psm_to_csv(target_casanovo_result, target_csv_path)
    extract_psm_to_csv(decoy_casanovo_result, decoy_csv_path)
    # decoy FDR calculation
    dfs, dfs_fdr, score_list, count_list = calculate_FDR(target_csv = target_csv_path, 
                                                         decoy_csv= decoy_csv_path, 
                                                         engine_score = engine_score_column,
                                                         fdr_list=[float(required_fdr)]
                                                         )
    # filter casanovo results with fdr cutoff
    casanovo_filter(target_csv_path, filtered_target_path, score_list[0])
    # compare with database
    output_csv = os.path.join(output_dir, "similarity_check.csv")
    database_compare(filtered_target_path, db, output_csv)
    # database_compare(filtered_target_path, db3, db3_path)
    # filter again to get only novel neuropeptides
    match_update(output_csv,filtered_target_path,output_dir)
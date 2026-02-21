#!/bin/bash

# help function
help_message() {
    echo "HELP PAGE..."
    echo
    echo
    echo "Syntax: bash make_dag.sh -s <sample_basename> -j <round1_json> -k <round2_json> -a <endogenius_fasta> -m <motif_csv> -f <fdr> -d <database_fasta> -o <output_filename>"
    echo "These arguments can be provided in any order, but all arguments are required. Options are case-sensitive."
    echo
    echo "Options:"
    echo "  -s    (required) Sample basename: Base name of your sample files (without extension). Example: sample"
    echo "                   The script will automatically look for sample.ms2, sample.mgf, and sample.mzML"
    echo "  -j    (required) Round1 JSON: Configuration file for EndoGenius round 1. Example: round1.json"
    echo "  -k    (required) Round2 JSON: Configuration file for EndoGenius round 2. Example: round2.json"
    echo "  -a    (required) EndoGenius FASTA database. Example: EG_database.fasta"
    echo "  -m    (required) EndoGenius Motif database (.csv). Example: motif.csv"
    echo "  -f    (required) FDR threshold. Example: 0.05"
    echo "  -d    (required) Database FASTA file for Casanovo result similarity comparison. Example: filter.fasta"
    echo "  -o    (required) DAG output file name: Desired name for DAG file. Example: ms_pipeline"
    echo
    echo "Note: Please convert your MS2 file to MGF and mzML formats using msconvert before running this pipeline."
    echo "      The script expects all three files (sample.ms2, sample.mgf, sample.mzML) to exist."
    echo
    echo "Example usage: bash make_dag.sh -s sample -j round1.json -k round2.json -f 0.05 -d filter.fasta -o ms_pipeline"
    exit 1
}

# Initialize variables
SAMPLE_BASE=""
ROUND1_JSON=""
ROUND2_JSON=""
ENDOGENIUS_FASTA=""
MOTIF_CSV=""
FDR_VALUE=""
DATABASE_FASTA=""
FILENAME=""

# Fixed parameters (hardcoded)
CASANOVO_YAML="casanovo.yaml"
CASANOVO_CKPT="casanovo_v4_2_0.ckpt"

# Parse flags
while getopts "s:j:k:a:m:f:d:o:h" flag; do
    case "${flag}" in
        s) SAMPLE_BASE=${OPTARG} ;;
        j) ROUND1_JSON=${OPTARG} ;;
        k) ROUND2_JSON=${OPTARG} ;;
	a) ENDOGENIUS_FASTA=${OPTARG} ;;
	m) MOTIF_CSV=${OPTARG} ;;
	f) FDR_VALUE=${OPTARG} ;;
	d) DATABASE_FASTA=${OPTARG} ;;
        o) FILENAME=${OPTARG} ;;
        h) help_message ;;
        *) help_message ;;
    esac
done

# Check that all required inputs are provided
if [[ -z "$SAMPLE_BASE" || -z "$ROUND1_JSON" || -z "$ROUND2_JSON" || -z "$ENDOGENIUS_FASTA" || -z "$MOTIF_CSV" || -z "$FDR_VALUE" || -z "$DATABASE_FASTA" || -z "$FILENAME" ]]; then
    echo "Please check that you have provided all required inputs."
    help_message
    exit 1
fi

# Construct file names from sample basename
MS2_FILE="${SAMPLE_BASE}.ms2"
MGF_FILE="${SAMPLE_BASE}.mgf"
MZML_FILE="${SAMPLE_BASE}.mzML"

# Derive output file names
DECOY_MGF="${SAMPLE_BASE}_decoy.mgf"
TARGET_MZTAB="${SAMPLE_BASE}.mztab"
DECOY_MZTAB="${SAMPLE_BASE}_decoy.mztab"
FILTERED_FASTA="results/novel_peptides.fasta"

# Create logs directory if it doesn't exist
mkdir -p logs

# Main script
rm -f "$FILENAME.dag"
echo "Creating DAG named $FILENAME.dag"
touch "$FILENAME.dag"

# Write job definitions (msconvert removed - users convert manually)
echo "JOB ENDOGENIUS_R1 endogenius_round1.sub" >> "$FILENAME.dag"
echo "VARS ENDOGENIUS_R1 json_config=\"$ROUND1_JSON\" mzml_file=\"$MZML_FILE\" ms2_file=\"$MS2_FILE\" fasta_db=\"$ENDOGENIUS_FASTA\" motif_db=\"$MOTIF_CSV\"" >> "$FILENAME.dag"

echo "JOB GENERATE_DECOY generate_decoy.sub" >> "$FILENAME.dag"
echo "VARS GENERATE_DECOY mgf_file=\"$MGF_FILE\"" >> "$FILENAME.dag"

echo "JOB CASANOVO_TARGET casanovo_target.sub" >> "$FILENAME.dag"
echo "VARS CASANOVO_TARGET mgf_file=\"$MGF_FILE\" yaml_file=\"$CASANOVO_YAML\" ckpt_file=\"$CASANOVO_CKPT\"" >> "$FILENAME.dag"

echo "JOB CASANOVO_DECOY casanovo_decoy.sub" >> "$FILENAME.dag"
echo "VARS CASANOVO_DECOY mgf_file=\"$DECOY_MGF\" yaml_file=\"$CASANOVO_YAML\" ckpt_file=\"$CASANOVO_CKPT\"" >> "$FILENAME.dag"

echo "JOB ANALYZE_FDR filter_casanovo.sub" >> "$FILENAME.dag"
echo "VARS ANALYZE_FDR target_mztab=\"$TARGET_MZTAB\" decoy_mztab=\"$DECOY_MZTAB\" fdr=\"$FDR_VALUE\" database=\"$DATABASE_FASTA\"" >> "$FILENAME.dag"

echo "JOB ENDOGENIUS_R2 endogenius_round2.sub" >> "$FILENAME.dag"
echo "VARS ENDOGENIUS_R2 json_config=\"$ROUND2_JSON\" mzml_file=\"$MZML_FILE\" ms2_file=\"$MS2_FILE\" filtered_fasta=\"$FILTERED_FASTA\" motif_db=\"$MOTIF_DB\"" >> "$FILENAME.dag"

# Write dependencies (msconvert removed)
echo "PARENT ENDOGENIUS_R1 CHILD CASANOVO_TARGET" >> "$FILENAME.dag"
echo "PARENT GENERATE_DECOY CHILD CASANOVO_DECOY" >> "$FILENAME.dag"
echo "PARENT CASANOVO_TARGET CASANOVO_DECOY CHILD ANALYZE_FDR" >> "$FILENAME.dag"
echo "PARENT ANALYZE_FDR CHILD ENDOGENIUS_R2" >> "$FILENAME.dag"

echo "Created $FILENAME.dag"

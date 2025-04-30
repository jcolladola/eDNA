#!/bin/bash
#SBATCH --job-name=blast_array
#SBATCH --output=/home/meg/rgallego/reports/blast_%A_%a.out
#SBATCH --error=/home/meg/rgallego/reports/blast_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

# First, outside this script,  you have to split the huge input files with
# awk '/^>/{n++} {print > sprintf("split_part_%04d.fasta", int(n/4000))}' input.fasta
# Second, count the number of input files with ls split_part_* | wc -l
# usage  sbatch --array=1-{numfiles}%6 blast_array.sh [splitdir] [OUTPUTDIR] [db]

# Set variables
INPUT_DIR=$1
DB=$3
OUTPUT_DIR=$2
SPLIT_FILES=($(ls $INPUT_DIR/split_part_*))  # Array of input files
QUERY_FILE=${SPLIT_FILES[$SLURM_ARRAY_TASK_ID - 1]}  # Select file based on job ID

# Run BLAST
blastn -query "$QUERY_FILE" -db "$DB" -num_threads $SLURM_CPUS_ON_NODE -outfmt "6 std staxids qlen" -out "$OUTPUT_DIR/$(basename $QUERY_FILE).out"

